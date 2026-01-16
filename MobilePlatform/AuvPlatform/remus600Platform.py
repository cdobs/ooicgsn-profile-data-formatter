"""
class: remus600Platform

description: Driver class for translating Remus 600 Subset data file data
into NetCDF files for import into GliderDac repository

history:
09/21/2021 ppw created
01/16/2026 cdobs - adds cdom correction factor
"""
from MobilePlatform.AuvPlatform.auvPlatform import auvPlatform
from FileReader.jsonCfgReader import jsonCfgReader
from FileReader.AuvReader.remus600SubsetDataReader import remus600SubsetDataReader
from DataProcessor.AuvProcessor.remus600Processor import remus600Processor
from FileWriter.NetCDFWriter.dacNetCDFWriter import dacNetCDFWriter
from FileWriter.NetCDFWriter.dataExplorerNetCDFWriter import dataExplorerNetCDFWriter
import common.constants as cc
import tempfile
import FileReader.AuvReader.remus600SubsetData as r600data
import FileReader.AuvReader.remus600SubsetMsgData as r600msgdata
import os
import logging
import json
import datetime
import pandas
import numpy as np

class remus600Platform( auvPlatform ) :

    def __init__( self ) :
        super().__init__()

        self.cfgReader = jsonCfgReader()
        self.dataFileReader = remus600SubsetDataReader()
        self.dataProcessor = remus600Processor()
        self.outputFileWriter = dacNetCDFWriter()

        # ** extract platform specific args here **

        # Any overrides of default file readers/writers/processors goes here


    def isValidFile(thePath, theFile ):
        """
        Verifies validity and existence of theFile at thePath
        :param thePath
        :param theFile:
        :return: True or False
        """

        ret = True
        if not os.path.isfile( os.path.join( thePath, theFile) ):
            logging.error(
                'Expected file {:s} not found at {:s}'.format(
                    theFile, thePath ) )
            ret = False
        return ret


    def validateSettings( self ):
        """
        validate existence of config files
        REMUS 600 requires deployment, global_attributes,
        instruments and sensor_defs as json format config files
        :return: 0 - valid, -1 - invalid
        """

        ret = remus600Platform.isValidFile( self.cfgPath, 'deployment.json' )
        ret = ret and remus600Platform.isValidFile( self.cfgPath, 'global_attributes.json' )
        ret = ret and remus600Platform.isValidFile(self.cfgPath, 'instruments.json')
        ret = ret and remus600Platform.isValidFile(self.cfgPath, 'sensor_defs.json')

        # caller expects 0 for valid, -1 invalid
        if ret:
            return 0
        else:
            logging.error('Missing required configuration file(s)')
            return -1


    def getInstrumentFromCfg(instrumentsCfg, instrumentName):
        """
        Search instrumentCfg dictionary for instrument with field
        nc_var_name set to instrumentName
        :param instrumentCfg
        :param instrumentName:
        :return: instr (dictionary), else None
        """

        for instr in instrumentsCfg:
            if 'nc_var_name' in instr:
                if instr['nc_var_name'] == instrumentName:
                    return instr

        logging.warning( 'No configuration found for instrument ' + instrumentName )
        return None

    def getSensorDefFromCfg(sensorCfg, sensorName):
        """
        Finds sensor in sensorCfg (dictionary) having
        nc_var_name set to sensorName
        :param sensorCfg
        :param sensorName:
        :return: sensor (dictionary), else None
        """

        for sensor, sensorAttrs in sensorCfg.items():
            if 'nc_var_name' in sensorAttrs:
                if sensorAttrs['nc_var_name'] == sensorName:
                    return sensorAttrs

        logging.warning( 'No configuration found for sensor ' + sensorName )
        return None

    def adjust_ctd_msg_id( self, data ) :

        # Some old deployments used an old CTD with a different msg_id.
        # If no data returned, try the alternate message id
        
        ctdCfg = remus600Platform.getInstrumentFromCfg( self.instrumentsCfg, 'instrument_ctd')
        ctdData = data.getDataForMessageId( int(ctdCfg['attrs']['subset_msg_id']) )
        if ctdData is None:
            OLD_CTD_MSG_ID = 1181
            if data.getDataForMessageId( OLD_CTD_MSG_ID ) is not None:
                ctdCfg['attrs']['subset_msg_id'] = OLD_CTD_MSG_ID
            else :
                logging.error('FATAL: No CTD message type found in input file')
                sys.exit()
            

    def useCtdDataToComputeProfiles(self, data):
        """
        Wrap profile computation to force hi-res ctd data
        to go out of scope when done
        :param data: remus subset file data
        :return: list of profile start, end times
        """
        ctdCfg = remus600Platform.getInstrumentFromCfg(
            self.instrumentsCfg, 'instrument_ctd')
        ctdData = data.getDataForMessageId( int(ctdCfg['attrs']['subset_msg_id']) )

        allProfileBounds = self.dataProcessor.computeProfiles(
            ctdData, "timestamp", "missionTime", "depth")

        if allProfileBounds is None or len(allProfileBounds) == 0:
            logging.warning('No valid profiles found in data')

        return allProfileBounds

    def sensorAttrMatches(sensorDef, attrName, match):
        """
        Determines if sensorDef dictionary has an attribute named 'attrName'
        with a value that matches 'match'
        :param attrName:
        :param match:
        :return: True - has matching named attribute, False - does not
        """

        if 'attrs' in sensorDef:
            if attrName in sensorDef['attrs']:
                if sensorDef['attrs'][attrName] == match:
                    return True
        return False

    def sensorHasAttr(sensorDef, attrName):
        """
        Determines if sensorDef has an attribute with the passed name
        :param sensorDef
        :param attrName:
        :return: True - has named attribute, False - does not
        """

        if 'attrs' in sensorDef:
            return attrName in sensorDef['attrs']


    def isInstrument(instrCfg, name ):
        """
        Determines if instrCfg list has an item with a key of
        'nc_var_name' with a value matching 'name'
        :param instrCfg
        :param name:
        :return: True - match, False - no match
        """

        for instr in instrCfg:
            if 'nc_var_name' in instr and instr['nc_var_name'] == name:
                return True

        return False


    def getDataSlice( dataset, startTime, endTime ):
        """
        Slices a Pandas DataFrame having a 'timestamp' column to
        within the start and end times
        :param dataset:
        :param startTime:
        :param endTime:
        :return: Pandas DataFrame slice w/i time bounds
        """

        # return data within passed time window

        sliceRange = dataset['timestamp'].between( startTime, endTime )

        return dataset[ sliceRange ]


    def setupFormatting(self):
        """
        Virtual method for performing setup work prior to formatting data for output
        :return: None
        """

        # read configuration settings

        self.globalsCfg = self.readCfgFile( self.cfgPath, 'global_attributes.json')
        self.deploymentCfg = self.readCfgFile( self.cfgPath, 'deployment.json')
        self.instrumentsCfg = self.readCfgFile( self.cfgPath, 'instruments.json')
        self.sensorsCfg = self.readCfgFile( self.cfgPath, 'sensor_defs.json')

        # pass settings to output file writer

        self.outputFileWriter.outputPath = self.outputPath
        self.outputFileWriter.overwriteExistingFiles = self.replaceOutputFiles
        self.outputFileWriter.outputCompressionLevel = self.outputCompression
        self.outputFileWriter.writeFormat = self.outputFormat

    def FormatData(self ):
        """
        For each data file passed, use the configuration settings to drive the
        generation of one or more output files (1 / profile).
        :return: 0
        """

        ret = 0

        # If debug mode, go no further
        if self.suppressOutput:
            logging.info('Output suppression indicated, terminating processing')
            return 0

        # for each data file

        for dataFile in self.dataFiles:

            # Establish temp directory for splitting Remus
            # subset files into individual messages (instruments)

            logging.debug( 'Processing ' + dataFile )

            outputFiles = []

            with tempfile.TemporaryDirectory() as tempPath:

                # read in the subset data file

                data = self.dataFileReader.read( dataFile, tempPath )
                if data is None:
                    logging.error("Bad data file encountered {:s}".format(dataFile))
                    ret = -1
                    continue

                # If necessary, remap the CTD subset_message_id
                self.adjust_ctd_msg_id( data )

                # Compute 1 sec res. gps data once for whole data file
                # Gps data is 1 second cadence at surface, gaps during dives.

                gpsCfg = remus600Platform.getInstrumentFromCfg(
                    self.instrumentsCfg, 'instrument_gps' )
                gpsData = data.getDataForMessageId( int( gpsCfg['attrs']['subset_msg_id'] ))
                gpsDataNoGaps = self.dataProcessor.interpolateGpsData( gpsData )

                #dumpfile = open( '/tmp/gps_04052019.csv', 'w')
                #dumpfile.write( 'timestamp,lat,lon\n')
                #for iii in range( len(gpsDataNoGaps.timestamp) ):
                #    dumpfile.write(str(gpsDataNoGaps.timestamp[iii]) + ',' +
                #                   str(gpsDataNoGaps.latitude[iii]) + ',' +
                #                   str(gpsDataNoGaps.longitude[iii]) + '\n')
                #dumpfile.close()

                # compute profile bounds using data from CTD

                allProfileBounds = self.useCtdDataToComputeProfiles( data )
                if allProfileBounds is None or len(allProfileBounds) == 0:
                    logging.warning('No valid profiles found in data file, skipping.')
                    ret = -1
                    continue

                # for each profile (id unique w/i trajectory 1..n)

                profileId = 1
                for profileBounds in allProfileBounds:

                    logging.debug('Processing profile ' + str( profileId ))

                    # Each profile requires a separate output file for GliderDac
                    # Re-init the output writer for each profile

                    filename = self.deploymentCfg['glider'] + "_" + \
                               datetime.datetime.fromtimestamp(
                                   profileBounds[0] ).strftime('%Y%m%dT%H%M') + "_" + \
                               self.deploymentCfg['global_attributes']['mode'] + ".nc"
                    outputFiles.append( filename )

                    self.outputFileWriter.resetAll()
                    self.outputFileWriter.fileName = filename
                    self.outputFileWriter.profileId = profileId
                    self.outputFileWriter.profileStartTime = profileBounds[0]
                    self.outputFileWriter.profileEndTime = profileBounds[1]
                    self.outputFileWriter.trajectory = \
                        self.deploymentCfg['trajectory_name']
                    self.outputFileWriter.trajectoryDatetime = \
                        self.deploymentCfg['trajectory_datetime']
                    self.outputFileWriter.sourceFile = dataFile

                    try:
                        self.formatProfileData( profileId, profileBounds[0], profileBounds[-1],
                                                allProfileBounds, data, gpsDataNoGaps )

                        # Generate an output file

                        self.outputFileWriter.setupOutput()
                        self.outputFileWriter.writeOutput()
                        self.outputFileWriter.cleanupOutput()

                    except Exception as e:
                        logging.warning( "Profile " + str(profileId) + " invalid, ignored ")

                    profileId = profileId + 1


            # If output target is OOI Explorer, feed the output
            # files formatted for GliderDAC to the OOI Explorer
            # file writer for reformatting.

            if self.targetHost == cc.OOI_EXPLORER_TARGET:
                if len(outputFiles) > 0:
                    deWriter = dataExplorerNetCDFWriter()
                    deWriter.outputPath = self.outputPath
                    deWriter.overwriteExistingFiles = self.replaceOutputFiles
                    deWriter.outputCompressionLevel = self.outputCompression
                    deWriter.writeFormat = self.outputFormat
                    deWriter.deploymentId = 'R' + \
                       self.deploymentCfg['global_attributes']['deployment_number']
                    deWriter.trajectoryName = self.deploymentCfg['trajectory_name']
                    deWriter.trajectoryDateTime = self.deploymentCfg['trajectory_datetime']
                    deWriter.sourceFile = dataFile
                    deWriter.inputFiles = outputFiles

                    deWriter.setupOutput()
                    deWriter.writeOutput()
                    deWriter.cleanupOutput()
                else:
                    logging.warning('No output NetCDF files produced, conversion to OOI format skipped.')
                    ret = -1

        return ret

    def cleanupFormatting(self):
        """
        Virtual method for performing post data formatting cleanup activities
        :return: 0
        """

        return 0

    def formatProfileData( self, profileId, profileStartTime,
                           profileEndTime, allProfileBounds, data, gpsData ):
        """
        Generate the output attributes and variables, then write output file for one profile
        :param profileId:
        :param profileStartTime:
        :param profileEndTime:
        :param allProfileBounds:
        :param data:
        :param gpsData:
        :return: 0
        """

        # Add time invariant attributes and vars

        self.formatGlobalAttributes( profileId )
        self.formatInstrumentVars()

        # Perform any custom variable calculations before creating profile variables

        calculatedVars = self.computeCalculatedVars(
            profileId, profileStartTime, profileEndTime,
            allProfileBounds, data, gpsData )

        # Create profile variables

        for sensorName, sensorDef in self.sensorsCfg.items():

            dimension = sensorDef.get('dimension')
            if dimension is not None and dimension == "time":

                # get data for sensor's instrument within profile bounds

                profileData = self.getProfileData( sensorDef, data, gpsData,
                                                   profileStartTime, profileEndTime ).copy()

                # combine time fields to get time at finest available resolution

                dataTimesMs = data.timesInMillisecs( profileData.get('timestamp'),
                                                     profileData.get('missionTime') )

                # Purge out of range values, based on sensor configured
                # valid_min, valid_max and __FillValue settings

                self.purgeOutOfRangeData( profileId, profileData, sensorDef, dataTimesMs )

            # use profile data directly for measured sensors

            if remus600Platform.sensorAttrMatches( sensorDef, 'observation_type', 'measured'):

                self.formatMeasuredVar( sensorDef, profileData, dataTimesMs )

            # use calculated data for calculated sensors

            elif remus600Platform.sensorAttrMatches( sensorDef, 'observation_type', 'calculated'):

                self.formatCalculatedVar( sensorDef, calculatedVars, dataTimesMs )

            # profile id is an outlier informational var with a value

            elif sensorName == 'profile_id':

                self.formatInformationalVar( sensorDef, profileId )

            # some entries in sensor defs require no data (informational only)

            else:
                self.formatInformationalVar( sensorDef, None )

        return 0


    def computeCalculatedVars( self, profileId, profileStartTime, profileEndTime,
                               allProfileBounds, data, gpsData):
        """
        Custom calculations and calibrations for all supported sensors and GliderDac
        output parameters. The current list includes profile avg time, latitude and longitude;
        depth averaged current time, latitude longitude, current_north and current east,
        as required by the GliderDac. The following are calculated if the sensor
        is present in the configured sensor defs: density, irradiance, dissolved oxygen,
        current north and current east.

        :param profileStartTime:
        :param profileEndTime:
        :param data:
        :param gpsData:
        :return: dictionary of calculated vars: { 'varname': {} }
        """

        calculatedVars = {}

        #
        # current components (north, east)
        #

        sensorDef = remus600Platform.getSensorDefFromCfg( self.sensorsCfg, 'current_eastward' )
        if sensorDef:
            instrData = self.getProfileData( sensorDef, data, gpsData,
                                             profileStartTime, profileEndTime)
            currentEast, currentNorth = self.dataProcessor.calculateCurrentComponents(
                instrData['averageCurrent'], instrData['averageDirection'] )
            dataTimesMs = data.timesInMillisecs( instrData.get('timestamp'),
                                                 instrData.get('missionTime'))
            calculatedVars['current_eastward'] = {'values': currentEast, 'times': dataTimesMs}
            calculatedVars['current_northward'] = {'values': currentNorth, 'times': dataTimesMs}
        else:
            logging.warning('Missing sensor current_eastward in sensor_defs config,' +
                            ' required for current calculations')
            raise Exception("Invalid profile, see log file for details")

        #
        # Profile avg time, latitude, longitude
        #

        sensorDef = remus600Platform.getSensorDefFromCfg( self.sensorsCfg, 'profile_time' )
        if sensorDef:
            profileData = remus600Platform.getDataSlice( gpsData, profileStartTime, profileEndTime )
            if not profileData.empty:

                profileTime, profileLat, profileLon = self.dataProcessor.findMidpointTimeLatLon(
                    profileData.get('timestamp'), profileData.get('latitude'), profileData.get('longitude'))

                calculatedVars['profile_time'] = {'values': profileTime, 'times': None}
                calculatedVars['profile_lat'] = {'values': profileLat, 'times': None}
                calculatedVars['profile_lon'] = {'values': profileLon, 'times': None}
            else:
                logging.warning('No gps data bounding profile, assuming abbreviated profile')
                raise Exception("Invalid profile, see log file for details")

        else:
            logging.warning('Missing sensor profile_time in sensor_defs config,' +
                            ' required for profile calculations')
            raise Exception("Invalid profile, see log file for details")

        #
        # Depth avg current time, lat, lon, east current, north current
        #

        sensorDef = remus600Platform.getSensorDefFromCfg( self.sensorsCfg, 'time_uv' )
        if sensorDef:
            adcpData = self.getProfileData( sensorDef, data, gpsData,
                profileStartTime, profileEndTime)

            # compute over full dive (may precede or follow profile),
            # handle endpoint cases

            profileIndex = profileId - 1
            previousProfileIndex = profileIndex - 1
            if previousProfileIndex < 0:
                previousProfileIndex = 0
            nextProfileIndex = profileIndex + 1
            if nextProfileIndex == len( allProfileBounds ):
                nextProfileIndex = profileIndex

            adcpCfg = remus600Platform.getInstrumentFromCfg(
                    self.instrumentsCfg, sensorDef['attrs']['instrument'])

            uvTime, uvLat, uvLon, u, v = self.dataProcessor.calculateUV(
                data, int( adcpCfg['attrs']['subset_msg_id'] ), 0.0,
                allProfileBounds[ profileIndex ],
                allProfileBounds[ previousProfileIndex ],
                allProfileBounds[ nextProfileIndex])

            calculatedVars['time_uv'] = {'values': uvTime, 'times': None }
            calculatedVars['lat_uv'] = {'values': uvLat, 'times': None }
            calculatedVars['lon_uv'] = {'values': uvLon, 'times': None }
            calculatedVars['u'] = {'values': u, 'times': None }
            calculatedVars['v'] = {'values': v, 'times': None }
        else:
            logging.warning('Missing sensor time_uv in sensor_defs config,' +
                            ' required for depth-avg current calculations')
            raise Exception("Invalid profile, see log file for details")

        #
        # some CTDs used on Remus AUV's have depth, not pressure
        # so always compute pressure from depth
        #

        sensorDef = remus600Platform.getSensorDefFromCfg(self.sensorsCfg, 'pressure' )
        if sensorDef:
            instrData = self.getProfileData( sensorDef, data, gpsData,
                                             profileStartTime, profileEndTime )
            pressure = self.dataProcessor.depthToPressure( instrData['depth'], instrData['latitude'] )
            dataTimesMs = data.timesInMillisecs( instrData.get('timestamp'),
                                                 instrData.get('missionTime') )
            calculatedVars['pressure'] = { 'values': pressure, 'times': dataTimesMs }
        else:
            logging.warning('Missing sensor pressure in sensor_defs config,' +
                            ' required for pressure calculation')
            raise Exception("Invalid profile, see log file for details")
        

        #
        # density
        #

        sensorDef = remus600Platform.getSensorDefFromCfg(self.sensorsCfg, 'density' )
        if sensorDef:
            instrData = self.getProfileData( sensorDef, data, gpsData,
                                             profileStartTime, profileEndTime )

            # Note: some older Remus platforms use a CTD that has depth instead of pressure
            # If no pressure column exists, compute it from depth. - ppw09212023
            
            if not 'pressure'  in instrData.columns :
                logging.warning('Old-style CTD detected, computing pressure from depth')
                instrPressure = self.dataProcessor.depthToPressure( instrData['depth'],
                                                               instrData['latitude'] )
            else :
                instrPressure = instrData['pressure']
            
            density = self.dataProcessor.calculateDensity(
                instrData['salinity'], instrData['temperature'], instrPressure,
                instrData['latitude'], instrData['longitude'] )
            dataTimesMs = data.timesInMillisecs( instrData.get('timestamp'),
                                                 instrData.get('missionTime') )
            calculatedVars['density'] = { 'values': density, 'times': dataTimesMs }
        else:
            logging.warning('Missing sensor density in sensor_defs config,' +
                            ' required for density calculations')
            raise Exception("Invalid profile, see log file for details")

        #
        # irradiance
        #

        sensorDef = remus600Platform.getSensorDefFromCfg( self.sensorsCfg, 'PAR' )
        if sensorDef:
            instrData = self.getProfileData( sensorDef, data, gpsData,
                                             profileStartTime, profileEndTime )
            par = self.dataProcessor.processPARData(
                instrData['sensorVoltage'],
                sensorDef['attrs']['calibration_dark_offset'],
                sensorDef['attrs']['calibration_scale_factor'] )
            dataTimesMs = data.timesInMillisecs( instrData.get('timestamp'),
                                                 instrData.get('missionTime') )
            calculatedVars['PAR'] = { 'values': par, 'times': dataTimesMs }
        else:
            logging.warning('Missing sensor PAR in sensor_defs config,' +
                            ' required for irradiance calculations')
            raise Exception("Invalid profile, see log file for details")

        #
        # CDOM
        #

        sensorDef = remus600Platform.getSensorDefFromCfg( self.sensorsCfg, 'cdom' )
        if sensorDef:
            instrData = self.getProfileData( sensorDef, data, gpsData,
                                             profileStartTime, profileEndTime )
            corrCDOM = self.dataProcessor.processCDOMData(
                instrData[ sensorDef['attrs']['subset_field'] ],
                sensorDef['attrs']['calibration_dark_offset'],
                sensorDef['attrs']['calibration_scale_factor'],
                sensorDef['attrs'].get('correction_scale_factor', 1.0))

            dataTimesMs = data.timesInMillisecs( instrData.get('timestamp'),
                                                 instrData.get('missionTime') )
            calculatedVars[ 'cdom' ] = { 'values': corrCDOM, 'times': dataTimesMs }
        else:
            logging.warning('Missing sensor CDOM in sensor_defs config,' +
                            ' required for CDOM calculations')
            raise Exception("Invalid profile, see log file for details")

        
        #
        # Clorophyll
        #

        sensorDef = remus600Platform.getSensorDefFromCfg( self.sensorsCfg, 'chlorophyll_a' )
        if sensorDef:
            instrData = self.getProfileData( sensorDef, data, gpsData,
                                             profileStartTime, profileEndTime )
            corrChl = self.dataProcessor.processChlorophyllData(
                instrData[ sensorDef['attrs']['subset_field'] ],
                sensorDef['attrs']['calibration_dark_offset'],
                sensorDef['attrs']['calibration_scale_factor'] )
            dataTimesMs = data.timesInMillisecs( instrData.get('timestamp'),
                                                 instrData.get('missionTime') )
            calculatedVars[ 'chlorophyll_a' ] = { 'values': corrChl, 'times': dataTimesMs }
        else:
            logging.warning('Missing sensor chlorophyll_a in sensor_defs config,' +
                            ' required for chlorophyll_a calculations')
            raise Exception("Invalid profile, see log file for details")

        
        #
        # dissolved oxygen
        #

        sensorDef = remus600Platform.getSensorDefFromCfg( self.sensorsCfg, 'dissolved_oxygen' )
        if sensorDef:
            instrData = self.getProfileData( sensorDef, data, gpsData,
                                             profileStartTime, profileEndTime )
            o2 = self.dataProcessor.processOxygenData(
                instrData['concentration'], instrData['salinity'],
                instrData['depth'], instrData['temperature'],
                instrData['latitude'], instrData['longitude'] )
            dataTimesMs = data.timesInMillisecs( instrData.get('timestamp'),
                                                 instrData.get('missionTime') )
            calculatedVars['dissolved_oxygen'] = { 'values': o2, 'times': dataTimesMs }
        else:
            logging.warning('Missing sensor dissolved_oxygen in sensor_defs config,' +
                            ' required for oxygen calculations')
            raise Exception("Invalid profile, see log file for details")

        return calculatedVars

    def getProfileData( self, sensorDef, data, gpsData, profileStartTime, profileEndTime):
        """
        If a sensor has an associated instrument, returns the data for that instrument
        that falls within the profile time bounds. If the instrument is GPS, the data
        is retrieved from the previously recomputed 1 second resolution gps cache instead.
        :param sensorDef:
        :param data:
        :param gpsData:
        :param profileStartTime:
        :param profileEndTime:
        :return: data (DataFrame) within the profile time bounds
        """

        profileData = None

        if remus600Platform.sensorHasAttr(sensorDef, 'instrument'):

            # Get instrument data within the profile time bounds
            # Read input data unless gps, which was recomputed at 1 sec cadence

            instrCfg = remus600Platform.getInstrumentFromCfg(
                self.instrumentsCfg, sensorDef['attrs']['instrument'] )

            if instrCfg['nc_var_name'] != 'instrument_gps':

                profileData = data.getDataSliceForMessageId(
                    int( instrCfg['attrs']['subset_msg_id'] ),
                    profileStartTime, profileEndTime )

            else:
                profileData = remus600Platform.getDataSlice(
                    gpsData, profileStartTime, profileEndTime )

        return profileData


    def formatGlobalAttributes( self, profileId ):
        """
        Create output attributes for all configured global attributes
        :return: None
        """

        # Create output attributes for all entries in global_attributes configuration

        for name, value in self.globalsCfg.items():
            self.outputFileWriter.addGlobalAttr( name, value )

        # Some data and attributes found in deployment cfg

        trajectoryName = self.deploymentCfg['trajectory_name']
        trajectoryDateTime = self.deploymentCfg['trajectory_datetime']
        vehicle = self.deploymentCfg['glider']
        if 'global_attributes' in self.deploymentCfg:
            comment = self.deploymentCfg['global_attributes']['comment']
            deploymentNumber = self.deploymentCfg['global_attributes']['deployment_number']
            wmoId = self.deploymentCfg['global_attributes']['wmo_id']
        else:
            comment = " "
            deploymentNumber = " "
            wmoId = " "

        self.outputFileWriter.addGlobalAttr( 'comment', comment )
        self.outputFileWriter.addGlobalAttr( 'id', trajectoryName + "_" + str(profileId) )
        self.outputFileWriter.addGlobalAttr( 'title', trajectoryName )
        self.outputFileWriter.addGlobalAttr('wmo_id', wmoId )

        nowUtc = datetime.datetime.utcnow()
        nowUtcString = nowUtc.strftime('%Y-%m-%dT%H:%M:%SZ')
        self.outputFileWriter.addGlobalAttr('history', nowUtcString + ': ./profileDataFormatter.py')

        self.outputFileWriter.addGlobalAttr('date_created', nowUtcString )
        self.outputFileWriter.addGlobalAttr('date_modified', nowUtcString )
        self.outputFileWriter.addGlobalAttr('date_issued', nowUtcString )


    def formatInstrumentVars( self ):
        """
        Create an output variable for each instrument in instrument configuration
        :return: None
        """

        # Create an output variable for each instrument

        for instrCfg in self.instrumentsCfg:
            self.outputFileWriter.addVariable(
                instrCfg["nc_var_name"], instrCfg["type"],
                None, instrCfg['attrs'], 0, None)


    def formatMeasuredVar(self, sensorDef, profileData, profileTimes):
        """
        Create an output variable for a measured sensor definition
        :param sensorDef:
        :param profileData:
        :param profileTimes:
        :return: None
        """

        self.outputFileWriter.addVariable(
            sensorDef['nc_var_name'],
            sensorDef['type'],
            sensorDef['dimension'],
            sensorDef['attrs'],
            profileData[sensorDef['attrs']['subset_field']],
            profileTimes )

        # Variables need a corresponding quality control indicator variable
        self.formatQCVar( sensorDef, profileTimes )

    def formatCalculatedVar( self, sensorDef, calculatedData, profileTimes ):
        """
        Create an output variable for a calculated sensor def and calculation results
        :param sensorDef:
        :param calculatedData:
        :param profileTimes:
        :return: None
        """

        # Find calculated value(s) for passed sensor

        if sensorDef['nc_var_name'] in calculatedData:
            self.outputFileWriter.addVariable(
                sensorDef['nc_var_name'],
                sensorDef['type'],
                sensorDef['dimension'],
                sensorDef['attrs'],
                calculatedData[ sensorDef['nc_var_name'] ]['values'],
                calculatedData[ sensorDef['nc_var_name'] ]['times'] )

            # Variables need a corresponding quality control indicator variable
            self.formatQCVar( sensorDef, calculatedData[ sensorDef['nc_var_name'] ]['times'] )

        else:
            logging.warning('Encountered unexpected calculated sensor: ' +
                            sensorDef['nc_var_name'] + ' ignored.')


    def formatInformationalVar(self, sensorDef, sensorValue ):
        """
        Create an output variable for a sensor definition having no data or time dimension
        :param sensorDef:
        :return: None
        """

        # The "platform" sensor is an outlier, in that there exists
        # some "platform" attributes within the deployment config
        # Add those non-duplicates attributes to the attributes
        # in the "platform" sensor def four output

        sensorAttrs = sensorDef['attrs'].copy()
        if sensorDef['nc_var_name'] == 'platform':
            for key, value in self.deploymentCfg['platform'].items():
                if not key in sensorAttrs:
                    sensorAttrs[key] = value

        self.outputFileWriter.addVariable(
            sensorDef['nc_var_name'], sensorDef['type'], None,
            sensorAttrs, sensorValue, None )

        # Scalar variables need a corresponding quality control indicator variable
        if remus600Platform.sensorAttrMatches( sensorDef, 'type', 'platform') == False and \
            remus600Platform.sensorAttrMatches( sensorDef, 'type', 'instrument') == False and \
            sensorDef['nc_var_name'] not in ['crs', 'profileId' ]:
               self.formatQCVar( sensorDef, None )

    def formatQCVar(self, sensorDef, qcTimes ):
        '''
        Generate a quality control variable for passed sensor.
        Value zero indicates no quality control performed
        :param sensorDef:
        :param qcTimes:
        :return: None
        '''

        attrs = { '_FillValue': -127,
                  'flag_meanings': 'no_qc_performed good_data probably_good_data' +
                                   ' bad_data_that_are_potentially_correctable' +
                                   ' bad_data value_changed not_used not_used' +
                                   ' interpolated_value missing_value',
                  'flag_values': np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9],dtype=np.int8),
                  'long_name': sensorDef['nc_var_name'] + ' quality flag',
                  'valid_max': np.int8(9),
                  'valid_min': np.int8(0) }

        if qcTimes is None:
            values = 0
        else:
            values = np.zeros( len( qcTimes ))

        self.outputFileWriter.addVariable(
            sensorDef['nc_var_name'] + '_qc',
            'byte',
            sensorDef['dimension'],
            attrs,
            values,
            qcTimes )

    def purgeOutOfRangeData( self, profileId, profileData, sensorDef, dateTimeMs ):
        '''
        If profile data contains out of range values, log them and
        replace them with the sensor's fill value
        :param profileId (for logging)
        :param profileData
        :param sensorDef
        :param dateTimeMs (for logging)
        :return: None
        '''

        if not remus600Platform.sensorHasAttr(sensorDef, '_FillValue'):
            return

        profileDataColumn = sensorDef['attrs']['subset_field']
        if not profileDataColumn in profileData.columns:
            return

        instrument = "instrument unspecified"
        if remus600Platform.sensorHasAttr(sensorDef, 'instrument'):
            instrument = sensorDef['attrs']['instrument']

        if remus600Platform.sensorHasAttr(sensorDef, 'valid_min'):

            outOfRangeLow = profileData[profileDataColumn] < sensorDef['attrs']['valid_min']

            if profileData[profileDataColumn][outOfRangeLow].size > 0:

                # TBD - eventually stub out logging, leave in to get a handle on bad data
                replacedStr = ''
                for index, value in profileData[profileDataColumn][outOfRangeLow].items():
                    replacedStr = replacedStr + str(dateTimeMs[index]) + ' ' + \
                    str(profileData[profileDataColumn][index]) + '\n'
                    
                logging.warning('Profile ' + str(profileId) + ', ' + instrument +
                                ', sensor ' + profileDataColumn +
                                ' contains low out of range data. ' +
                                'Replacing with _FillValue ' +
                                str(sensorDef['attrs']['_FillValue']) +
                                ' at: \n' + replacedStr)

                # Conditionally replacing dataframe values requires odd syntax below
                profileData.loc[ outOfRangeLow, profileDataColumn] = sensorDef['attrs']['_FillValue']

        if remus600Platform.sensorHasAttr(sensorDef, 'valid_max'):

            outOfRangeHigh = profileData[ profileDataColumn ] > sensorDef['attrs']['valid_max']

            if profileData[profileDataColumn][outOfRangeHigh].size > 0:

                # TBD - eventually stub out logging, leave in to get a handle on bad data
                replacedStr = ''
                for index, value in profileData[profileDataColumn][outOfRangeHigh].items():
                    replacedStr = replacedStr + str(dateTimeMs[index]) + ' ' + \
                    str(profileData[profileDataColumn][index]) + '\n'
                    
                logging.warning('Profile ' + str(profileId) + ', ' + instrument +
                                ', sensor ' + profileDataColumn +
                                ' contains high out of range data. ' +
                                'Replacing with _FillValue ' +
                                str(sensorDef['attrs']['_FillValue']) + ' at: \n' +
                                replacedStr)

                # Conditionally replacing dataframe values requires odd syntax below
                profileData.loc[ outOfRangeHigh, profileDataColumn] = sensorDef['attrs']['_FillValue']
