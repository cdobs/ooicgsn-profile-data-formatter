"""
class: remusXYProcessor

description:

history:
09/21/2021 ppw created
01/16/2026 cdobs - adds cdom correction factor
"""
import logging
import numpy as np
import numexpr as ne
from pandas import DataFrame, Series
from scipy import signal
from DataProcessor.AuvProcessor.auvProcessor import auvProcessor
from common.constants import OCEAN_DEPTH_M
from gsw import SP_from_C, SA_from_SP, CT_from_t, rho, z_from_p, p_from_z


class remus600Processor( auvProcessor ) :

    def __init__( self ) :
        super().__init__()

        # constants (move to config?)
        self.MIN_PROFILE_DEPTH_METERS = 10
        self.MIN_PROFILE_TIME_SECONDS = 60
        self.MAX_TIME_GAP_SECONDS = 30
        self.SMOOTHING_WINDOW_SIZE = 150

    def computeProfiles(self, trajectoryData, timeIndex, missionTimeIndex, depthIndex):
        """
        Given trajectory data for any instrument containing time, mission time
        and depth columns, computes a series of profile start and end times.

        :param trajectoryData: pandas dataframe containing time, mission time and depth
        :param timeIndex: index identifier of the dataframe time column
        :param missionTimeIndex: index identifier of the dataframe mission time column
        :param depthIndex: index identifier of the dataframe depth column
        :return: pandas dataframe containing 'startTime' and 'endTime' columns
        """

        # convert pandas dataframe series into ndarrays for processing

        times = np.asarray( trajectoryData[timeIndex] )
        missionTimes = np.asarray( trajectoryData[missionTimeIndex] )
        depths = np.asarray( trajectoryData[depthIndex] )

        # Combine time fields into millisecs since trajectory start
        # (Time in seconds, missionTime in millisecs into day)

        startTime = times[0]
        timeOffsets = self.computeTimeOffsets(
            startTime, times, missionTimes )

        # Create clean depths ( 0 <= depth <= ocean depth )

        times, depths = self.cleanDepthData(
            timeOffsets, depths )

        # Smooth depths

        smoothingSize = 10
        times, depths = self.smoothDepths(
            times, depths, smoothingSize)

        # compute inflection points (profile bounds)

        profileBounds = self.computeInflectionTimes( times,
                                                     depths,
                                                     startTime,
                                                     self.MIN_PROFILE_DEPTH_METERS,
                                                     self.MIN_PROFILE_TIME_SECONDS,
                                                     self.MAX_TIME_GAP_SECONDS)

        return profileBounds

    def computeTimeOffsets(self, startTime, timestamps, dayOffsets):
        """
        Combine timestamps, dayOffsets into a single offset from start of mission
        :param startTime: unixtime defining start of mission
        :param timestamps: unixtime (secs) for each sample
        :param dayOffsets: offset (msecs) into day
        :return: offsets (ms) from start of mission
        """

        timeOffsets = (1000 * (timestamps - startTime)) + \
            np.fmod( dayOffsets, 1000 )

        return timeOffsets


    def cleanDepthData(self, timeOffsets, rawDepths ):
        """
        Excise samples where depth < 0 or unreasonably large
        :param timeOffsets:
        :param rawDepths:
        :return: updated times, depths
        """

        # Create clean depths ( 0 <= depth <= ocean depth )
        cleanDepths = rawDepths
        cleanDepths = np.where( cleanDepths > OCEAN_DEPTH_M, np.nan, cleanDepths)
        cleanDepths = np.where( cleanDepths < 0.0, np.nan, cleanDepths)

        # Filter bad depths out of time and depth
        depthIndices = np.isfinite( cleanDepths )
        timeOffsets = timeOffsets[ depthIndices ]
        cleanDepths = cleanDepths[ depthIndices ]

        return timeOffsets, cleanDepths

    def smoothDepths(self, times, depths, windowSize):
        """
        Perform boxcar smoothing on depth data
        :param times:
        :param depths:
        :param windowSize:
        :return: times and depths with smoothing applied
        """

        window = signal.windows.boxcar(windowSize)
        smoothedDepths = signal.convolve(depths, window, 'same') / windowSize

        # remove the extra points with filter edge effects
        smoothedDepths = smoothedDepths[windowSize:-(windowSize+1)]
        timesOut = times[windowSize:-(windowSize+1)]

        return timesOut, smoothedDepths

    def computeInflectionTimes(self, times, depths, missionStartTime, minDeltaD,
                               minDeltaTSecs, maxTimeGapSecs):
        """
        Find profile bounds, eliminating periods of level flight, short profiles and time gaps
        :param times: offset in milliseconds from trajectory start
        :param depths: meters
        :param missionStartTime: unixtime UTC
        :param minDeltaD: min depth change for a valid profile (meters)
        :param minDeltaTSecs: min time duration for valid profile (seconds)
        :param maxTimeGapSecs: max time gap valid within a profile
        :return: array of valid profile [start, end] times in UnixTime UTC seconds
        """

        logging.info("computeInflectionTimes")
        logging.info( "Mission start time: " + str(missionStartTime))

        minDeltaTMillisecs = 1000 * minDeltaTSecs

        # compute rate of depth change and direction
        # note: have seen consecutive, duplicate times; handle it
        # to avoid divide by zero
        dt = np.diff( times )
        dt = np.where( dt == 0, 1, dt )
        dZdT = np.diff( depths ) / dt
        updownlevel = np.sign( dZdT )

        # init return array of [start, end] items
        profileBounds = []

        # traverse time/depth looking for direction changes and time gaps
        start = 0
        end = 1
        while end < len(updownlevel):

            if (times[end] - times[end - 1] > 1000 * self.MAX_TIME_GAP_SECONDS):
                logging.info(" gap between " + str(end-1) + ' and ' + str(end))

            # if time gap too large or direction changes, end current profile here

            if (times[end] - times[end-1] > 1000 * self.MAX_TIME_GAP_SECONDS) or \
                    (updownlevel[end] != 0 and updownlevel[end] != updownlevel[end - 1]):
                if np.fabs(depths[end-1] - depths[start]) >= minDeltaD and \
                    times[end-1] - times[start] >= minDeltaTMillisecs:
                    profileBounds.append( [ times[start], times[end-1] ] )
                    logging.info("profile between " + str(start) + ' and ' + str(end-1))

                start = end
            end = end + 1

        # if incomplete profile at end of data, complete it

        if start < len(depths)-1:
            if np.fabs(depths[end - 1] - depths[start]) >= minDeltaD and \
                    times[end - 1] - times[start] >= minDeltaTMillisecs:
                profileBounds.append( [ times[start], times[end-1] ] )
                logging.info("profile between " + str(start) + ' and ' + str(end - 1))

        # Convert profile bounds back to unixtime UTC

        profileBoundsRA = np.asarray(profileBounds,dtype=int)
        return self.offsetTimesToEpochSecs( missionStartTime, profileBoundsRA )


    def offsetTimesToEpochSecs(self, missionStartSecs, profileBoundsRA):
        """
        convert millisec offset from mission start to unix time utc
        :param missionStartSecs:
        :param profileBoundsRA:
        :return:
        """
        return profileBoundsRA / 1000 + missionStartSecs

    def interpolateGpsData(self, gpsData ):
        """
        gps data at 1 sec resolution while surfaced, blank when not;

        It turns out that this gps instrument continues to report the lat/lon values
        for the most recent surface fix when submerged, simply updating the FixAge variable.
        Soooo, we now have to filter the values we interpolate from with (FixAge == 0)

        fill to 1 sec resolution through interpolation
        :param gpsData:
        :return: filledGpsData
        """

        # Filter incoming gps data to those values where FixAge == 0

        rows = gpsData.loc[ gpsData['fixAge'] == 0]
        filteredTime = rows['timestamp']
        filteredLat = rows['latitude']
        filteredLon = rows['longitude']

        #  create new 1 second res timestamp array (no gaps)

        allTimestamps = np.arange( int(gpsData['timestamp'][0]),
                                   int(gpsData['timestamp'].iloc[-1] + 1))

        # Interpolate lat, lon to 1 second resolution in gaps

        #allLatitudes = np.interp( allTimestamps, gpsData['timestamp'], gpsData['latitude'] )
        #allLongitudes = np.interp( allTimestamps, gpsData['timestamp'], gpsData['longitude'] )
        allLatitudes = np.interp( allTimestamps, filteredTime, filteredLat )
        allLongitudes = np.interp( allTimestamps, filteredTime, filteredLon )

        filledGpsData = { 'timestamp': allTimestamps,
                          'latitude': allLatitudes,
                          'longitude': allLongitudes }

        return DataFrame( filledGpsData )


    def depthToPressure(self, depth, lat) :
        """
        Convert depth (in meters) to pressure (in dB)
        """

        raDepth = depth.to_numpy()
        raLat = lat.to_numpy()

        # depth to pressure
        pressure = p_from_z( -raDepth , raLat)

        return Series( pressure )


    def processOxygenData(self, rawO2Concentration, salinity, depth, temp, lat, lon, pref=0):
        """
        Note:
        This salinity correction was swiped directly from OOI CI processing functions in:
        https://github.com/oceanobservatories/ion-functions/blob/master/ion_functions/data/do2_functions.py.
        If anything else there is used, that package should be made a dependency and this function
        replaced with a call to it. (Too many nested dependencies of the ion package that are not
        required by this function alone to justify at this time.)

        :param rawO2Concentration - uncorrected O2 (uM)
        :param salinity - practical salinity
        :param depth - depth (m)
        :param temp - temperature (deg C)
        :param lat - latitude (deg)
        :param lon - longitude (deg)
        :param pref=0     - pressure reference, default to 0
        :return: correctedO2 - (uM)

        Implemented by:
            2013-04-26: Stuart Pearce. Initial Code.
            2015-08-04: Russell Desiderio. Added Garcia-Gordon reference.
        References:
            OOI (2012). Data Product Specification for Oxygen Concentration
                from "Stable" Instruments. Document Control Number
                1341-00520. https://alfresco.oceanobservatories.org/ (See:
                Company Home >> OOI >> Controlled >> 1000 System Level
                >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)
            "Oxygen solubility in seawater: Better fitting equations", 1992,
            Garcia, H.E. and Gordon, L.I. Limnol. Oceanogr. 37(6) 1307-1312.
            Table 1, 5th column.
        """

        # Issues with running under python 3.8, & compatible libraries, gsw hates Dataframe/Series. Use np.array

        raRawO2Concentration = rawO2Concentration.to_numpy()
        raSalinity = salinity.to_numpy()
        raDepth = depth.to_numpy()
        raTemp = temp.to_numpy()
        raLat = lat.to_numpy()
        raLon = lon.to_numpy()

        # depth to pressure
        pressure = p_from_z( -raDepth , raLat)

        # density calculation from GSW toolbox
        SA = SA_from_SP( raSalinity, pressure, raLon, raLat)
        CT = CT_from_t(SA, raTemp, pressure)
        pdens = rho(SA, CT, pref)  # potential referenced to p=0

        # Convert from volume to mass units:
        DO = ne.evaluate('1000*raRawO2Concentration/pdens')

        # Pressure correction:
        DO = ne.evaluate('(1 + (0.032*pressure)/1000) * DO')

        # Salinity correction (Garcia and Gordon, 1992, combined fit):
        S0 = 0
        ts = ne.evaluate('log((298.15-raTemp)/(273.15+raTemp))')
        B0 = -6.24097e-3
        B1 = -6.93498e-3
        B2 = -6.90358e-3
        B3 = -4.29155e-3
        C0 = -3.11680e-7
        Bts = ne.evaluate('B0 + B1*ts + B2*ts**2 + B3*ts**3')
        DO = ne.evaluate('exp((raSalinity-S0)*Bts + C0*(raSalinity**2-S0**2)) * DO')

        # convert back to volume
        # DO = ne.evaluate('DO*pdens/1000')

        return DO

    def processPARData(self, sensorVoltage, calibratedDarkOffset, calibratedScaleFactor):
        """
        PAR data needs calculation from voltage using calibration constants
        :param sensorVoltage: output by Biospherical 2150
        :param calibratedDarkOffset - in sensor_defs.json from calibration certificate
        :param calibratedScaleFactor - in sensor_defs.json from calibration certificate
        :return: calculatedPar (umol m-2 s-1)
        """

        # taken from legacy/ooidac/processing.py recalc_par(), same instrument
        calculatedPar = (sensorVoltage - calibratedDarkOffset) / calibratedScaleFactor

        return calculatedPar

    def processCDOMData(self, sensorCounts, calibratedDarkOffset, calibratedScaleFactor, correctionScaleFactor=1.0):
        """
        CDOM data needs calculation from counts using calibration constants
        :param sensorCounts: output by FLBBCD
        :param calibratedDarkOffset - in sensor_defs.json from calibration certificate
        :param calibratedScaleFactor - in sensor_defs.json from calibration certificate
        :param correctionScaleFactor - optional, in sensor_defs.json cdom correction scale factor
        :return: correctedCDOM
        """

        correctedCDOM = (sensorCounts - calibratedDarkOffset) * calibratedScaleFactor * correctionScaleFactor

        return correctedCDOM

    def processChlorophyllData(self, sensorCounts, calibratedDarkOffset, calibratedScaleFactor):
        """
        Clorophyll data needs calculation from counts using calibration constants
        :param sensorCounts: output by FLBBCD
        :param calibratedDarkOffset - in sensor_defs.json from calibration certificate
        :param calibratedScaleFactor - in sensor_defs.json from calibration certificate
        :return: correctedChlorophyll
        """

        correctedChlorophyll = (sensorCounts - calibratedDarkOffset) * calibratedScaleFactor

        return correctedChlorophyll

    def calculateDensity(self, salinity, temperature, pressure, latitude, longitude):
        """Calculates density given practical salinity, temperature, pressure, latitude,
        and longitude using Gibbs gsw SA_from_SP and rho functions.

        Parameters:
            temperature (C), pressure (dbar), salinity (psu PSS-78),
            latitude (decimal degrees), longitude (decimal degrees)

        Returns:
            density (kg/m**3),
        """

        # dBar_pressure = pressure * 10

        absolute_salinity = SA_from_SP(
            salinity,
            pressure,
            longitude,
            latitude
        )

        conservative_temperature = CT_from_t(
            absolute_salinity,
            temperature,
            pressure
        )

        density = rho(
            absolute_salinity,
            conservative_temperature,
            pressure
        )

        return density

    def calculateCurrentComponents(self, currentSpeeds, currentDirections):

        currentEast = currentSpeeds * np.sin( np.deg2rad( currentDirections ) )
        currentNorth = currentSpeeds * np.cos( np.deg2rad( currentDirections ) )

        return currentEast, currentNorth

    def calculateUV(self, data, adcpMsgId, minDepth, currentProfile, previousProfile, nextProfile):
        """
        Computes depth averaged current components w/ mean time, lat and lon
        :param data:        data object containing all trajectory data
        :param adcpMsgId:   msg id associated with adcp data
        :param minDepth:    min depth (m) at which to calculate current
        :param currentProfile: profile start, end times for profile to compute
        :param previousProfile: start and end times for preceding profile
        :param nextProfile: start and end times for succeeding profile
        :return:            time, lat, lon, u, v
        """

        # Get data over the full dive. Either previous profile + current if ascending
        # or current profile + next if descending

        profileData = data.getDataSliceForMessageId( adcpMsgId, currentProfile[0], currentProfile[1])
        if profileData['depth'].iloc[0] < profileData['depth'].iloc[-1]:
            profileData = data.getDataSliceForMessageId(adcpMsgId, currentProfile[0], nextProfile[1])
        else:
            profileData = data.getDataSliceForMessageId(adcpMsgId, previousProfile[0], currentProfile[1])

        atDepthData = np.where( profileData['depth'] >= minDepth)

        currents = profileData['averageCurrent'].iloc[atDepthData]
        directions = profileData['averageDirection'].iloc[atDepthData]
        U, V = self.calculateCurrentComponents( currents, directions )
        u = U.mean()
        v = V.mean()

        times = profileData['timestamp'].iloc[atDepthData]
        latitudes = profileData['latitude'].iloc[atDepthData]
        longitudes = profileData['longitude'].iloc[atDepthData]
        uvTime, uvLat, uvLon = self.findMidpointTimeLatLon( times, latitudes, longitudes )

        return uvTime, uvLat, uvLon, u, v


    def findMidpointTimeLatLon(self, times, latitudes, longitudes):
        """
        Finds mean time, closest latitude and longitude
        :param times:
        :param latitudes:
        :param longitudes:
        :return: time, lat, lon
        """

        midpointTime = times.mean()
        timeIndex = times.sub( midpointTime ).abs().idxmin()
        midpointLat = latitudes[ timeIndex ]
        midpointLon = longitudes[ timeIndex ]

        return midpointTime, midpointLat, midpointLon
