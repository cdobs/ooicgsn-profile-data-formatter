# README #

New notes:
01/16/2026 - Added an option to adjust CDOM values affected by a vendor-identified issue that caused bias due to improperly prepared calibration standards. The issue is described in further detail at https://oceanobservatories.org/2024/12/sbs-issues-notice-for-certain-cdom-fluorometers/. To adjust CDOM values, update the correction_scale_factor value in the flort_cdom section of sensor_defs.cfg file. The default is 1.0.  

09/25/2023 - Several early AUV deployments use the old Seabird CTD only. No Neil Brown CTD was installed. For those deployments, the subset_message_id field for the CTD instrument in the instruments.cfg file must be changed from 1107 to 1181. Additionally, in the sensor_defs.cfg file, the observation_type for the pressure sensor def must be changed from "measured" to "calculated", as this CTD's message data does not contain a pressure column.

05/21/2023 - To calculate CDOM using dark offset and scale factors from the sensor_defs.cfg, you must also change the "observation_type" attribute to "calculated", as shown below. This note also holds true for chlorophyll_a, which is now also calculated in the same way.

![sensor_defs.json](sensor_defs_cdom.png "sensor_defs") 


### What is this repository for? ###

* This is the ProfileDataFormatter software development repo, containing the ProfileDataFormatter python application. The app is used to to format data from profiling mobile platforms (such as gliders and AUVs) into netCDF files for import into oceanographic instrument data repositories, such as the IOOS-DAC (https://ioos.github.io/ioosngdac/) and OOI Data Explorer (https://dataexplorer.oceanobservatories.org/). The application is derived from work done by Stuart Pearce at OSU <https://github.com/s-pearce/gliderdac/wiki> processing Slocum Glider data for import into IOOS-DAC. That gliderdac application is imported as a git subtree under legacy in this repository. The gliderdac code is used "as is" to implement formatting for the Slocum glider.

* The ProfileDataFormatter is an enhancement to gliderdac in that supports both multiple input and output data formats. It currently accepts either Slocum 2.0 Glider or Remus 600 subset data files as input and produces either or both GliderDAC compatible NetCDF files and Data Explorer compatible NetCDF files as output. Using the older gliderdac software, only GliderDAC compatible NetCDF files can be produced directly. To generate Data Explorer compatible files, output needed to be ingested into the GliderDAC then subsequently exported from GliderDAC into NetCDF files compatible with Data Explorer.

### Installation ###

    > Assumption: miniconda python 3.8, 64 bit environment
    > git clone https://github.com/WHOIGit/ooicgsn-profile-data-formatter.git
    > cd ooicgsn-profile-data-formatter
    > conda env create --prefix ./pdfenv -f environment.yml
    > conda activate ./pdfenv

### Usage ###

* The ProfileDataFormatter is a command prompt driven application supporting a number of passed parameters, as shown below:

-h  
   help (optional),
   displays a list of command-line options
   
-c {path}  
   configuration path (required)  
   directory in which configuration files are to be found

-d {data file(s)}  
   One or more data files to be processed (required)  
   Wildcards are supported using glob syntax

-m "Slocum Glider 2.0" or "Remus 600 AUV"  
   Mobile platform that is the source of input data  
   (optional, default is "Slocum Glider 2.0")

-t "IOOS-DAC" or "OOI-EXPLORER"  
   Target repository for which to format data file(s)  
   (optional, default is "IOOS-DAC")
   Note that selection of OOI-EXPLORER will produce output files for both OOI and DAC. Smaller, profile specific output files of the format "{trajectoryName}{profileTime}_delayed.nc" or "{trajectoryName}{profileTime}_rt.nc" are intended for DAC. The larger, full trajectory NetCDF file "{trajectoryName}{trajectoryTime}.nc" is intended for OOI Explorer.

-p { mobile platform specific parameters }  
   Platform specific arguments (optional)
   
   Syntax: JSON dictionary string '{"paramName": "paramValue", ...}'  
   ie:  "{'ctd_sensor_prefix' : 'sci', 'start_profile_id' : 15}"
   
   For Slocum Glider 2.0, the following are supported:

   - 'ctd_sensor_prefix' : 'sci' or 'm'  [default 'sci']

   - 'start_profile_id' : n  [default: 0, implies use unix timestamp

-o {path}  
   output path (optional, default is '.')  
   Path into which output files are written
   
-k
   Clobber flag  
   Indicates that existing output files for the same trajectory  
   are to be overwritten

-f "NETCDF3_CLASSIC" or "NETCDF4_CLASSIC" or "NETCDF4"  
   NetCDF file format to be written (optional, default is NETCDF4_CLASSIC)

-cl {0,1,2,3,4,5,6,7,8,9}  
   Compression level for output NetCDF file (optional, default is 1)

-s  
   Suppress output.  
   Verifies configuration files without writing output.

-l {debug,info,warning,error,critical}  
   Log level (optional, default is info)  
   Log file is ProfileDataFormatter.log, written to the current working directory. The file is appended for each new run, with newest log entries at the end of the file.

In addition to the above command-line options, the ProfileDataFormatter utilizes the same configuration files and settings as its gliderdac predecessor: sensor_defs.json, global_attributes.json, instruments.json and deployment.json. Contents of these configuration files are described as follows:

*sensor_defs.json* - contains items referred to as sensor_defs that correspond to NetCDF output variables and their attributes. The list of NetCDF variables, their definitions and required attributes and values are defined by the IOOS GliderDAC repository at the previously defined link.

   * Note that for support of Remus AUV Subset data files as input, a sensor_def attribute, "subset_field" has been added to the list of required attributes for sensor_defs having an observation_type attribute with values of either "measured" or "calculated". This defines the column name in the subset data file identifying the value to be used in the message type corresponding to the sensor of interest.

*global_attributes.json* - contains the static global attributes that are written to the global attributes in the NetCDF file. The list of required global attributes are defined by the IOOS GliderDAC repository at the link above.

*instruments.json* - contains the set of scientific instruments aboard the mobile platform that is the source of the data files to be processed. All metadata, such as instrument calibration, serial numbers and attributes are included in this file.

   * Note that for support of Remus AUV Subset data files as input, an instrument attribute named "subset_msg_id" has been added to the attribute list for each instrument reporting data in the file. The value of the subset_msg_id attribute is to be taken from the most recently available version of the "Data Subsetting Manual", published and distributed by Hydroid (Kongsberg Maritime).

*deployment.json* - contains a dictionary of attributes describing the specific deployment of the mobile platform being processed


Note that Slocum 2.0 Glider files are expected to be in the form of merged, ascii file format as described in the gliderdac usage instructions at the link above.

### Examples ###

Example files for validating the installation of the ProfileDataFormatter are provided in the repository. Calling parameters, configuration settings and input data files from both Slocum Glider and Remus AUV mobile platforms are supplied. All test files can be found under {installation directory}/tests, as follows:

* tests/test_parameters.txt - calling sequences for each of the test cases. Paste and run each from the ProfileDataFormatter's root installation directory

* tests/config/slocum200, tests/config/remus600 - each directory contains the four required json configuration files for a particular mobile platform deployment. These are deployment.json, global_attributes.json, instruments.json and sensor+defs.json

* tests/data - contains Slocum200 ascii merge data file used as input to the Slocum200 test cases

* tests/auvdata - contains Remus600 subset data files used to drive the Remus600 test cases

* tests/output, tests/output/auv - empty paths to which test case output gets written

When each test is run, in addition to NetCDF output files being generated, a log file, ProfileDataFormatter.log will be written to the directory from which the application is run. The file should be examined for errors and warnings.


### Support ###

* For questions, please email paul.whelan@whoi.edu

### Notes ###

 The soft links "configuration.py", "profile_filters.py" and "ooidac" in the root of this project
 point at ./legacy/gliderdac code and directories due to the architecture of the legacy gliderdac application. 
 They are referenced from nested packages within gliderdac, as though it were at the 
 root of the PYTHONPATH. These were required to utilize the gliderdac code without modification.
