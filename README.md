# MUSES_GPPNPP_estimation_programme

MUSES 500m GPP/NPP Estimation Program Running Instructions

1.	Input Data
   
The input data include remote sensing data, meteorology data, elevation data and land cover data.
Detailed description of all input data:
1). Average air temperature (AAT)： Average air temperature was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: ℃. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of the AAT is “2m_temperature2001001.tif”.
2). Surface solar radiation downwards (SSRD): SSRD was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: MJm-2d-1. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of SSRD is “surface_solar_radiation_downwards2001001.tif”.
3). Dew point temperature: Dew point temperature was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: ℃; Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of the dew point temperature is “2m_dewpoint_temperature2001001.tif”.
4). Surface_net_thermal_radiation (SNTR): SNTR was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: wm-2. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of SNTR is “surface_net_thermal_radiation2001001.tif”.
5). Surface_net_solar_radiation (SNSR): SNSR was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: wm-2. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of SNSR is “surface_net_solar_radiation2001001.tif”.
6). FPAR (fraction of absorbed photosynthetically active radiation): spatial resolution of the FPAR data is 500m, and the temporal resolution is 8 days. Data type: byte. Unit: none. Scale factor = 0.004. Fill value: 255. Projection: Sinusoidal. HDF format. One example file name of FPAR is “MUSES.A2001001.H05V10.MODIS.FAPAR500M.C05.HDF”.
7). LAI (leaf and index): spatial resolution of the LAI data is 500m, and the temporal resolution is 8 days. Data type: integer. Unit: m2/m2. Scale factor = 0.01. Fill value: 2000, 2500. Projection: Sinusoidal. HDF format. One example file name of LAI is “MUSES.A2001001.H00V08.MODIS.LAI500M.C05.HDF”.
8). LAImaxy (Maximum LAI in a year): spatial resolution is 500m. Data type: integer. Unit: m2/m2. Scale factor = 0.01. Fill value: 2000, 2500. Projection: Sinusoidal. Raw format. One example file name of LAImaxy is “LAImaxA2001H00V08.tif”.
9). MODIS landcover: spatial resolution is 500m, Data type: Byte. Unit: none; Projection: Sinusoidal. Raw format. One example file name of landcover is “MCD12Q1.A2001001.h00v08.tif”.
10). DEMH (Digital Elevation Model): spatial resolution is 500m, Data type: Unsigned Integer. Unit: m. Projection: Sinusoidal. TIF format. One example file name of DEMH is “DEM_500M_H00V08. tif”. This is obtained by resampling the original 1km DEM data.
11). DEML (Digital Elevation Model): spatial resolution is 0.25 degrees, Data type: float. Unit: m. Geographic coordinate system: GCS_WGS_1984. TIF format. The example file name of DEML is “GlobalDEM025_1441_721.tif”. This is obtained by aggregating the original 1km DEM data.
12). Longitude and latitude data: spatial resolution is 500m. there are 2 bands in a file. The first band is latitude data, and the second band is longitude data. Data type: float. Unit: degree. Projection, Sinusoidal. TIF format. One example file name of Longitude and latitude data is “latitude500m_H00V08.tif”.

Note: Data at 500m resolution are tiled according to NASA MODIS products, and the image size of each granule is 2400 rows  2400 columns.

2.	Output Data
   
The output data was GPP, NPP and LUE every 8 days.
Output data descriptions：
GPP(gross primary productivity): valid data range: 0-65535. Unit: gCm-2d-1; Projection: Sinusoidal. TIF format; The temporal resolution is 8 days. Scalefactor: 1000. 16-bit unsigned integer.
NPP(net primary productivity): valid data range: 0-65535. Unit: gCm-2d-1; Projection: Sinusoidal. TIF format. The temporal resolution is 8 days. Scalefactor: 1000. 16-bit unsigned integer.
LUE(light use efficiency): Unit: gC MJ-1; Projection: Sinusoidal. TIF format; The temporal resolution is 8 days. Type: 32-bit floating-point. 

3.	Preparation of The Program
   
Precautions before running the program:
1). Read the "Detailed Format Instructions for All Input Data" in detail.
2). Check the valid value range, scale factor, and default value of each image.
3). The spatial range of FPAR, LAI, and LANDCOVER images must be the same, and the number of rows and columns of the three images must be the same.
4). The classification codes for LANDCOVER data need to be converted to match the classification criteria of the IGBP.
5). Check that the naming format of all data is consistent with the format defined in the program. (Location: NPP.cpp module).
6). The GDAL library is required to run the program.
7). Edit the example configuration file “param_2001001.xml” (storage location: project folder) to set the path for reading data and the path for outputting files.

Here is an example of a configuration file.xml. The specific files to be configured are as follows: comments // followed by a description of each parameter.

<?xml version="1.0" encoding="UTF-8"?> 
<parameters> 
	<param name="YEAR">2001</param>//The calculation 
	<param name="N8doy">0</param>//Nth 8day of the year 
	<param name="LAI">F:\global500m_data\MUSES_LAI\2001\</param>//File storage directory of LAI data.
	<param name="FPAR">F:\global500mdata\MUSES_FPAR\2001\</param>//File storage directory of FPAR data.
	<param name="LAI_MAX">F:\global500mdata\MUSES_LAImax\2001\</param>//File storage directory of LAImaxy data.
	<param name="LANDCOVER">F:\global500mdata\MCD12Q1\2001\</param>//File storage directory of landcover data.
	<param name="ERA5">F:\global500mdata\ERA5\</param>//File storage directory of ERA5 data.
	<param name="LL_HV">F:\global500mdata\LL_HV\</param>//File storage directory or longitude and latitude data for each granule.
	<param name="SRTM_DEMH">F:\global500mdata\HV_DEMTIF\</param>//File storage directory of 500 m resolution DEM. 
	<param name="SRTM_DEML">F:\global500mdata\GlobalDEM025_1441_721.tif</param>//File storage directory of 0.25 degree resolution DEM. 
	<param name="Output">F:\global500mdata\result\2001\</param>//File storage directory of output data.
	<param name="root directory">F:\global500mdata\result\2001\</param>//File storage directory of the log file of program runing.
</parameters> 

4.	Program Execution
   
Enter the DOS environment through the CMD command, and then run the command line: “NPP.exe param_2001001.xml”
// “NPP.exe” is the executable file and “param_2001001.xml” is the configuration file.
// One run will execute all the image granule for the same time period.

