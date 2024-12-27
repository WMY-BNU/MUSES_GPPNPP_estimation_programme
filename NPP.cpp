// NPP.cpp Main program: Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "IMGALG.h"
#include "npp_fun.h"
#include <sstream>
#include "ConfigMgr.h"
#include <gdal_priv.h>
#include<io.h>
#include "Markup.h"
#include <stdlib.h>

#define CONFIG_CON 13
#define SCALE_FACT 1000


using namespace std;



//////all description added by ZHL on 2020.1.1 is for 500m GPP / NPP/ LUE

//Detailed description of all input data:
//1. Average air temperature (AAT)： Average air temperature was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: ℃. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of the AAT is “2m_temperature2001001.tif”.
//2. Surface solar radiation downwards (SSRD): SSRD was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: MJm-2d-1. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of SSRD is “surface_solar_radiation_downwards2001001.tif”.
//3. Dew point temperature: Dew point temperature was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: ℃; Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of the dew point temperature is “2m_dewpoint_temperature2001001.tif”.
//4. Surface_net_thermal_radiation (SNTR): SNTR was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: wm-2. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of SNTR is “surface_net_thermal_radiation2001001.tif”.
//5. Surface_net_solar_radiation (SNSR): SNSR was calculated from ERA5 data and has a spatial resolution of 0.25 degrees and temporal resolution of 8 days. Data type: float. Unit: wm-2. Geographic coordinate system: GCS_WGS_1984. TIF format. One example file name of SNSR is “surface_net_solar_radiation2001001.tif”.
//6. FPAR (fraction of absorbed photosynthetically active radiation): spatial resolution of the FPAR data is 500m, and the temporal resolution is 8 days. Data type: byte. Unit: none. Scale factor = 0.004. Fill value: 255. Projection: Sinusoidal. HDF format. One example file name of FPAR is “MUSES.A2001001.H05V10.MODIS.FAPAR500M.C05.HDF”.
//7. LAI (leaf and index): spatial resolution of the LAI data is 500m, and the temporal resolution is 8 days. Data type: integer. Unit: m2/m2. Scale factor = 0.01. Fill value: 2000, 2500. Projection: Sinusoidal. HDF format. One example file name of LAI is “MUSES.A2001001.H00V08.MODIS.LAI500M.C05.HDF”.
//8. LAImaxy (Maximum LAI in a year): spatial resolution is 500m. Data type: integer. Unit: m2/m2. Scale factor = 0.01. Fill value: 2000, 2500. Projection: Sinusoidal. Raw format. One example file name of LAImaxy is “LAImaxA2001H00V08.tif”.
//9. MODIS landcover: spatial resolution is 500m, Data type: Byte. Unit: none; Projection: Sinusoidal. Raw format. One example file name of landcover is “MCD12Q1.A2001001.h00v08.tif”.
//10. DEMH (Digital Elevation Model): spatial resolution is 500m, Data type: Unsigned Integer. Unit: m. Projection: Sinusoidal. TIF format. One example file name of DEMH is “DEM_500M_H00V08.TIF”.
//11. DEML (Digital Elevation Model): spatial resolution is 0.25 degrees, Data type: float. Unit: m. Geographic coordinate system: GCS_WGS_1984. TIF format. The example file name of DEML is “GlobalDEM025_1441_721.tif”.
//12. Longitude and latitude data: spatial resolution is 500m. there are 2 bands in a file. The first band is latitude data, and the second band is longitude data. Data type: float. Unit: degree. Projection, Sinusoidal. TIF format. One example file name of Longitude and latitude data is “latitude500m_H00V08.tif”.

//Output data descriptions
//GPP(gross primary productivity): valid data range: 0-65535. Unit: gCm-2d-1; Projection: Sinusoidal. TIF format; The temporal resolution is 8 days. Scalefactor: 1000. 16-bit unsigned integer.
//NPP(net primary productivity): valid data range: 0-65535. Unit: gCm-2d-1; Projection: Sinusoidal. TIF format. The temporal resolution is 8 days. Scalefactor: 1000. 16-bit unsigned integer.
//LUE(light use efficiency): Unit: gC MJ-1; Projection: Sinusoidal. TIF format; The temporal resolution is 8 days. Type: 32-bit floating-point.
///


//Precautions before running the program:
//1. Read the "Detailed Format Instructions for All Input Data" in detail
//2. Check the valid value range, scalefactor, and default value of each image.
//3. The spatial range of FPAR, LAI, AND LANDCOVER images must be the same, and the number of rows and columns of the three images must be the same

string num2str(double i)
{
	stringstream ss;
	ss << i;
	return ss.str();
}
//Read the configuration file
vector<string> readConfig(const char* path)
{
	ifstream fin(path);  
	string s; 
	vector<string> Path;
	getline(fin,s);
	getline(fin,s);
	/*getline(fin,s);*/

	for(int i=0;i<CONFIG_CON;i++)
	{
		getline(fin,s);
		getline(fin,s);
		Path.push_back(s);
	}
	return Path;
}

string getfileName(string filestring)
{
        int sstart=filestring.find_last_of('\\');
        int send=filestring.find_last_of('.');
        return filestring.substr(sstart+1,send-sstart-1);
}
bool IsRSParExist=true;
bool IsRSSMExist=true;



int _tmain(int argc, _TCHAR* argv[])
{	
	if (argc < 2)
	{
		cout << "parameters are missing" << endl;
		return 0;
	}

	//注册GDAL
	CPLSetConfigOption("CPL_DEBUG",NULL);
	GDALAllRegister();
	//CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	IsRSParExist=false;
	IsRSSMExist=true;
	
	
	//===============Read the configuration file===================================
	//for(int fd=1982001;fd<=1982353;fd+=16)
	//{
	//*****************Modify the year in which the data was run****************************************
	string fpath;
	fpath = argv[1];
	string configPath;
	configPath = fpath;
	//configPath = "F:\\ZHL\\global_GPP_500m\\NPP_CONFIG_TEST.txt";
	CConfigMgr clsConfig;
	if (!clsConfig.Parse(configPath.c_str()))
	{
		cout << "Failed to interpret the configuration file." << endl;
		return 0;
	}//add by zhl interpret the xml configuration file
	vector<string> Path = readConfig(configPath.c_str());
	string  yearnum = clsConfig.GetPathBy("YEAR");
	int year = atoi(yearnum.c_str());//zhl

	for (int index_year = year; index_year <= year; index_year++)
	{
		string  d_numvalue = clsConfig.GetPathBy("N8doy");
		int start_num = atoi(d_numvalue.c_str());//zhl
		for (int d = start_num; d <= start_num; d++)//d<46
		{
			int year = index_year;
			char startDaychar[10];
			itoa(year, startDaychar, 10);
			string stryear(startDaychar);
			int FirstDay = 1000 * year + 1;//fd;
		/* ****************************************************************** */

			for (int H_index = 0;H_index<=35;H_index++)
			{
				for(int V_index=0;V_index<18;V_index++)
				{
					string H;
					string V;
					if(H_index<10)
					{
						H ="0"+to_string(H_index);
					}
					else
					{
						H = to_string(H_index);
					}
					if(V_index<10)
					{
						V ="0"+to_string(V_index);
					}
					else
					{
						V = to_string(V_index);
					}
					string HV = "H"+H+"V"+V;
					string HV2= "h"+H+"v"+V;


					//******************Modify the configuration file path*****************************************************

					string root = clsConfig.GetPathBy("Output");
					const char* logoPath = root.append("NPP_logo.txt").c_str();
					makeLogo(logoPath, "Program running");
					string temp = "Configuration file：" + configPath;
					makeLogo(logoPath, temp.c_str());
					//string temp_dir = Path[5];
					string temp_dir = clsConfig.GetPathBy("Output");




				/* ****************************************************************** */
	
				//The date on which the output file name is obtained, which is also used to find the name of the file in the folder that needs to be preprocessed.

				//string lai = Path[0];//500m resolution LAI file path + file name
				//string fpar = Path[1];//500m resolution FPAR file path + file name
				string lai = clsConfig.GetPathBy("LAI");
				string fpar = clsConfig.GetPathBy("FPAR");
				string laimax = clsConfig.GetPathBy("LAI_MAX");
				int startDayint=FirstDay+d*8;
				itoa(startDayint,startDaychar,10);
				string startDay(startDaychar);

				ifstream fin(laimax + "LAImaxA" + to_string(year) + HV + ".tif");//Check whether there is LAImax data of the row and column numbers and the date zhl20200926
				if (!fin)
				{
					cout << "Can not open LAImax" << endl;
					continue;
				}
				ifstream fin1(fpar + "MUSES.A" + startDay + "." + HV + ".MODIS.FAPAR500M.C05.HDF");//Check whether there is FPAR data of the row and column numbers and the date zhl20200926
				if (!fin1)
				{
					cout << "Can not open fpar" << endl;
					continue;
				}
				ifstream fin2(lai + "MUSES.A" + startDay + "." + HV + ".MODIS.LAI500M.C05.HDF");//Check whether there is LAI data of the row and column numbers and the date zhl20200926
				if (!fin2)
				{
					cout << "Can not open LAI" << endl;
					continue;
				}

				string lai300m = lai + "MUSES.A" + startDay + "." + HV + ".MODIS.LAI500M.C05.HDF";
				string fpar1000m = fpar + "MUSES.A" + startDay + "." + HV + ".MODIS.FAPAR500M.C05.HDF";

					cout << '\r' <<"The start time to be calculated: " <<startDay;
					cout << '\n';

					const char* gldasSourcePath = clsConfig.GetPathBy("ERA5").c_str();//Meteorological data input path
					string gldasResultPath = clsConfig.GetPathBy("Output");//Meteorological data output path

					////========preprocessing remote sensing inversion PAR and soil moisture data, HDF format=======

					int PARnodata = 0;   //Determined based on the filled value of the data

					float FPARnodata = 255;  //Determined based on the filled value of the data

					float LAInodata = 1000;//25.5;//Determined based on the filled value of the data //Determined based on the filled value of the data, GLASS 500M LAI valid range=[0,1000],it is changed to greater than 250 if invalid in the core_fun.cpp

					int landcovernodata = 255;  //Determined based on the filled value of the data
					////========Generate a latitude and longitude file to generate a corresponding latitude and longitude image based on a 500-meter landcover===================
	
					string landcover1000m = clsConfig.GetPathBy("LANDCOVER") + "MCD12Q1.A" + to_string(year) + "001." + HV2 + ".tif";
					string latitudePathFile = clsConfig.GetPathBy("LL_HV") + "latitude500m_" + HV + ".tif";


			
					//===============Open the file with GDAL===================================
					GDALDataset  *FPAR_DT = (GDALDataset*)GDALOpen(fpar1000m.c_str(), GA_ReadOnly);
					if(FPAR_DT==NULL)
					{
						cout << "Can not open: "<<fpar1000m<< endl;
						temp="Can not open file: "+ fpar1000m;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}

					//LAI(Leaf Area Index 500m)tiff///added by MJ Wang 2019.3.28
					GDALDataset  *LAI_DT = (GDALDataset*)GDALOpen(lai300m.c_str(), GA_ReadOnly);
					if(LAI_DT==NULL) 
					{
						cout << "Can not open: "<<lai300m<< endl;
						temp="Can not open file: "+ lai300m;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}

					//LAImaxy 500M TIFF  //Since the global quantity of modis 500m is inconsistent with the LAI fpar, it is based on LAI's, and if there is no document or the HV number is invalid, it will go directly to the next cycle
					string LAIMAXY500 = clsConfig.GetPathBy("LAI_MAX") + "LAImaxA" + to_string(year) + HV + ".tif";
					GDALDataset  *LAImaxy_DT = (GDALDataset*)GDALOpen(LAIMAXY500.c_str(), GA_ReadOnly);
					if (LAImaxy_DT == NULL)
					{
						cout << "Can not open: " << LAIMAXY500 << endl;
						temp = "Can not open file: " + LAIMAXY500;
						makeLogo(logoPath, temp.c_str());
						makeLogo(logoPath, "The programme does not successfully run");
						return 4;

					}

					GDALDataset *LANDCOVER_DT=(GDALDataset*)GDALOpen(landcover1000m.c_str(), GA_ReadOnly);
					if(LANDCOVER_DT==NULL)
					{
						cout <<"Can not open: "<<landcover1000m<< endl;
						temp="Can not open file: "+ landcover1000m;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}

					string gldasAatPathfile = clsConfig.GetPathBy("ERA5") + "2m_temperature" + startDay + ".tif"; //0.25D ERA5 ZHL
					GDALDataset *GLDAS_Tav_DT=(GDALDataset*)GDALOpen(gldasAatPathfile.c_str(), GA_ReadOnly);
					if(GLDAS_Tav_DT==NULL)
					{
						cout <<"Can not open: "<<gldasAatPathfile<< endl;
						temp="Can not open file: "+ gldasAatPathfile;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}
					string gldasParPathfile = clsConfig.GetPathBy("ERA5") + "surface_solar_radiation_downwards" + startDay + ".tif"; //0.25D ERA5 ZHL
					GDALDataset *GLDAS_par_DT=(GDALDataset*)GDALOpen(gldasParPathfile.c_str(), GA_ReadOnly);
					if(GLDAS_par_DT==NULL)
					{
						cout << "Can not open: "<<gldasParPathfile<< endl;
						temp="Can not open file: "+ gldasParPathfile;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}


					string gldasSWnetPathfile = clsConfig.GetPathBy("ERA5") + "surface_net_solar_radiation" + startDay + ".tif"; //0.25D ERA5 ZHL
					GDALDataset *GLDAS_swnet_DT=(GDALDataset*)GDALOpen(gldasSWnetPathfile.c_str(), GA_ReadOnly);
					if(GLDAS_swnet_DT==NULL)
					{
						cout << "Can not open: "<<gldasSWnetPathfile<< endl;
						temp="Can not open file: "+ gldasSWnetPathfile;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}


					string gldasLWnetPathfile = clsConfig.GetPathBy("ERA5") + "surface_net_thermal_radiation" + startDay + ".tif"; //0.25D ERA5 ZHL
					GDALDataset *GLDAS_lwnet_DT=(GDALDataset*)GDALOpen(gldasLWnetPathfile.c_str(), GA_ReadOnly);
					if(GLDAS_lwnet_DT==NULL)
					{
						cout << "Can not open: "<<gldasLWnetPathfile<< endl;
						temp="Can not open file: "+ gldasLWnetPathfile;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}

					string gldasDTPathfile = clsConfig.GetPathBy("ERA5") + "2m_dewpoint_temperature" + startDay + ".tif"; //0.25D ERA5 ZHL
					GDALDataset *GLDAS_dt_DT=(GDALDataset*)GDALOpen(gldasDTPathfile.c_str(), GA_ReadOnly);
					if(GLDAS_dt_DT==NULL)
					{
						cout << "Can not open: "<<gldasDTPathfile<< endl;
						temp="Can not open file: "+ gldasDTPathfile;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}
					string DEM500 = clsConfig.GetPathBy("SRTM_DEMH") + "DEM_500M_" + HV + ".TIF";
					GDALDataset *DEMM_DT = (GDALDataset*)GDALOpen(DEM500.c_str(), GA_ReadOnly);
					if(DEMM_DT==NULL)
					{
						cout <<"Can not open: "<<DEM500<< endl;
						temp="Can not open file: "+DEM500;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}

					//SRTM DEML(Digital Elevation Model 0.25degree) TIFF  ////added by MJ Wang 2019.3.28
					GDALDataset  *DEML_DT = (GDALDataset*)GDALOpen(clsConfig.GetPathBy("SRTM_DEML").c_str(), GA_ReadOnly);
					if(DEML_DT==NULL)
					{
						cout << "Can not open: " << clsConfig.GetPathBy("SRTM_DEML") << endl;
						temp = "Can not open file: " + clsConfig.GetPathBy("SRTM_DEML");
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}

					string latPathFile = clsConfig.GetPathBy("LL_HV") + "latitude500m_" + HV + ".tif";//500 m latitude and longitude file path + file name
					GDALDataset  *Lati_DT = (GDALDataset*)GDALOpen(latPathFile.c_str(), GA_ReadOnly);
					if(Lati_DT==NULL)
					{
						cout <<"Can not open: "<<latPathFile<< endl;
						temp="Can not open file: "+ latPathFile;
						makeLogo(logoPath,temp.c_str());
						makeLogo(logoPath,"The programme does not successfully run");
						return 4;
					}
			

					string root_Dir = clsConfig.GetPathBy("root directory");//Output file path

					//===============Read spatial coordinates and projection information===================================
					//GLDAS_GeoTransform[6]  The array GLDAS_GeoTransform holds some parameters from the affine transformation, which are explained below
					//GLDAS_GeoTransform[0]  X coordinates in the upper left corner
					//GLDAS_GeoTransform[1]  East-west resolution
					//GLDAS_GeoTransform[2]  Rotation angle, 0 indicates the image "North up"
					//GLDAS_GeoTransform[3]  Y coordinates in the upper left corner
					//GLDAS_GeoTransform[4]  Rotation angle, 0 indicates the image "North up"
					//GLDAS_GeoTransform[5]  North-south resolution
					double GLDAS_GeoTransform[6]; //0.25 degree, meteorological data  //added by MJ Wang 2019.3.28
					double STRMM_DEM_GeoTransform[6]; //500m DEM //added by MJ Wang 2019.3.28
					double Latitude_GeoTransform[6]; //500m latitude//added by MJ Wang 2019.3.28
					double DEML_GeoTransform[6];//0.25°DEM   //added by MJ Wang 2019.3.28

					///average air temperature
					GLDAS_Tav_DT->GetGeoTransform(GLDAS_GeoTransform);
					int x=GLDAS_GeoTransform[0];//X coordinates in the upper left corner
					//ERA5 data 0.25 degree
					int GLDAS_brc=GLDAS_Tav_DT->GetRasterCount();///Get the number of bands
					int GLDAS_xrc=GLDAS_Tav_DT->GetRasterXSize();///X-axis image size
					int GLDAS_yrc=GLDAS_Tav_DT->GetRasterYSize();///Y-axis image size
					/////DEM 500m
					DEMM_DT->GetGeoTransform(STRMM_DEM_GeoTransform);
			
					int DEMM_brc=DEMM_DT->GetRasterCount();
					int DEMM_xrc=DEMM_DT->GetRasterXSize();
					int DEMM_yrc=DEMM_DT->GetRasterYSize();
					/////500m 
					Lati_DT->GetGeoTransform(Latitude_GeoTransform);
					const char* projectionref=Lati_DT->GetProjectionRef();//The coordinate system of the latitude and longitude image
					//500m，input data, the range needs to be calculated输入的数据，需要计算的范围大小
					int Lati_brc=Lati_DT->GetRasterCount();
					int Lati_xrc=Lati_DT->GetRasterXSize();
					int Lati_yrc=Lati_DT->GetRasterYSize();
					////DEM 0.25 degrees of coarse spatial resolution
					DEML_DT->GetGeoTransform(DEML_GeoTransform);
					int DEML_xrc=DEML_DT->GetRasterXSize();
					int DEML_yrc=DEML_DT->GetRasterYSize();

	

					//===============Declares pointers for various stored images===================================
					//The naming method is the same as GDALDataset+paf, H stands for data with a resolution of 500M, M stands for data with a resolution of 500M,
					//ML stands for data with a resolution of 500m, and L stands for data with a resolution of 0.25 degrees


					//Define SWnet pointer with a resolution of 0.25 degrees// revised by MJ Wang 2019.4.19
					float* pafSWnetL=new float[sizeof(float)*GLDAS_xrc*GLDAS_yrc];//zhl

					//Define LWnet pointer with a resolution of 0.25 degrees//revised by MJ Wang 2019.4.19
					float* pafLWnetL=new float[sizeof(float)*GLDAS_xrc*GLDAS_yrc];//zhl
					//Define DTt pointer with a resolution of 0.25 degrees//revised by MJ Wang 20190812
					float* pafDTL=new float[sizeof(float)*GLDAS_xrc*GLDAS_yrc];//zhl
					//Define PAR pointer with a resolution of 0.25 degrees// revised by MJ Wang 2019.4.19
					float* pafGldasParL=new float[sizeof(float)*GLDAS_xrc*GLDAS_yrc];//zhl
					//Define AAT pointer with a resolution of 0.25 degrees (after DEM correction)// revised by MJ Wang 2019.4.19
					float* pafTavL=new float[sizeof(float)*GLDAS_xrc*GLDAS_yrc];//zhl
					//Define AAT pointer with a resolution of 0.25 degrees (before DEM correction)// revised by MJ Wang 2019.4.19
					float* pafTavGldasL=new float[sizeof(float)*GLDAS_xrc*GLDAS_yrc];//zhl
					//Define DEM pointer with a resolution of 0.25 degrees//revised by MJ Wang 2019.4.19
					float* pafDemL=new float[sizeof(float)*DEML_xrc*DEML_yrc];//zhl
					//Define DEM pointer with a resolution of 500m// revised by MJ Wang 2019.4.19
					float* pafDemM=new float[sizeof(float)*DEMM_xrc];//zhl
					//Define the longitude pointer (500m)
					float* pafLongH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the latitude pointer (500m)//revised by MJ Wang 2019.4.19
					float* pafLatH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the LUE pointer (500m)
					float* pafLUEH=new float[sizeof(float)*Lati_xrc];//zhl
					//Definition of average air temperature pointer (after DEM correction)（500m）revised by MJ Wang 2019.4.19
					float* pafTavH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the average air temperature pointer (before DEM correction)（500m）revised by MJ Wang 2019.4.19
					float* pafTavGldasH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the optimal temperature pointer（500m）revised by MJ Wang 2019.4.19
					float* pafToptH=new float[sizeof(float)*Lati_xrc];//zhl
					//Defines the pointer of the parameter k1（500m）revised by MJ Wang 2019.4.19
					double* pafK1H=new double[sizeof(double)*Lati_xrc];//zhl
					//Defines the pointer of the parameter k2（500m）revised by MJ Wang 2019.4.19
					double* pafK2H=new double[sizeof(double)*Lati_xrc];//zhl
					//Defines the pointer of the parameter TE2（500m）revised by MJ Wang 2019.4.19
					double* pafTE2H=new double[sizeof(double)*Lati_xrc];//zhl

					////////////////////////////////////////////////////////////

					//Define the pointer of the parameter Evaporation（500m）revised by MJ Wang 2019.4.19
					double* pafEvaporation=new double[sizeof(double)*Lati_xrc];//zhl
					//Define the pointer of the parameter Transporation（500m）
					double* pafTransporation=new double[sizeof(double)*Lati_xrc];//zhl
					//Define the pointer of the parameter Evapotransporation（500m）
					double* pafEvapotransporation=new double[sizeof(double)*Lati_xrc];//zhl
					//Define the pointer of the parameter Potential Evapotransporation（500m）
					double* pafPEvapotransporation=new double[sizeof(double)*Lati_xrc];//zhl

					//Define the pointer of the parameter f2（500m）
					float* paff2H=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the pointer of the parameter CI（500m）added by MJ Wang 2019.4.20
					float* CI=new float[sizeof(float)*Lati_xrc];//zhl
					//Defines the pointer of PAR after resampling（500m）revised by MJ Wang 2019.4.19
					float* pafGldasParH=new float[sizeof(float)*Lati_xrc];//zhl
					//Defines the pointer of SWnet after resampling（500m）revised by MJ Wang 2019.4.19
					float* pafSWnetH=new float[sizeof(float)*Lati_xrc];//zhl
					//Defines the pointer of LWnet after resampling（500m）revised by MJ Wang 2019.4.19
					float* pafLWnetH=new float[sizeof(float)*Lati_xrc];//zhl
					//Defines the pointer of DT after resampling（500m）revised by MJ Wang 20190812
					float* pafDTH=new float[sizeof(float)*Lati_xrc];//zhl
					//Defines the pointer of the remote sensing inverted PAR (500m)
					float* pafRsParH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the pointer to the PAR required to calculate the GPP(500m)revised by MJ Wang 2019.4.19
					float* pafParH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the pointer to the FPAR required to calculate the GPP(500m)revised by MJ Wang 2019.4.19
					DT_8U* pafFparH = new DT_8U[sizeof(DT_8U)*Lati_xrc];//zhl
					//Define the pointer to the LAI required to calculate the GPP(500m)revised by MJ Wang 2019.4.19
					DT_16U* pafLaiH = new DT_16U[sizeof(DT_16U)*Lati_xrc];//zhl
					//Define the pointer to the LANDCOVER required to calculate the GPP(500m)revised by MJ Wang 2019.4.19
					DT_8U* pafLandcoverH=new DT_8U[sizeof(DT_8U)*Lati_xrc];//zhl
					//Define the GPP pointer(500m)
					double* pafGppH=new double[sizeof(double)*Lati_xrc];//zhl
					//Define the output LUE pointer (500m)revised by MJ Wang 2019.6.1
					float* pafLUEoutp=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the pointer for calculating the maintenance respiration RM(500m)revised by MJ Wang 2019.4.19
					float* pafRmH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the NPP pointer(500m)revised by MJ Wang 2019.4.19
					double* pafNppH=new double[sizeof(double)*Lati_xrc];//zhl
					//Define the DEM pointer(500m)revised by MJ Wang 2019.4.19
					float* pafDEMH=new float[sizeof(float)*Lati_xrc];//zhl
					//Define the LAImaxy pointer(500m)revised by MJ Wang 2019.4.19
					double* pafLaimaxH = new double[sizeof(double)*Lati_xrc];//zhl
					
					//Define the pointer to output GPP, NPP (16-bit unsigned integer)(500m)revised by MJ Wang 2019.4.19
					DT_16U* pafGPPShortIntH=new DT_16U[sizeof(DT_16U)*Lati_xrc];//zhl
					DT_16U* pafNPPShortIntH=new DT_16U[sizeof(DT_16U)*Lati_xrc];//zhl


					//===============Initialize the pointer===================================
					initPoint(pafGldasParL,GLDAS_xrc*GLDAS_yrc,0);
					initPoint(pafTavL,GLDAS_xrc*GLDAS_yrc,0);
					initPoint(pafTavGldasL,GLDAS_xrc*GLDAS_yrc,0);
					initPoint(pafDemL,DEML_xrc*DEML_yrc,0);
					initPoint(pafSWnetL,GLDAS_xrc*GLDAS_yrc,0);
					initPoint(pafLWnetL,GLDAS_xrc*GLDAS_yrc,0);
					initPoint(pafDTL,GLDAS_xrc*GLDAS_yrc,0);//added by MJ Wang 20190812
					initPoint(pafLongH,Lati_xrc,0);
					initPoint(pafLatH,Lati_xrc,0);
					initPoint(pafLUEH,Lati_xrc,0);
					initPoint(pafTavH,Lati_xrc,0);
					initPoint(pafTavGldasH,Lati_xrc,0);
					initPoint(pafToptH,Lati_xrc,0);
					initPoint(pafK1H,Lati_xrc,0);
					initPoint(pafK2H,Lati_xrc,0);
					initPoint(pafTE2H,Lati_xrc,0);
					initPoint(paff2H,Lati_xrc,0);
					initPoint(CI,Lati_xrc,0);
					initPoint(pafGldasParH,Lati_xrc,0);
					initPoint(pafSWnetH,Lati_xrc,0);
					initPoint(pafLWnetH,Lati_xrc,0);
					initPoint(pafDTH,Lati_xrc,0);///added by MJ Wang 20190812
					initPoint(pafDemM, DEMM_xrc, 0);
					initPoint(pafRsParH,Lati_xrc,0);
					initPoint(pafParH,Lati_xrc,0);
					initPoint(pafFparH,Lati_xrc,0);
					initPoint(pafLaiH,Lati_xrc,0);
					initPoint(pafLandcoverH,Lati_xrc,0);
					initPoint(pafGppH,Lati_xrc,0);
					initPoint(pafLUEoutp,Lati_xrc,0);
					initPoint(pafRmH,Lati_xrc,0);
					initPoint(pafNppH,Lati_xrc,0);
					initPoint(pafLaimaxH,Lati_xrc,0);
					initPoint(pafDemM, DEMM_xrc, 0);
					initPoint(pafGPPShortIntH,Lati_xrc,0);
					initPoint(pafNPPShortIntH,Lati_xrc,0);
					initPoint(pafDEMH,Lati_xrc,0);
					initPoint(pafEvaporation,Lati_xrc,0);
					initPoint(pafTransporation,Lati_xrc,0);
					initPoint(pafEvapotransporation,Lati_xrc,0);
					initPoint(pafPEvapotransporation,Lati_xrc,0);

					//============According to the latitude and longitude, the coordinates of 4 points were found from the image of 0.25° (500m).=====================
					//=============Because the coordinates of the center point of the cell are required, they are processed on latitude and longitude=====================

					//The latitude and longitude coordinates and pixel size of the pixel center of the upper left point of the ERA5S image (0.25 degrees).
					//ERA5 data (0.25 degrees) Pixel size/latitude and longitude coordinates at the center of the pixel at the upper left point
					float GLDAS_pix_size = GLDAS_GeoTransform[1];
					float GLDAS_Lup_long = GLDAS_GeoTransform[0] + GLDAS_pix_size/2;
					float GLDAS_Lup_lat = GLDAS_GeoTransform[3] - GLDAS_pix_size/2;

					//The data prepared by yourself (500m) is the latitude and longitude coordinates of the center of the pixel at the upper left point, and the size of the cell
					//High-resolution DEM data (500m) with pixel size/latitude and longitude coordinates at the center of the upper left point pixel
					float STRMM_DEM_pix_size = STRMM_DEM_GeoTransform[1];
					float STRMM_DEM_Lup_long = STRMM_DEM_GeoTransform[0] + STRMM_DEM_pix_size/2;
					float STRMM_DEM_Lup_lat = STRMM_DEM_GeoTransform[3] - STRMM_DEM_pix_size/2;	

					//The latitude and longitude coordinates and pixel size of the pixel center of the upper left point of 500m data after preprocessing
					//The latitude and longitude coordinates and pixel size of the pixel center of the upper left point of 500m data
					float HResolution_pix_size=Latitude_GeoTransform[1];
					float HResolution_Lup_long=Latitude_GeoTransform[0]+HResolution_pix_size/2;
					float HResolution_Lup_Lat=Latitude_GeoTransform[3]-HResolution_pix_size/2;

					////============Output files GPP, NPP=====================================
					//The final result is saved in TIF format and imported into the result line by line
					GDALDataset *gppDataset;//Defines the generated GPP file pointer
					GDALDataset *nppDataset;//Defines the generated NPP file pointer
					GDALDataset *lueDataset;

					//The driver used to create a new file
					GDALDriver *poDriver;
					//Create a new file
					const char *fomat="GTiff";
					poDriver = GetGDALDriverManager()->GetDriverByName(fomat);
					//Get the file format type
					//char **papszMetadata = poDriver->GetMetadata();//gggggggggggggggggggggggggggggggggggggg
					char **papszMetadata = NULL;

					string gppsavePath_temp = temp_dir +"GPP.500m."+ startDay +HV+".V01.tif" ;
					string nppsavePath_temp = temp_dir +"NPP.500m."+ startDay+ HV+".V01.tif";
					string luesavePath_temp = temp_dir +"LUE.500m."+ startDay+ HV+".V01.tif";

					const char* gppsavePath = gppsavePath_temp.c_str();//The path + file name of the GPP is generated
					const char* nppsavePath = nppsavePath_temp.c_str();//The path + file name of the NPP is generated
					const char* luesavePath = luesavePath_temp.c_str();//The path + file name of the LUE is generated

					////write GPP data
					gppDataset=poDriver->Create(gppsavePath,Lati_xrc,Lati_yrc,1,GDT_Int16,papszMetadata);
					if(gppDataset==NULL)
					{
						cout<<"Dataset creation failed!"<<endl;
						Clean_Dir(temp_dir.c_str(),0,".");
						return 9;
					}
					gppDataset->SetGeoTransform(Latitude_GeoTransform);
					gppDataset->SetProjection(projectionref);
					//write NPP data
					nppDataset=poDriver->Create(nppsavePath,Lati_xrc,Lati_yrc,1,GDT_Int16,papszMetadata);
					if(nppDataset==NULL)
					{
						cout<<"Dataset creation failed!"<<endl;
						Clean_Dir(temp_dir.c_str(),0,".");
						return 9;
					}
					nppDataset->SetGeoTransform(Latitude_GeoTransform);
					nppDataset->SetProjection(projectionref);

					//write LUE data
			
					lueDataset=poDriver->Create(luesavePath,Lati_xrc,Lati_yrc,1,GDT_Float32,papszMetadata);
					if(lueDataset==NULL)
					{
						cout<<"Dataset creation failed!"<<endl;
						Clean_Dir(temp_dir.c_str(),0,".");
						return 9;
					}
					lueDataset->SetGeoTransform(Latitude_GeoTransform);
					lueDataset->SetProjection(projectionref);
		

					////============Outside the loop, the entire 500m resolution of the entire frame is read, reducing the calculation time=====================================

					//Read the incident shortwave radiation downward of ERA5 (0.25 degrees)//revised by MJ Wang 2019.4.19
					GDALRasterBand* rasterBand=GLDAS_par_DT->GetRasterBand(1);
					GDALDataType dataType =rasterBand->GetRasterDataType();
					rasterBand->RasterIO(GF_Read,0,0,GLDAS_xrc,GLDAS_yrc,pafGldasParL,GLDAS_xrc,GLDAS_yrc,dataType,0,0);

					
					//Read the net shortwave radiation LWNET of ERA5 (0.25 degrees)//revised by MJ Wang 2019.4.19
					rasterBand=GLDAS_swnet_DT->GetRasterBand(1);
					dataType =rasterBand->GetRasterDataType();
					rasterBand->RasterIO(GF_Read,0,0,GLDAS_xrc,GLDAS_yrc,pafSWnetL,GLDAS_xrc,GLDAS_yrc,dataType,0,0);
					
					//Read the net longwave radiation LWNET of ERA5 (0.25 degrees)//revised by MJ Wang 2019.4.19
					rasterBand=GLDAS_lwnet_DT->GetRasterBand(1);
					dataType =rasterBand->GetRasterDataType();
					rasterBand->RasterIO(GF_Read,0,0,GLDAS_xrc,GLDAS_yrc,pafLWnetL,GLDAS_xrc,GLDAS_yrc,dataType,0,0);


					////added by MJ Wang 20190812
					//read Dewpoint Temperature（0.25 degree）//revised by MJ Wang 20190812
					rasterBand=GLDAS_dt_DT->GetRasterBand(1);
					dataType =rasterBand->GetRasterDataType();
					rasterBand->RasterIO(GF_Read,0,0,GLDAS_xrc,GLDAS_yrc,pafDTL,GLDAS_xrc,GLDAS_yrc,dataType,0,0);

					
					//Read the average air temperature of the EVA (0.25 degrees)//revised by MJ Wang 2019.4.19
					rasterBand=GLDAS_Tav_DT->GetRasterBand(1);
					dataType =rasterBand->GetRasterDataType();
					rasterBand->RasterIO(GF_Read,0,0,GLDAS_xrc,GLDAS_yrc,pafTavL,GLDAS_xrc,GLDAS_yrc,dataType,0,0);
					rasterBand->RasterIO(GF_Read,0,0,GLDAS_xrc,GLDAS_yrc,pafTavGldasL,GLDAS_xrc,GLDAS_yrc,dataType,0,0);

					
					//Read low-resolution DEML (0.25 degrees)//revised by MJ Wang 2019.4.19
					rasterBand=DEML_DT->GetRasterBand(1);
					dataType =rasterBand->GetRasterDataType();
					rasterBand->RasterIO(GF_Read,0,0,DEML_xrc,DEML_yrc,pafDemL,DEML_xrc,DEML_yrc,dataType,0,0);

					////deleted by MJ Wang 2019.4.14
					//The average EVA temperature is corrected to the sea level according to the geographical coordinates
					float DEML_pix_size=DEML_GeoTransform[1];
					float DEML_Lup_Long=DEML_GeoTransform[0]+DEML_pix_size/2;
					float DEML_Lup_Lat=DEML_GeoTransform[3]-DEML_pix_size/2;
					float Long_dis2=GLDAS_Lup_long-DEML_Lup_Long;
					float Lat_dis2=DEML_Lup_Lat-GLDAS_Lup_lat;
					int x_Low2=int(Long_dis2/GLDAS_pix_size);
					int y_Low2=int(Lat_dis2/GLDAS_pix_size);

					for (int j = 0; j < GLDAS_yrc; j++)
					{
						for (int i = 0; i < GLDAS_xrc; i++)
						{
							float dt=pafDemL[(y_Low2+j)*DEML_xrc+x_Low2+i]*AIR_COR;
					
							pafTavL[j*GLDAS_xrc+i]=pafTavL[j*GLDAS_xrc+i]+dt;
					
						}
					}
			
								
					//There are cases where PAR data for remote sensing inversion exists/////
					if(IsRSParExist)
					{
						//Remote Sensing PAR(5km)TIFF
						string rsParPathfile =  Path[5] + "rs_par"+HV+".tif"; //5 km PAR file path + file name for Remote Sensing inversion
						GDALDataset *RS_par_DT=(GDALDataset*)GDALOpen(rsParPathfile.c_str(), GA_ReadOnly);
						if(RS_par_DT==NULL)
						{
							cout <<"Can not open: "<<rsParPathfile<< endl;
							temp="Can not open file: "+ rsParPathfile;
							makeLogo(logoPath,temp.c_str());
							makeLogo(logoPath,"The programme does not successfully run");
							return 4;
						}

						double RsPar_GeoTransform[6]; //5 km, PAR for remote sensing inversion
						RS_par_DT->GetGeoTransform(RsPar_GeoTransform);

						//Remote sensing inversion of the PAR 5 km
						int RsPar_brc=RS_par_DT->GetRasterCount();
						int RsPar_xrc=RS_par_DT->GetRasterXSize();
						int RsPar_yrc=RS_par_DT->GetRasterYSize();


						//Defines the PAR pointer for remote sensing inversion with a resolution of 5 km
						//float* pafRsParML=(float*)CPLMalloc(sizeof(float)*RsPar_xrc*RsPar_yrc);
						float* pafRsParML=new float[sizeof(float)*RsPar_xrc*RsPar_yrc];//zhl
						initPoint(pafRsParML,RsPar_xrc*RsPar_yrc,0);

						//The latitude and longitude coordinates and pixel size of the center of the pixel at the upper left point of the PAR image (5km) inverted by remote sensing
						float RsPar_pix_size = RsPar_GeoTransform[1];
						float RsPar_Lup_long = RsPar_GeoTransform[0] + RsPar_pix_size/2;
						float RsPar_Lup_lat = RsPar_GeoTransform[3] - RsPar_pix_size/2;

						//Read PAR for remote sensing inversion (5 km)
						rasterBand=RS_par_DT->GetRasterBand(1);
						dataType =rasterBand->GetRasterDataType();
						rasterBand->RasterIO(GF_Read,0,0,RsPar_xrc,RsPar_yrc,pafRsParML,RsPar_xrc,RsPar_yrc,dataType,0,0);
		

						for(int j=0;j<Lati_yrc;j++)
						{
							comPar(Lati_DT, j, IsRSParExist, Lati_xrc, GLDAS_xrc, RsPar_xrc, GLDAS_GeoTransform, RsPar_GeoTransform,
								pafLongH, pafLatH, pafGldasParL, pafGldasParH, pafRsParML, pafRsParH, pafParH, PARnodata);

							////compute clearness index CI
							comCI(Lati_DT, j, d,Lati_xrc, GLDAS_xrc, GLDAS_GeoTransform,pafLatH, pafParH,CI);

							//Calculate the maximum light use efficiency based on the type of land use
							computerLUE(LANDCOVER_DT, j, Lati_xrc, pafLandcoverH, pafLUEH, CI, FPAR_DT, LAI_DT, pafFparH, pafLaiH);

							//Calculate the temperature stress factor
							comTempStress(Lati_DT,DEMM_DT, j, Lati_xrc, GLDAS_xrc,pafLandcoverH, DEMM_xrc, GLDAS_GeoTransform, STRMM_DEM_GeoTransform,
											pafLongH,  pafLatH,  pafTavH,  pafTavL, pafToptH,
											pafDemM,  pafDEMH,  pafK1H,  pafK2H,  pafTE2H);

							comET(Lati_DT, j, Lati_xrc, GLDAS_xrc, GLDAS_GeoTransform, pafLongH, pafLatH,pafTavGldasH,pafTavGldasL ,pafSWnetH, pafSWnetL,pafLWnetH,pafLWnetL,pafDTH,pafDTL,
								pafDEMH,FPAR_DT,LAI_DT,LANDCOVER_DT,pafFparH,pafLaiH,pafLandcoverH,pafEvaporation,pafTransporation,pafEvapotransporation,pafPEvapotransporation,paff2H,pafParH);

							comGPP(FPAR_DT, j, Lati_xrc, pafFparH, FPARnodata, pafLandcoverH, pafGppH,pafLUEoutp, pafLUEH, paff2H, pafTE2H, pafParH);//reviosed by MJ Wang 20190601
				
							comNPP(LAI_DT, LAImaxy_DT,j, Lati_xrc, pafLaiH, pafLaimaxH, pafLandcoverH, pafRmH, pafTavH, pafGppH, pafNppH, LAInodata);
					
							for(int k = 0; k < Lati_xrc; k++)
							{
								pafGPPShortIntH[k] = (DT_16U)(pafGppH[k] * SCALE_FACT);

								pafNPPShortIntH[k] = (DT_16U)(pafNppH[k] * SCALE_FACT);

							}

							//Write one row of the processed GPP image to the dataset
							gppDataset->RasterIO(GF_Write,0,j,Lati_xrc,1,pafGPPShortIntH,Lati_xrc,1,GDT_Int16,1,0,0,0,0);
							//Write one row of the processed NPP image to the dataset
							nppDataset->RasterIO(GF_Write,0,j,Lati_xrc,1,pafNPPShortIntH,Lati_xrc,1,GDT_Int16,1,0,0,0,0);
							//Write one row of the processed LUE image to the dataset
							lueDataset->RasterIO(GF_Write,0,j,Lati_xrc,1,pafLUEoutp,Lati_xrc,1,GDT_Float32,1,0,0,0,0);
							cout << flush << '\r' <<"The calculations completed: "+HV <<(j+1)*100/Lati_yrc<<'%';

						}
						GDALClose(gppDataset);
						GDALClose(nppDataset);
						GDALClose(lueDataset);
						delete[] pafRsParML;
						pafRsParML=NULL;//zhl
					}
					else
					{
						int RsPar_xrc=0;
						double RsPar_GeoTransform[6];
						float* pafRsParML = new float[sizeof(float)];//zhl
		
		
						for(int j=0;j<Lati_yrc;j++)
						{
							comPar(Lati_DT, j, IsRSParExist, Lati_xrc, GLDAS_xrc, RsPar_xrc, GLDAS_GeoTransform, RsPar_GeoTransform,
								pafLongH, pafLatH, pafGldasParL, pafGldasParH, pafRsParML, pafRsParH, pafParH, PARnodata);
							////compute clearness index CI
							comCI(Lati_DT, j, d,Lati_xrc, GLDAS_xrc, GLDAS_GeoTransform, pafLatH, pafParH,CI);
					
							//Calculate the maximum light use efficiency based on the type of land use
							computerLUE(LANDCOVER_DT, j, Lati_xrc, pafLandcoverH, pafLUEH, CI, FPAR_DT, LAI_DT, pafFparH, pafLaiH);

							//Calculate the temperature stress factor
							comTempStress(Lati_DT,DEMM_DT, j, Lati_xrc, GLDAS_xrc,pafLandcoverH, DEMM_xrc, GLDAS_GeoTransform, STRMM_DEM_GeoTransform,
											pafLongH,  pafLatH,  pafTavH,  pafTavL, pafToptH,
											pafDemM,  pafDEMH,  pafK1H,  pafK2H,  pafTE2H);

							comET(Lati_DT, j, Lati_xrc, GLDAS_xrc, GLDAS_GeoTransform, pafLongH, pafLatH,pafTavGldasH,pafTavGldasL ,pafSWnetH, pafSWnetL,pafLWnetH,pafLWnetL,pafDTH,pafDTL,
								pafDEMH,FPAR_DT,LAI_DT,LANDCOVER_DT,pafFparH,pafLaiH,pafLandcoverH,pafEvaporation,pafTransporation,pafEvapotransporation,pafPEvapotransporation,paff2H,pafParH);
					
							comGPP(FPAR_DT, j, Lati_xrc, pafFparH, FPARnodata, pafLandcoverH, pafGppH,pafLUEoutp, pafLUEH, paff2H, pafTE2H, pafParH);

							comNPP(LAI_DT, LAImaxy_DT,j, Lati_xrc, pafLaiH, pafLaimaxH, pafLandcoverH, pafRmH, pafTavH, pafGppH, pafNppH, LAInodata);
					

							for(int k = 0; k < Lati_xrc; k++)
							{
								pafGPPShortIntH[k] = (DT_16U)(pafGppH[k] * SCALE_FACT);
								pafNPPShortIntH[k] = (DT_16U)(pafNppH[k] * SCALE_FACT);
					
							}
							//Write one row of the processed GPP image to the dataset
							gppDataset->RasterIO(GF_Write,0,j,Lati_xrc,1,pafGPPShortIntH,Lati_xrc,1,GDT_Int16,1,0,0,0,0);
							//Write one row of the processed NPP image to the dataset
							nppDataset->RasterIO(GF_Write,0,j,Lati_xrc,1,pafNPPShortIntH,Lati_xrc,1,GDT_Int16,1,0,0,0,0);
							//Write one row of the processed LUE image to the dataset
							lueDataset->RasterIO(GF_Write,0,j,Lati_xrc,1,pafLUEoutp,Lati_xrc,1,GDT_Float32,1,0,0,0,0);
							cout << flush << '\r' <<"Calculating"+HV+": Completed" <<(j+1)*100/Lati_yrc<<'%';

						}
						GDALClose(gppDataset);
						GDALClose(nppDataset);
						GDALClose(lueDataset);
						delete[] pafRsParML;
						pafRsParML=NULL;//ZHL
					}

					//Close the file pointer

					GDALClose(FPAR_DT);
					GDALClose(LAI_DT);
					GDALClose(LAImaxy_DT);
					GDALClose(LANDCOVER_DT);
					GDALClose(GLDAS_Tav_DT);
					GDALClose(GLDAS_par_DT);
					GDALClose(DEMM_DT);
					GDALClose(DEML_DT);
					GDALClose(Lati_DT);
					GDALClose(GLDAS_swnet_DT);
					GDALClose(GLDAS_lwnet_DT);
					GDALClose(GLDAS_dt_DT);
					//Free memory, free the pointer
					delete[] pafSWnetL;
					pafSWnetL = NULL;
					delete[] pafSWnetL;//zhl
					pafSWnetL = NULL;
					delete[] pafLWnetL;//zhl
					pafLWnetL = NULL;
					delete[] pafDTL ;//zhl
					pafDTL = NULL;
					delete[] pafGldasParL;//zhl
					pafGldasParL = NULL;
					delete[] pafTavL;//zhl
					pafTavL = NULL;
					delete[] pafTavGldasL;//zhl
					pafTavGldasL = NULL;
					delete[] pafDemL;//zhl
					pafDemL = NULL;
					delete[] pafDemM;//zhl
					pafDemM = NULL;
					delete[] pafLongH;//zhl
					pafLongH = NULL;
					delete[] pafLatH;//zhl
					pafLatH = NULL;
					delete[] pafLUEH;//zhl
					pafLUEH = NULL;
					delete[] pafTavH;//zhl
					pafTavH = NULL;
					delete[] pafTavGldasH;//zhl
					pafTavGldasH = NULL;
					delete[] pafToptH;//zhl
					pafToptH = NULL;
					delete[] pafK1H;//zhl
					pafK1H = NULL;
					delete[] pafK2H;//zhl
					pafK2H = NULL;
					delete[] pafTE2H;//zhl
					pafTE2H = NULL;
					delete[] pafEvaporation;//zhl
					pafEvaporation = NULL;
					delete[] pafTransporation;//zhl
					pafTransporation = NULL;
					delete[] pafEvapotransporation;//zhl
					pafEvapotransporation = NULL;
					delete[] pafPEvapotransporation;//zhl
					pafPEvapotransporation = NULL;
					delete[] paff2H;//zhl
					paff2H = NULL;
					delete[] CI;//zhl
					CI = NULL;
					delete[] pafGldasParH;//zhl
					pafGldasParH = NULL;
					delete[] pafSWnetH;//zhl
					pafSWnetH = NULL;
					delete[] pafLWnetH;//zhl
					pafLWnetH = NULL;
					delete[] pafDTH;//zhl
					pafDTH = NULL;
					delete[] pafRsParH;//zhl
					pafRsParH = NULL;
					delete[] pafParH;//zhl
					pafParH = NULL;
					delete[] pafFparH;//zhl
					pafFparH = NULL;
					delete[] pafLaiH;//zhl
					pafLaiH = NULL;
					delete[] pafLandcoverH;//zhl
					pafLandcoverH = NULL;
					delete[] pafGppH;//zhl
					pafGppH = NULL;
					delete[] pafLUEoutp;//zhl
					pafLUEoutp = NULL;
					delete[] pafRmH;//zhl
					pafRmH = NULL;
					delete[] pafNppH;//zhl
					pafNppH = NULL;
					delete[] pafLaimaxH;//zhl
					pafLaimaxH = NULL;
					delete[] pafGPPShortIntH;//zhl
					pafGPPShortIntH = NULL;
					delete[] pafNPPShortIntH;//zhl
					pafNPPShortIntH = NULL;

		}
		}
		}
		}

	//system("pause");//The cmd window exits automatically
	return 0;
}