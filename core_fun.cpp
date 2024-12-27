#include "stdafx.h"
#include "IMGALG.h"
#include <gdal_priv.h>
#include <iostream>
#include <fstream>
#include <cpl_conv.h>
#include <string.h>
#include "stdlib.h"
#include"ogr_spatialref.h"
#include "npp_fun.h"
#include "math.h"

#define PI 3.1415926

using namespace std;

/**
* The core algorithm assembly for GPP and NPP calculations
**/

///added by MJ Wang 2019.4.20
void comCI(GDALDataset *Lati_DT, int j, int d,int Lati_xrc, int GLDAS_xrc, double GLDAS_GeoTransform[6], float* pafLatH, float* pafParH,float* CI)				   
{
		
	    
	    //Read the latitude, with the latitude being the first band
		GDALRasterBand* rasterBand = Lati_DT->GetRasterBand(1);
		GDALDataType dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1,pafLatH,Lati_xrc,1, dataType,0,0);
	     
        
        for(int i=0;i<Lati_xrc;i++)
		{
			int day_y;
			if (d<45)
				day_y = 8*d+4;
			else
				day_y = 8*d+3;
			pafLatH[i] = pafLatH[i]>89.5?89.5:pafLatH[i];
			pafLatH[i] = pafLatH[i]<-89.5?-89.5:pafLatH[i];
			pafLatH[i] =pafLatH[i]* PI /180;
			float omg;///Time angle, unit: rad, at this time, it is assumed that the altitude angle at sunrise and sunset is α=0, so sinα=0 (after the meridian circle of the observation point coincides with the sun, that is, the angle of rotation of the earth at noon, from 0° to 360° every day, and the hour angle at noon is 0°)
			float Q;
			float Qday;
			
			float delta =0.409*sin(2*PI*day_y/365+4.89);///The declination of the sun, unit: rad (the angle between the center line of the sun and the earth and the plane of the equator), -23°27'~+23°27' radians
			float dm = 1+0.0344*cos(2*PI*day_y/365);///Correction factor for distances between the sun and the earth
			float cosomg = tan(pafLatH[i])*tan(delta);
			if (cosomg<1 && cosomg>-1)/////Non-polar day and night
			{
				omg = acos(-1*cosomg);///Time angle, unit: rad, at this time, it is assumed that the altitude angle at sunrise and sunset is α=0, so sinα=0 (after the meridian circle of the observation point coincides with the sun, that is, the angle of rotation of the earth at noon, from 0° to 360° every day, and the hour angle at noon is 0°)
				Q = omg*sin(pafLatH[i])*sin(delta)+cos(pafLatH[i])*cos(delta)*sin(omg);///The total amount of solar radiant energy received per unit area in the upper bound of the atmosphere
				Qday = 24*3600*1367*dm*dm*Q/(PI*1000000);///S0=1367, unit: W/m2//unit: MJ/(m2*day)	
				CI[i]=(pafParH[i])/Qday; //Replace it with ERA5 SSRD and remove *2
				CI[i]=CI[i]>1?1:CI[i];

			}
			else if (cosomg>=1)/////Polar Day
			{
				omg = PI;///Time angle, unit: rad, at this time, it is assumed that the altitude angle at sunrise and sunset is α=0, so sinα=0 (after the meridian circle of the observation point coincides with the sun, that is, the angle of rotation of the earth at noon, from 0° to 360° every day, and the hour angle at noon is 0°)
				Q = omg*sin(pafLatH[i])*sin(delta)+cos(pafLatH[i])*cos(delta)*sin(omg);///The total amount of solar radiant energy received per unit area in the upper bound of the atmosphere
				Qday = 24*3600*1367*dm*dm*Q/(PI*1000000);///S0=1367, unit: W/m2//unit: MJ/(m2*day)	
				CI[i]=(pafParH[i])/Qday; //Replace it with ERA5 SSRD and remove *2
				CI[i]=CI[i]>1?1:CI[i];

			}
			else ///Polar Night
			{
				omg = 0;///Time angle, unit: rad, at this time, it is assumed that the altitude angle at sunrise and sunset is α=0, so sinα=0 (after the meridian circle of the observation point coincides with the sun, that is, the angle of rotation of the earth at noon, from 0° to 360° every day, and the hour angle at noon is 0°)
				Q = omg*sin(pafLatH[i])*sin(delta)+cos(pafLatH[i])*cos(delta)*sin(omg);///The total amount of solar radiant energy received per unit area in the upper bound of the atmosphere
				Qday = 24*3600*1367*dm*dm*Q/(PI*1000000);///S0=1367, unit: W/m2//unit: MJ/(m2*day)	
				CI[i]=(pafParH[i])/Qday; //Replace it with ERA5 SSRD and remove *2
				CI[i]=CI[i]>1?1:CI[i];

			}
			
			
		}
}



//LUE is calculated based on land use type and CI////revised by MJ Wang 20190531
void computerLUE(GDALDataset *LANDCOVER_DT, int j, int Lati_xrc, DT_8U* pafLandcoverH, float* pafLUE, float* CI, GDALDataset *FPAR_DT, GDALDataset *LAI_DT, DT_8U*pafFparH, DT_16U*pafLaiH)
{
//Add FPAR and LAI to determine if FPAR=0 or LAI=0, LUE=0， added by zhl 20200910
	//FPAR500M
	    GDALRasterBand*	rasterBand = FPAR_DT->GetRasterBand(1);
	    GDALDataType dataType = rasterBand->GetRasterDataType();
	    rasterBand->RasterIO(GF_Read, 0, j, Lati_xrc, 1, pafFparH, Lati_xrc, 1, dataType, 0, 0);
	////GLASS的LAI500m
	    rasterBand = LAI_DT->GetRasterBand(1);
	    dataType = rasterBand->GetRasterDataType();
	    rasterBand->RasterIO(GF_Read, 0, j, Lati_xrc, 1, pafLaiH, Lati_xrc, 1, dataType, 0, 0);
		rasterBand = LANDCOVER_DT->GetRasterBand(1);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1,pafLandcoverH,Lati_xrc,1,dataType ,0,0);
		for(int i=0;i<Lati_xrc;i++)
		{

			if (CI[i]>=0 && CI[i]<=1)////Polar day or non-polar night polar day
			{

			switch(pafLandcoverH[i])
			{

				//revised by liqi 20190921
				case 12:pafLUE[i]=0.94278*CI[i]+2.46575*(1-CI[i]);break;		//12croplands	//wang20181008				
				case 4: pafLUE[i] =0.5757344*CI[i] + 2.5648967*(1 - CI[i]); break;		//4 DBF revised by liqi 20190921
				case 2: pafLUE[i] =0.5974848*CI[i] + 2.6058105*(1 - CI[i]); break;		//2 EBF  revised by liqi 20190921
				case 3: pafLUE[i] = 0.5737743*CI[i] + 2.5563368*(1 - CI[i]); break; //3DNF zhl20201002
				case 1: pafLUE[i] =0.5737743*CI[i] + 2.5563368*(1 - CI[i]); break;		//1ENF Evergreen Needleleaf Forests  revised by liqi 20190921
				case 5: pafLUE[i] =0.6725481*CI[i] + 2.469193*(1 - CI[i]); break;		//5Mixed Forests 	revised by liqi 20190921
				case 10:pafLUE[i] =0.5108*CI[i] + 2.39796*(1 - CI[i]); break;		//10Grasslands	
				case 8: pafLUE[i] =0.4757012*CI[i] + 2.3785058*(1 - CI[i]); break;		//8Woody Savannas 	 revised by liqi 20190921
				case 9: pafLUE[i] =0.5204*CI[i] + 2.4662*(1 - CI[i]); break;		//9Savannas	
				case 6: pafLUE[i] =0.4673*CI[i] + 2.0701*(1 - CI[i]); break;		//6Closed Shrublands 	
				case 7: pafLUE[i] =0.4673*CI[i] + 2.0701*(1 - CI[i]); break; //7 OSH  zhl20201002
				case 11:pafLUE[i] = 0.5108*CI[i] + 2.39796*(1 - CI[i]); break;		//11Permanent Wetlands 		
				case 17: pafLUE[i]=0;break;			//0Water Bodies 						
				case 13:pafLUE[i] =0.5108*CI[i] + 2.39796*(1 - CI[i]); break;		//13Urban and Built-up Lands 		
				case 16:pafLUE[i] =0.5108*CI[i] + 2.39796*(1 - CI[i]); break;			//16Barren		
				case 15:pafLUE[i]=0;break;			//15Permanent Snow and Ice 					
				case 255:pafLUE[i]=0;break;			//255Unclassified		 
				case 14:pafLUE[i] =0.94278*CI[i] + 2.46575*(1 - CI[i]); break;		//14Cropland/Natural Vegetation Mosaics 	
				default:pafLUE[i]=0;break;			//		Default value
			}
			}
			else
			{
				
				pafLUE[i] = 0;////In the case of polar night, LUE=0

			    
			}
			if (pafFparH[i]==0 )
			{
				pafLUE[i] = 0;
			}
			if (pafFparH[i] == 255)
			{
				pafLUE[i] = 0;
			}
			if (pafLaiH[i] > 1000 || pafLaiH[i] ==0) 
			{
				pafLUE[i] = 0;
			}//added by zhl 20200910
			
		}
}




//Calculate the temperature stress factor
//All calculations are based on latitude and longitude data
//revised by MJ Wang
void comTempStress(GDALDataset *Lati_DT,GDALDataset *DEMM_DT, int j, int Lati_xrc, int GLDAS_xrc,DT_8U* pafLandcoverH,int DEMM_xrc, double GLDAS_GeoTransform[6],double STRMM_DEM_GeoTransform[6],
				   float* pafLongH, float* pafLatH, float* pafTavH, float* pafTavL,float* pafToptH,
				   float* pafDemM, float* pafDEMH, double* pafK1H, double* pafK2H, double* pafTE2H)
{
		//Read the longitude, the longitude is the second band
		GDALRasterBand* rasterBand = Lati_DT->GetRasterBand(2);
		GDALDataType dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafLongH,Lati_xrc,1, dataType,0,0);
		//Read the latitude, with the latitude being the first band
		rasterBand = Lati_DT->GetRasterBand(1);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1,pafLatH,Lati_xrc,1, dataType,0,0);
		//Read high-resolution DEMM data
	    rasterBand=DEMM_DT->GetRasterBand(1);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafDEMH,Lati_xrc,1, dataType,0,0);

		//============According to the latitude and longitude, the coordinates of 4 points were found from the image of 0.25° (5km, 500m).=====================
		//=============Because the coordinates of the center point of the cell are required, they are processed on latitude and longitude=====================

		
		//The latitude and longitude coordinates and pixel size of the pixel center of the upper left point of ERA5 data (0.25 degrees).//revised by MJ Wang 2019.4.19
		float GLDAS_pix_size = GLDAS_GeoTransform[1];
		float GLDAS_Lup_long = GLDAS_GeoTransform[0] + GLDAS_pix_size/2;
		float GLDAS_Lup_lat = GLDAS_GeoTransform[3] - GLDAS_pix_size/2;
		

		//The latitude and longitude coordinates and pixel size of the pixel center of the upper left point of the high-resolution DEM (500m).//revised by MJ Wang 2019.4.19
		float STRMM_DEM_pix_size = STRMM_DEM_GeoTransform[1];
		float STRMM_DEM_Lup_long = STRMM_DEM_GeoTransform[0] + STRMM_DEM_pix_size/2;
		float STRMM_DEM_Lup_lat = STRMM_DEM_GeoTransform[3] - STRMM_DEM_pix_size/2;


		//Resampling Air Average Temperature - Bilinear (0.25°C to 500m)//revised by MJ Wang 2019.4.19
		resampleL2H(pafLongH,pafLatH,GLDAS_pix_size,GLDAS_Lup_long,GLDAS_Lup_lat,pafTavH,1,Lati_xrc,pafTavL,GLDAS_xrc);
		

		//Temperature correction, using high-resolution image correction to correct the resampled ERA5 temperature data and convert the unit to K
		//T=Tav-(DEMH-DEML)*0.0065
		//The corrected temperature pafTavH is generated, where the average air temperature is Celsius and the optimal temperature unit is K, so the average temperature does not need to be subtracted K_value
		//The resampled pafTavH (500m) is corrected using a high-resolution DEM (500m).
		for(int i=0;i<Lati_xrc;i++)
		{
			float dt=pafDEMH[i]*AIR_COR;
			pafTavH[i]=pafTavH[i]-dt;
			//pafTavH[i]=pafTavH[i]-K_value;
			//pafToptH[i]=pafToptH[i]-dt-K_value;
			switch(pafLandcoverH[i])
			{
				case 12: pafToptH[i]=25;break;		//12croplands	wang20181008
				case 4: pafToptH[i]=20;break;		//4Deciduous Broadleaf Forests	
				case 2: pafToptH[i]=25;break;		//2Evergreen Broadleaf Forests 	
				case 3: pafToptH[i]=15;break;		//3Deciduous Needleleaf Forests 	
				case 1: pafToptH[i]=15;break;		//1Evergreen Needleleaf Forests 	
				case 5: pafToptH[i]=17;break;		//5Mixed Forests 	
				case 10: pafToptH[i]=18;break;		//10Grasslands	
				case 8: pafToptH[i]=19;break;		//8Woody Savannas 	
				case 9: pafToptH[i]=20;break;		//9Savannas	
				case 6: pafToptH[i]=20;break;		//6Closed Shrublands 	
				case 7: pafToptH[i]=16;break;		//7Open Shrublands 	//zhl20201002
				case 11: pafToptH[i]=18;break;		//11Permanent Wetlands 		
				case 17: pafToptH[i]=0;break;			//0Water Bodies 		
				case 13: pafToptH[i]=18;break;		//13Urban and Built-up Lands 		
				case 16: pafToptH[i]=18;break;			//16Barren		
				case 15: pafToptH[i]=0;break;			//15Permanent Snow and Ice 	
				case 255: pafToptH[i]=0;break;			//255Unclassified		 
                case 14: pafToptH[i]=25;break;		//14Cropland/Natural Vegetation Mosaics 	
				default:pafToptH[i]=0;break;			//		default value
			}



		}
		//Calculate the three parameters of K1, K2, and TE2
		
		for(int i=0;i<Lati_xrc;i++)
		{
			pafK1H[i]=1.0+exp(0.2*(pafToptH[i]-10-pafTavH[i]));
			pafK2H[i]=1.0+exp(0.3*(pafToptH[i]*(-1)-10+pafTavH[i]));
			pafTE2H[i]=1.1814/(pafK1H[i]*pafK2H[i]);
			if (pafTE2H[i]<0)
			{
				pafTE2H[i]=0;
			}
			
			
			//added by Sun Rui 2019.09.23
			if(pafLandcoverH[i] != 1 )
			{
				if (pafTavH[i]<(-5)) pafTE2H[i] =0;
			}
			if(pafLandcoverH[i] == 1 )
			{
				if (pafTavH[i]<(-10)) pafTE2H[i] =0;
			}
			

		}
		
}


//The water stress factor was calculated, and the actual evapotranspiration and potential evapotranspiration were calculated by using ERA5 meteorological data, LAI and FPAR, and then the water stress factor za was obtained

void comET(GDALDataset *Lati_DT, int j, int Lati_xrc, int GLDAS_xrc,  double GLDAS_GeoTransform[6],
				   float* pafLongH, float* pafLatH, float* pafTavH, float* pafTavGldasL,float* pafSWnetH, float* pafSWnetL, float* pafLWnetH, float* pafLWnetL, float* pafDTH, float* pafDTL,float* pafDEMH,
				   GDALDataset *FPAR_DT, GDALDataset *LAI_DT, GDALDataset *LANDCOVER_DT, DT_8U* pafFparH, DT_16U* pafLaiH, DT_8U* pafLandcoverH, double* pafEvaporation, double* pafTransporation, double* pafEvapotransporation, double* pafPEvapotransporation, float* paff2H, float* pafParH)/////revised by MJ Wang 2019.4.24

				   
{
	    ////FPAR DATA of GLASS
        GDALRasterBand*	rasterBand=FPAR_DT->GetRasterBand(1);
		GDALDataType dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafFparH,Lati_xrc,1, dataType,0,0);
		////GLASS LAI500m
        rasterBand=LAI_DT->GetRasterBand(1);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafLaiH,Lati_xrc,1, dataType,0,0);
		///MODIS Landcover
        rasterBand=LANDCOVER_DT->GetRasterBand(1);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafLandcoverH,Lati_xrc,1, dataType,0,0);

        //Read the longitude, the longitude is the second band
		rasterBand = Lati_DT->GetRasterBand(2);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafLongH,Lati_xrc,1, dataType,0,0);
		//Read the latitude, with the latitude being the first band
		rasterBand = Lati_DT->GetRasterBand(1);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1,pafLatH,Lati_xrc,1, dataType,0,0);
		//============According to the latitude and longitude, the coordinates of 4 points were found from the image of 0.25° (5km, 500m).=====================
		//=============Because the coordinates of the center point of the cell are required, they are processed on latitude and longitude=====================

		//The latitude and longitude coordinates and pixel size of the pixel center of the upper left point of ERA5 data (0.25 degrees).////added by MJ Wang 2019.4.20
		float GLDAS_pix_size = GLDAS_GeoTransform[1];
		float GLDAS_Lup_long = GLDAS_GeoTransform[0] + GLDAS_pix_size/2;
		float GLDAS_Lup_lat = GLDAS_GeoTransform[3] - GLDAS_pix_size/2;

		//added by Sun Rui  20190103
		float clumping_index;

        //Resampling Net Shortwave Radiation - Bilinear (0.25 degrees to 500m) // added by MJ Wang 2019.3.28
		resampleL2H(pafLongH,pafLatH,GLDAS_pix_size,GLDAS_Lup_long,GLDAS_Lup_lat,pafSWnetH,1,Lati_xrc,pafSWnetL,GLDAS_xrc);
        //Resampling Net Longwave Radiation - Bilinear (0.25 degrees to 500m)  //added by MJ Wang 2019.3.28
		resampleL2H(pafLongH,pafLatH,GLDAS_pix_size,GLDAS_Lup_long,GLDAS_Lup_lat,pafLWnetH,1,Lati_xrc,pafLWnetL,GLDAS_xrc);
		//Resampling Dewpoint Temperature - Bilinear (0.25 degrees to 500m)    added by MJ Wang 20190812
		resampleL2H(pafLongH,pafLatH,GLDAS_pix_size,GLDAS_Lup_long,GLDAS_Lup_lat,pafDTH,1,Lati_xrc,pafDTL,GLDAS_xrc);
		//Resampling uncorrected temperature data // Resampling uncorrected temperature (0.25 degrees to 500m) // added by MJ Wang 2019.3.28
		resampleL2H(pafLongH,pafLatH,GLDAS_pix_size,GLDAS_Lup_long,GLDAS_Lup_lat,pafTavH,1,Lati_xrc,pafTavGldasL,GLDAS_xrc);//2015.9.15增加

        for(int i=0;i<Lati_xrc;i++)
		{

		    float TavK=pafTavH[i]+K_value;
		    float TavC=pafTavH[i];//2015.9.16
			float a=(pafSWnetH[i]+pafLWnetH[i]);///revised by sunrui 20191022,

			//added by Sun Rui
			// calculate fraction of vegetation cover
			 switch(pafLandcoverH[i])
			{
				case 12: clumping_index=0.9;break;		//cropland
				case 4: clumping_index=0.8;break;		//DBF
				case 2: clumping_index=0.8;break;		//EBF
				case 3: clumping_index=0.6;break;		//DNF
				case 1: clumping_index=0.6;break;		//ENF
				case 5: clumping_index=0.7;break;		//mixed forest
				case 10: clumping_index=0.9;break;		//grasslands
				case 8: clumping_index=0.8;break;		//wooded savanas
				case 9: clumping_index=0.8;break;		//savanas
				case 6: clumping_index=0.8;break;		//Closed Shrublands
				case 7: clumping_index=0.8;break;		//Open Shrublands
				case 11: clumping_index=0.9;break;		//Permanent Wetlands
				case 17: clumping_index=0.9;break;			//water bodies
				case 13: clumping_index=0.9;break;		//urban and built-up
				case 16: clumping_index=0.9;break;			//Barren
				case 15: clumping_index=0.9;break;			//Permanent Snow and Ice
				case 255: clumping_index=0.9;break;			//background
				default:clumping_index=0.9;break;			//default value
				
			}
			///calculate energy components for the canopy[acanopy] and soil surface [asoil].  //added by MJ Wang 20190510
			float pafFcoverH=1-exp(-0.5*clumping_index*pafLaiH[i]*0.01);//LAI SCALE FACTOR 1000,FPAR SCALE FACTOR 1 zhl20200925 MUSES LAI为100
			float acanopy=a*pafFcoverH;//unit:W/m2 //added by MJ Wang 20190510
			float asoil=a*(1-pafFcoverH);//unit:W/m2 //added by MJ Wang 20190510

			///calculate saturation water vapor pressure at air temperature unit: Pa  //added by MJ Wang 20190510
		    float esat_temp=7.5*TavC/(237.3+TavC);
		    double esat=611*pow(10,esat_temp);////saturation water vapor pressure at air temperature unit: Pa  //added by MJ Wang 20190510
			
			///calculate the slope of the curve relating saturated water vapor pressure (Pa) to air temperature
		    double slope=esat*(6463/(273+TavC)-3.927)/(273+TavC);  //unit: Pa/K  //added by MJ Wang 20190510
			/////revised by MJ Wang 20190812 this method uses dewpoint temperature to caculate actual water vapor pressure (Pa)
			float eact_temp=7.5*pafDTH[i]/(237.3+pafDTH[i]);
		    double eact=611*pow(10,eact_temp);////actual water vapor pressure unit: Pa
			float p_temp=(293-0.0065*pafDEMH[i])/293;//unit: K///
		    double p=101300*pow(p_temp,5.26);///unit: Pa


			///calculate air density[rho] unit: Kg/m^3  //added by MJ Wang 20190510
			float rho=p/(287.05*TavK*(1+0.378*eact/p));///wang20181008//unit: Kg/m^3 
			///calculate vapor pressure deficit[vpd] unit: Pa  //added by MJ Wang 20190510
		    float vpd=esat-eact;
			if (vpd < 0)
			{
				vpd=0;
			}

		    float lammda=2.5-0.002361*TavC;///the latent heat of vaporization 
		    float gamma=1.013*p/(0.622*lammda)/1000;///psychrometric constant  
			///The values of KQ are from [Leuning, R., Zhang, Y. Q., Rajaud, A., Cleugh, H., & Tu, K. (2008). A simple surface conductance model to estimate regional evaporation using MODIS leaf area index and the Penman‐Monteith equation. Water Resources Research, 44(10).]
            float KQ=0.6;//revised by MJ Wang 2019.5.8
			//The values of Q50,D50 are from [Zhang, Y. Q., Chiew, F. H. S., Zhang, L., Leuning, R., & Cleugh, H. A. (2008). Estimating catchment evaporation and runoff using MODIS leaf area index and the Penman‐Monteith equation. Water Resources Research, 44(10).]
			int Q50=2.6;
			int D50=800;////revised by MJ Wang 2019.4.24
			
            float gsx;
            float ga;
			float gas;
            float gtot;
			
			float gch;
            int k;

			////revised by MJ Wang 2019.5.10
			////the values of gsx, ga (canopy), gtot,k, gch are from [Zhang, K., Kimball, J. S., Nemani, R. R., & Running, S. W. (2010). A continuous satellite‐derived global record of land surface evapotranspiration from 1983 to 2006. Water Resources Research, 46(9).]
			switch(pafLandcoverH[i])
			{


				//revised by Sun Rui 20190925
				case 12: gsx=0.0069; ga=0.001;gtot=0.003;k=500;gch=0.03;break;		//	cropland
				case 4: gsx=0.0047; ga=0.02;gtot=0.002;k=200;gch=0.01;break;		 //	DBF
				case 2: gsx=0.0047; ga=0.02;gtot=0.002;k=200;gch=0.01;break;		//	EBF
				case 3: gsx=0.0047; ga=0.03;gtot=0.002;k=150;gch=0.08;break;		//	DNF
				case 1: gsx=0.0047; ga=0.03;gtot=0.002;k=150;gch=0.08;break;		//	ENF
				case 5: gsx=0.0047; ga=0.025;gtot=0.002;k=200;gch=0.045;break;		//	mixed forest
				case 10: gsx=0.0040; ga=0.002;gtot=0.002;k=450;gch=0.04;break;		//	grasslands
				case 8: gsx=0.005; ga=0.002;gtot=0.0018;k=300;gch=0.04;break;		//	wooded savanas
				case 9: gsx=0.004; ga=0.002;gtot=0.001;k=400;gch=0.04;break;		//	savanas
				case 6: gsx=0.0045; ga=0.01;gtot=0.003;k=400;gch=0.04;break;		//	Closed Shrublands
				case 7: gsx=0.0045; ga=0.01;gtot=0.004;k=200;gch=0.04;break;		//  Open Shrublands
				case 11: gsx=0.0050; ga=0.001;gtot=0.003;k=450;gch=0.04;break;		//	Permanent Wetlands
				case 17: gsx=0.0050; ga=0.001;gtot=0.003;k=450;gch=0.04;break;		//	water bodies
				case 13: gsx=0.0050; ga=0.001;gtot=0.003;k=450;gch=0.04;break;		//	urban and built-up
				case 14: gsx=0.0069; ga=0.001;gtot=0.003;k=500;gch=0.03;break;		//	crop mosaic
				case 16: gsx=0.0050; ga=0.001;gtot=0.003;k=450;gch=0.04;break;		//	Barren
				case 15: gsx=0.0050; ga=0.001;gtot=0.003;k=450;gch=0.04;break;		//	Permanent Snow and Ice
				case 255: gsx=0.0050; ga=0.001;gtot=0.003;k=450;gch=0.04;break;		//	background
				default:gsx=0.0050; ga=0.001;gtot=0.003;k=450;gch=0.04;break;		//	default value


			}

			//Half of shortwave radiation///////////
			////calculate the canopy conductance[gc]. the following equations are from [Zhang, Y. Q., Chiew, F. H. S., Zhang, L., Leuning, R., & Cleugh, H. A. (2008). Estimating catchment evaporation and runoff using MODIS leaf area index and the Penman‐Monteith equation. Water Resources Research, 44(10).]
            float qh=pafParH[i]*0.48; //QH=PAR,The input data is the total radiation, so *0.48
            double gc_temp1 = (-1)*KQ*pafLaiH[i]*0.01;//LAI SCALE FACTOR 1000,FPAR SCALE FACTOR 1 zhl20200925 MUSES LAI为100
            double gc_temp2=(qh+Q50)/(qh*exp(gc_temp1)+Q50);
            double gc=gsx/KQ*log(gc_temp2)*(1/(1+vpd/D50));
			///calculate Transpitation of canopy, unit: mm
			pafTransporation[i]=(((slope*acanopy+rho*1004*vpd*ga)/(slope+gamma*(1+ga/gc)))/lammda)*24*60*60/1000000;//wang20181008
			pafTransporation[i] = pafTransporation[i]<0?0:pafTransporation[i];////added by MJ Wang 2019.4.24

            float RH=eact/esat;
			if (RH > 1)
			{
				RH=1;
			}

			///calculate Evaporation of soil, the equations are from [Zhang, K., Kimball, J. S., Nemani, R. R., & Running, S. W. (2010). A continuous satellite‐derived global record of land surface evapotranspiration from 1983 to 2006. Water Resources Research, 46(9).]
            float Gcorr=((273.15+TavK)/293.15)*(101300/p);
            float gtotc=gtot*Gcorr;//unit: m/s
            double t_temp=pow(RH,(vpd/k));/////moisture constraint on soil evaporation [t_temp]
			double sigma=5.67e-8;//unit: W/(m2·K4)
			double grh = (4.0*sigma*pow(TavK,3))/(rho*1004);//unit: m/s
			gas = grh+gch;          
			pafEvaporation[i]=((t_temp*((slope*asoil+rho*1013*vpd*gas))/(slope+gamma*gas/gtotc))/lammda)*24*60*60/1000000;
			pafEvaporation[i]=pafEvaporation[i]<0?0:pafEvaporation[i];////added by MJ Wang 2019.4.24
			

            pafEvapotransporation[i]=pafEvaporation[i]+pafTransporation[i];
            pafPEvapotransporation[i]=((1.26*a*slope/(slope+gamma))/lammda)*24*60*60/1000000;
			pafPEvapotransporation[i] = pafPEvapotransporation[i]<0?0:pafPEvapotransporation[i];////added by MJ Wang 2019.4.24

		    if(pafLandcoverH[i]==17||pafLandcoverH[i]==15||pafLandcoverH[i]==255)///revised by MJ Wang 0water bodies,15Permanent Snow and Ice,255Unclassified default ET/PET=-1000 
            {
                pafEvapotransporation[i]=-1000;
                pafPEvapotransporation[i]=-1000;
            }

			double paff2H_temp;

			if(pafLandcoverH[i]==11) paff2H_temp=1.0;//Wetland water limiting factors =1
			else if (pafEvapotransporation[i]==0)
				paff2H_temp = 0.5;
			else if((pafPEvapotransporation[i]==0) &&(pafLandcoverH[i]==1||pafLandcoverH[i]==3||pafLandcoverH[i]==5))
				paff2H_temp = 1.0;
			else if ((pafPEvapotransporation[i]==0) &&(pafLandcoverH[i]==17||pafLandcoverH[i]==2||pafLandcoverH[i]==4||pafLandcoverH[i]>=6))
				paff2H_temp = 0.5;
			else
				paff2H_temp = 0.5+0.5*(pafEvapotransporation[i]/pafPEvapotransporation[i]);

            paff2H[i]=(float) paff2H_temp;
			
			if(paff2H[i]<0.5)
            {
                paff2H[i]=0.5;
            }
            if(paff2H[i]>1)
            {
                paff2H[i]=1;
            }

			

		}

}
//If there is no valid value for the par of remote sensing inversion, the par data calculated by ERA5 is used instead
void comPar(GDALDataset *Lati_DT, int j, bool IsDataExist, int Lati_xrc, int GLDAS_xrc, int RsPar_xrc, double GLDAS_GeoTransform[6], double RsPar_GeoTransform[6],
	float* pafLongH, float* pafLatH, float* pafsGldasParL, float* pafGldasParH, float* pafRsParML, float* pafRsParH, float* pafParH, float nodata)
{
	//6 lines below added by Sun Rui 20210226
	//Read the longitude, the longitude is the second band
	GDALRasterBand* rasterBand = Lati_DT->GetRasterBand(2);
	GDALDataType dataType = rasterBand->GetRasterDataType();
	rasterBand->RasterIO(GF_Read, 0, j, Lati_xrc, 1, pafLongH, Lati_xrc, 1, dataType, 0, 0);
	//Read the latitude, with the latitude being the first band
	rasterBand = Lati_DT->GetRasterBand(1);
	dataType = rasterBand->GetRasterDataType();
	rasterBand->RasterIO(GF_Read, 0, j, Lati_xrc, 1, pafLatH, Lati_xrc, 1, dataType, 0, 0);
	//============According to the latitude and longitude, the coordinates of 4 points were found from the image of 0.25° (5km, 500m).=====================
	//=============Because the coordinates of the center point of the cell are required, they are processed on latitude and longitude=====================

	//The latitude and longitude coordinates and pixel size of the center of the pixel at the upper left point of the ERA5 image (0.25 degrees).
	float GLDAS_pix_size = GLDAS_GeoTransform[1];
	float GLDAS_Lup_long = GLDAS_GeoTransform[0] + GLDAS_pix_size / 2;
	float GLDAS_Lup_lat = GLDAS_GeoTransform[3] - GLDAS_pix_size / 2;

	//Resampling the PAR generated by ERA5 with a bilinear method (0.25 degrees to 500 m)
	resampleL2H(pafLongH, pafLatH, GLDAS_pix_size, GLDAS_Lup_long, GLDAS_Lup_lat, pafGldasParH, 1, Lati_xrc, pafsGldasParL, GLDAS_xrc);

	if (IsDataExist)
	{
		//The latitude and longitude coordinates and pixel size of the center of the pixel at the upper left point of the PAR image (5km) inverted by remote sensing
		float RsPar_pix_size = RsPar_GeoTransform[1];
		float RsPar_Lup_long = RsPar_GeoTransform[0] + RsPar_pix_size / 2;
		float RsPar_Lup_lat = RsPar_GeoTransform[3] - RsPar_pix_size / 2;

		//PAR resampled by remote sensing inversion by bilinear method (5km degree to 500m)
		resampleL2H(pafLongH, pafLatH, RsPar_pix_size, RsPar_Lup_long, RsPar_Lup_lat, pafRsParH, 1, Lati_xrc, pafRsParML, RsPar_xrc, nodata);

		//If there is no valid value for the par of the remote sensing inversion, the par data calculated by ERA5 is used instead
		for (int i = 0; i<Lati_xrc; i++)
		{
			//Note: Assuming that the range of valid values is less than 200, this value will need to be modified after obtaining the PAR of the final remote sensing inversion
			if (pafRsParH[i] < 200 && pafRsParH[i] >0)
			{
				pafParH[i] = pafRsParH[i];
			}
			else
			{
				pafParH[i] = pafGldasParH[i];
			}
		}
	}
	else
	{
		for (int i = 0; i<Lati_xrc; i++)
		{
			pafParH[i] = pafGldasParH[i];
		}
	}
}

//Calculate GPP
void comGPP(GDALDataset *FPAR_DT, int j, int Lati_xrc, DT_8U* pafFparH, float FPARnodata, DT_8U* pafLandcoverH, double* pafGppH, float* pafLUEoutp, float* pafLUEH, float* paff2H, double* pafTE2H, float* pafParH)///revised by MJ Wang 20190601
////revised by MJ Wang 2019.4.14
{
		GDALRasterBand*	rasterBand=FPAR_DT->GetRasterBand(1);
		GDALDataType dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafFparH,Lati_xrc,1, dataType,0,0);

		
		for(int i=0;i<Lati_xrc;i++)
		{
			switch(pafLandcoverH[i])
				{
					case 17:
						pafGppH[i]=0;//water bodies
						pafLUEoutp[i] = 0;break;
					case 15:
						pafGppH[i]=0;//Permanent Snow and Ice
						pafLUEoutp[i] = 0;break;
					case 255:
						pafGppH[i]=0;//background
						pafLUEoutp[i] = 0;break;
					default:
						pafFparH[i]==FPARnodata?pafGppH[i]=0:pafGppH[i]=pafLUEH[i]*paff2H[i]*pafTE2H[i]*(pafParH[i]*0.48)*(pafFparH[i]*0.004);
						pafLUEoutp[i]=pafLUEH[i]*paff2H[i]*pafTE2H[i];///added by MJ Wang 20190601 in order to output LUE
					    break;
			}

		}
}

//Calculate NPP
// revised by MJ Wang
void comNPP(GDALDataset *LAI_DT, GDALDataset *LAImaxy_DT,int j, int Lati_xrc, 
	DT_16U* pafLaiH, double* pafLaimaxH, DT_8U* pafLandcoverH, float* pafRmH, float* pafTavH,
			double* pafGppH, double* pafNppH,float LAInodata)
{

		//============According to the latitude and longitude, the coordinates of 4 points were found from the 500m image=====================
		//=============Because the coordinates of the center point of the cell are required, they are processed on latitude and longitude=====================


	    //Read LAImax data
		GDALRasterBand* rasterBand=LAImaxy_DT->GetRasterBand(1);
		GDALDataType dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafLaimaxH, Lati_xrc,1, dataType,0,0);

	    //Read LAI data
		rasterBand=LAI_DT->GetRasterBand(1);
		dataType =rasterBand->GetRasterDataType();
		rasterBand->RasterIO( GF_Read,0,j,Lati_xrc,1, pafLaiH, Lati_xrc,1, dataType,0,0);

		for(int i=0;i<Lati_xrc;i++)
		{
			float SLA,stema,stemb,rm1,rm2,rm31,rm32,frleaf,crstem,Q10;
			///revised by MJ Wang
            if(pafLaiH[i] > LAInodata)   //GLASS 500M LAI valid range=[0,1000],LAInodata is 1000，revised by zhl
			{
				pafRmH[i] = 0;
				pafNppH[i] = 0;
			}
			else if (pafLandcoverH[i]==17||pafLandcoverH[i]==15||pafLandcoverH[i]==255)
			{
				pafRmH[i] = 0;
				pafNppH[i] = 0;
			}
			else
			{
				//revised by MJ Wang 2019.4.2 
				//10grasslands,11permanent wetlands,13urban and build-up lands
				if(pafLandcoverH[i]==10||pafLandcoverH[i]==11||pafLandcoverH[i]==13||pafLandcoverH[i]==16)
				{
					float SLA_grass = 49;
					Q10 =2.0;
					float coef_leaf_respir = 0.006976;     //added by sunrui 20190926
					float coef_root_respir = 0.004152;     //added by sunrui 20190926
					//Calculation of leaf biomass gC/m2
					float biomass_leaf = 1000*pafLaiH[i]*0.01/SLA_grass;//revised by zhl，MuSyQLAI SCALE-FACTOR=1000,MUSES LAI为100
					//Calculate fine root biomass gC/m2
					float biomass_root = 2*1000*(pafLaimaxH[i]*0.01/SLA_grass);//revised by zhl，MuSyQLAI SCALE-FACTOR=1000，MUSES LAI为100
					float exponent = (pafTavH[i]-20)/10;

					pafRmH[i]=(biomass_root*coef_root_respir*pow(Q10,exponent))+(biomass_leaf*coef_leaf_respir*pow(Q10,exponent));
				}
				///8woody savannas = 45% DBF + 55% grass are weighted according to the proportion of components to maintenance respiration// revised by MJ Wang 20190516
				if(pafLandcoverH[i]==8)
				{
					////Calculate the maintenance respiration of the grass fraction
					float SLA_grass = 49;
					Q10 =2.0;
					float coef_leaf_respir = 0.00872;                 
					float coef_root_respir = 0.00519;                
					//Calculation of leaf biomass gC/m2
					float biomass_leaf = 1000*pafLaiH[i]*0.01/SLA_grass;//revised by zhl，MuSyQLAI SCALE-FACTOR=1000，MUSES LAI为100
					//Calculate fine root biomass gC/m2
					float biomass_root = 2*1000*(pafLaimaxH[i]*0.01/SLA_grass);//revised by zhl，MuSyQLAI SCALE-FACTOR=1000，MUSES LAI为100
					float exponent = (pafTavH[i]-20)/10;
					float pafRmH_grass = (biomass_root*coef_root_respir*pow(Q10,exponent))+(biomass_leaf*coef_leaf_respir*pow(Q10,exponent));

					////Calculate the maintenance respiration of the DBF fraction
					float SLA_DBF = 32.0;                                 //revised by MJ Wang 2019.4.16
					stema = -0.5189;                         // The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 15.649;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					
					rm1=0.0119;                              //Maintenance respiration coefficient of the leaves 
					rm2=0.00512;                               //Maintenance respiration coefficient of stems 
					rm31=0.00448;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.004992;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.2;                                  //ratio of fine root and leaf 
					crstem= 0.22;                                  //ratio of coarse root and stem
					                                      
					
					float pafRmH_DBF = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA_DBF,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
					
					pafRmH[i] = 0.55*pafRmH_grass + 0.45*pafRmH_DBF;
				}

				///9 Savanna = 20% DBF + 80% grass, eighted according to the proportion of components to maintenance respiration// revised by MJ Wang 20190516
				if(pafLandcoverH[i]==9)
				{
					////Calculate the maintenance respiration of the grass fraction
					float SLA_grass = 49;
					Q10 =2.0;
					float coef_leaf_respir = 0.00872;          
					float coef_root_respir = 0.00519;          
					//Calculation of leaf biomass gC/m2
					float biomass_leaf = 1000*pafLaiH[i]*0.01/SLA_grass;//revised by zhl，MuSyQLAI SCALE-FACTOR=1000MUSES LAI为100
					//Calculate fine root biomass gC/m2
					float biomass_root = 2*1000*(pafLaimaxH[i]*0.01/SLA_grass);//revised by zhl，MuSyQLAI SCALE-FACTOR=1000MUSES LAI为100
					float exponent = (pafTavH[i]-20)/10;
					float pafRmH_grass = (biomass_root*coef_root_respir*pow(Q10,exponent))+(biomass_leaf*coef_leaf_respir*pow(Q10,exponent));

					////Calculate the maintenance respiration of the DBF fraction
					float SLA_DBF = 32.0;                                 //revised by MJ Wang 2019.4.16
					stema = -0.5189;                         // The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 15.649;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					
					rm1=0.0119;                              //aintenance respiration coefficient of the leaves 
					rm2=0.00512;                               //Maintenance respiration coefficient of stems 
					rm31=0.00448;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.004992;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.2;                                  //ratio of fine root and leaf 
					crstem= 0.22;                                  //ratio of coarse root and stem
					                                     
					
					float pafRmH_DBF = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA_DBF,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
					
					pafRmH[i] = 0.8*pafRmH_grass + 0.2*pafRmH_DBF;
				}

				////12Croplands,14Cropland/Natural Vegetation Mosaics 
				if(pafLandcoverH[i]==12||pafLandcoverH[i]==14)
				{
					float SLA_grass = 49;
					Q10 =2.0;
					float coef_leaf_respir = 0.006976;     //added by sunrui 20190925
					float coef_root_respir = 0.004152;     //added by sunrui 20190925
					//Calculation of leaf biomass gC/m2
					float biomass_leaf = 1000*pafLaiH[i]*0.01/SLA_grass;//revised by zhl，MuSyQLAI SCALE-FACTOR=1000 
					//Calculate fine root biomass gC/m2
					float biomass_root = 2*1000*(pafLaiH[i]*0.01/SLA_grass);//revised by zhl，MuSyQLAI SCALE-FACTOR=1000 
					float exponent = (pafTavH[i]-20)/10;

					pafRmH[i]=(biomass_root*coef_root_respir*pow(Q10,exponent))+(biomass_leaf*coef_leaf_respir*pow(Q10,exponent));
				}

				//revised by MJ Wang
				//Deciduous Needleleaf Forests 
				if(pafLandcoverH[i]==3)
				{
					
					SLA=15.0;                                 //SLA is the specific leaf area, through which the conversion of biomass to LAI can be realized, reference: Earth Interactions x Volume 4 (2000) No.3  Page 1,Michael A.White
					stema = 0.1962;                         // The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 10.074;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					rm1=0.0087;                              //Maintenance respiration coefficient of the leaves
					rm2=0.004;                               //Maintenance respiration coefficient of stems 
					rm31=0.0032;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.0029;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.4;                                  //ratio of fine root and leaf 
					crstem= 0.29;                                  //ratio of coarse root and stem
					Q10=2.0;                                      //Q10 is set to 2.0
					
					pafRmH[i] = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
				}

				//Evergreen Needleleaf Forests
				if(pafLandcoverH[i]==1)
				{
					
					SLA=15.0;                                 //revised by MJ Wang 2019.4.16//2019.5.18
					stema = -2.761;                         // The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 50.564;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					rm1=0.0057;                              //Maintenance respiration coefficient of the leaves   
					rm2=0.004;                               //Maintenance respiration coefficient of stems 
					rm31=0.0036;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.0033;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.4;                                  //ratio of fine root and leaf 18
				    crstem= 0.29;                                  //ratio of coarse root and stem
					Q10=2.0;                                      //Q10 is set to 2.0
					
					pafRmH[i] = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
				}
				


				//revised by MJ Wang
				//Deciduous Broadleaf Forests 
				if(pafLandcoverH[i]==4)
				{
					SLA=32.0;                                 //revised by MJ Wang 2019.4.16
					stema = -1.0378;                         //The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 31.298;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					rm1=0.0119;                              //Maintenance respiration coefficient of the leaves
					rm2=0.00512;                               //Maintenance respiration coefficient of stems 
					rm31=0.00448;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.004992;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.2;                                  //ratio of fine root and leaf 
					crstem= 0.22;                                  //ratio of coarse root and stem
					Q10=2.0;                                      //Q10 is set to 2.0
					
					pafRmH[i] = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
				}

				//Evergreen Broadleaf Forests 
				if(pafLandcoverH[i]==2)
				{
					SLA=32.0;                                 //revised by MJ Wang 2019.4.16
					stema = -1.0378;                         // The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 31.298;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					rm1=0.011776;                              //Maintenance respiration coefficient of the leaves
					rm2=0.004736;                               //Maintenance respiration coefficient of stems 
					rm31=0.004736;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.00576;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.2;                                  //ratio of fine root and leaf 
					crstem= 0.22;                                  //ratio of coarse root and stem
					Q10=2.0;                                      //Q10 is set to 2.0
					
					pafRmH[i] = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
				}


				//revised by MJ Wang
				//mixed forest = mean(DBF + ENF)revised by MJ Wang 2019.4.16//2019.5.18
				if(pafLandcoverH[i]==5)
				{
					SLA=23.5;                                  //SLA is the specific leaf area
					stema = -1.59;                         // The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 41.412;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					rm1=0.0088;                              //Maintenance respiration coefficient of the leaves
					rm2=0.00456;                               //Maintenance respiration coefficient of stems 
					rm31=0.00404;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.004146;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.3;                                  //ratio of fine root and leaf
					crstem= 0.255;                                  //ratio of coarse root and stem
					Q10=2.0;                                      //Q10 is set to 2.0

					pafRmH[i] = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
				}

				//Other forest (closed shrubland, open shrubland)
				if(pafLandcoverH[i]==6||pafLandcoverH[i]==7)
				{
					SLA=17.0;                                 //specific leaf area
					stema = -1.18272;                         //The coefficient of the quadratic term when livestem C is calculated using LAImax
					stemb = 10.0048;                         // The coefficient of the primary term when livestem C is calculated using LAImax
					rm1=0.0068;                              //Maintenance respiration coefficient of the leaves
					rm2=0.0038;                               //Maintenance respiration coefficient of stems 
					rm31=0.0034;                              //Maintenance respiration coefficient of coarse roots
					rm32=0.0033;                                 //Maintenance respiration coefficient of fine roots
					frleaf= 1.2;                                  //ratio of fine root and leaf
					crstem= 0.22;                                  //ratio of coarse root and stem
					Q10=2.0;                                      //Q10 is set to 2.0

					pafRmH[i] = coutRM(pafLaiH[i],pafLaimaxH[i],pafTavH[i],
						SLA,stema,stemb,rm1,rm2,rm31,rm32,frleaf, crstem, Q10,Lati_xrc,20);
				}

				//Calculate growth respiration Rg
				if(pafGppH[i] >= 65.534)
				{
					pafNppH[i]=pafGppH[i];
				}
				else
				{
					
					float Rg;
					
				    if(pafGppH[i]>pafRmH[i])
					{
						Rg=0.25*(pafGppH[i]-pafRmH[i]);
						float Ra=pafRmH[i]+Rg;
						pafNppH[i]=pafGppH[i]-Ra;
					}
					else
					{ 
						pafNppH[i]=0;  
					}
				}
			}
		}
}




//Maintenance respiration is calculated based on land use type
//revised by Mengjia

float coutRM(double pafLAI,unsigned int pafLAImaxy,float pafTcor,
	float SLA,float stema, float stemb, float rm1,float rm2,float rm31,float rm32,
	float frleaf, float crstem, float Q10,	int Lati_xrc,float Tb)//revised by zhl，GLASS 500m LAI SCALE-FACTOR=0.01
{

		//Calculation of leaf biomass gC/m2
	float M1 = 1000 * pafLAI*0.01 / SLA; //revised by zhl，MuSyQLAI SCALE - FACTOR = 1000, MUSES LAI: 100
		//Calculate the maintenance respiration of the leaves	gC/(m2*d)	
		float leafmr=M1*rm1*pow(Q10,(pafTcor-Tb)/10);
		//LeafMR is non-negative
		leafmr=leafmr<0?0:leafmr;

		//Calculation of stem biomass was  gC/m2
		float M2 = stema*pafLAImaxy*0.01*pafLAImaxy*0.01 + stemb*pafLAImaxy*0.01; //revised by zhl，MuSyQLAI SCALE - FACTOR = 1000 MUSES LAI为100
		//Calculation the maintenance respiration of stem  gC/(m2*d)
		float stemmr=M2*rm2*pow(Q10,(pafTcor-Tb)/10);	
		//stemmris non-negative
		stemmr=stemmr<0?0:stemmr;

		//Calculate coarse root biomass gC/m2
		float M31 = crstem*M2;
		//Calculate fine root biomass gC/m2
		float M32 = frleaf*M1;
		
		//Calculate the maintenance respiration of the coarse roots gC/(m2*d)
		float crootmr=M31*rm31*pow(Q10,(pafTcor-Tb)/10);    
		crootmr=crootmr<0?0:crootmr;
		//Calculate the maintenance respiration of the fine roots gC/(m2*d)
		float frootmr=M32*rm32*pow(Q10,(pafTcor-Tb)/10);     
		frootmr=frootmr<0?0:frootmr;
		//Total maintenance respiration 
		float result=(crootmr+frootmr+stemmr+leafmr);
		return result;
}


