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
#include <gdal.h>
//#include <netcdf>
  
using namespace std;

/** 
* Preprocessing of input data, including climate data preprocessing and TIFF data preprocessing
**/
string GetStartDay(string fparPath)
{
	string startDay = "";//Define the start days to be returned
	int n=0;
	int ncRecount=0;
	int count=0;	

	//A character vector used to store the path of the PAR file
	string GLDASfilename =fparPath;
	int GLDASfilesize = GLDASfilename.size();

	for(int j=0;j<GLDASfilesize;j++)
	{
		if(GLDASfilename[j]=='.')
			count+=1;
		if(count==3)
		{
			n=j;
			break;
		}
	}
	int date_match;
	date_match=(GLDASfilename[n+1]-48)*1000000+(GLDASfilename[n+2]-48)*100000+(GLDASfilename[n+3]-48)*10000+(GLDASfilename[n+4]-48)*1000+(GLDASfilename[n+5]-48)*100+(GLDASfilename[n+6]-48)*10+1;
   
	//Convert an integer to a string
    char t[256];
    sprintf(t, "%d", date_match);
    startDay = t;

	return startDay;
}


extern bool IsRSParExist;
extern bool IsRSSMExist;


void CreateLL(string inPath,string outPath)
{
	GDALDataset *SMDataset;
	//The driver used to create a new file
	GDALDriver *poDriver;
	GDALDataset *pDataSet = (GDALDataset *) GDALOpen(inPath.c_str(), GA_ReadOnly );	
	int Lati_xrc=pDataSet->GetRasterXSize();
	int Lati_yrc=pDataSet->GetRasterYSize();
	double adfGeoTransform[6];
	const char * ProjectionStr;
	pDataSet->GetGeoTransform(adfGeoTransform);
	// Gets the coordinate system of the source image 
	const char *pszSrcWKT = NULL;
	char*  pszDstWKT = NULL;
    //pszSrcWKT = ProjectionStr;
	pszSrcWKT=GDALGetProjectionRef(pDataSet);
    CPLAssert( pszSrcWKT != NULL &&strlen(pszSrcWKT) > 0 );	
	const char *fomat="GTiff";
	poDriver = GetGDALDriverManager()->GetDriverByName(fomat);
	char **papszMetadata = poDriver->GetMetadata();
	//string resultPath=outPath.append("latitude.tif");
	string resultPath=outPath;
	SMDataset=poDriver->Create(resultPath.c_str(),Lati_xrc,Lati_yrc,2,GDT_Float32,NULL);
	SMDataset->SetGeoTransform(adfGeoTransform);
	SMDataset->SetProjection(pszSrcWKT);

	double* Lati=(double*)CPLMalloc(sizeof(double)*Lati_xrc);
	double* Long=(double*)CPLMalloc(sizeof(double)*Lati_xrc);
	//The latitude and longitude data of double is converted to the latitude and longitude data of float
	float* Lati_temp=(float*)CPLMalloc(sizeof(float)*Lati_xrc);
	float* Long_temp=(float*)CPLMalloc(sizeof(float)*Lati_xrc);

	for (int m = 0; m < Lati_yrc; m++)
	{
		//Geographic coordinates convert latitude and longitude
		OGRCoordinateTransformation *poTransform; 
		OGRSpatialReference oSRS,*poLatLong;
		char* pc=new char[512];
		int x=strlen(pszSrcWKT);
		strcpy(pc,pszSrcWKT);
		char* te=pc;
		oSRS.importFromWkt(&pc);
		delete[] te;

		poLatLong=oSRS.CloneGeogCS();

		poTransform = OGRCreateCoordinateTransformation( &oSRS,poLatLong ); 
		if(poTransform==NULL)
		{
			printf("Can't build poTransform ");
		}

		for(int i=0;i<Lati_xrc;i++)
		{
			ImageRowCol2Projection(adfGeoTransform,i,m,
				Long[i],Lati[i]);

		}
		poTransform-> Transform(Lati_xrc,Long,Lati);

		for(int i=0;i<Lati_xrc;i++)
		{
			Lati_temp[i] = (float)Lati[i];
			Long_temp[i] = (float)Long[i];
		}

		SMDataset->GetRasterBand(1)->RasterIO(GF_Write,0,m,Lati_xrc,1,Lati_temp,Lati_xrc,1,GDT_Float32,0,0);
		SMDataset->GetRasterBand(2)->RasterIO(GF_Write,0,m,Lati_xrc,1,Long_temp,Lati_xrc,1,GDT_Float32,0,0);

	}

	CPLFree(Lati);
	CPLFree(Long);
	CPLFree(Lati_temp);
	CPLFree(Long_temp);
	GDALClose(SMDataset);
	GDALClose(pDataSet);

}

void SoilMoisture(float* pafsoilM1,int xrc,int yrc)
{   
	float** a = new float*[yrc];//Create an array of pointers
	for(int i = 0;i<yrc;i++)
	{
		a[i] = new float[xrc];//Allocate space for each row
	}

	for(int i=0;i<yrc;i++)
	{
		for(int j=0;j<xrc;j++)
			a[i][j]=pafsoilM1[i*xrc+j];
	}
	//The first row of data is processed, starting from the second pixel. The first one and the last one are dealt separately. 
	//The first pixel is processed
	int *P = new int[5];        //Parameter operator, which is used to count the number of anomalies in the window
	for(int m=0;m<5;m++)
      P[m]=1;
   if(a[0][0]>NoData) //In the first pixel treatment, due to the abnormally large value of the water area, it is considered to be water area greater than NoData. The threshold section is adjusted according to the range of values taken in different regions.
   {
	   int num=0;
	   if(a[0][1]>NoData)
	   {
			num+=1;
			P[0]=0;
	   }
	   if(a[1][0]>NoData)
       {
			num+=1;
			P[1]=0;
	   }
	   if(a[1][1]>NoData)
       {
			num+=1;
			P[2]=0;	   
	   }
	   if(num==3)
		pafsoilM1[0] = a[0][0];
	   else
		pafsoilM1[0]=(a[0][1]*P[0]+a[1][0]*P[1]+a[1][1]*P[2])/(3-num);
	}
	else
		pafsoilM1[0] = a[0][0];

	//Intermediate pixel processing
    for (int j = 1; j < xrc-1; j++)
	{
		for(int m=0;m<5;m++)
        P[m]=1;
 
		int num=0;

		if(a[0][j]>NoData)
		{	
	        if(a[0][j-1]>NoData)
			{
				num+=1;
				P[0]=0;
			}
	        if(a[0][j+1]>NoData)
			{
				num+=1;
				P[1]=0;
			}
	        if(a[1][j-1]>NoData)
			{
				num+=1;
				P[2]=0;
			}
	        if(a[1][j]>NoData)
			{
				num+=1;
				P[3]=0;
			}
	        if(a[1][j+1]>NoData)
			{
				num+=1;
				P[4]=0;
			}
			if(num==5)
				pafsoilM1[j] = a[0][j];
			else
               pafsoilM1[j]=(a[0][j]*P[0]+a[0][j+1]*P[1]+a[1][j-1]*P[2]+a[1][j]*P[3]+a[1][j+1]*P[4])/(5-num);
			//   pafThreeLineWin[j]=pafOutputBuf[j];
		}
        else
        pafsoilM1[j] = a[0][j];  
	}

	//The last pixel of the first row is processed
	for(int m=0;m<3;m++)
		P[m]=1; 
	if(a[0][xrc-1]>NoData) 
    {
		int num=0;
	    if(a[0][xrc-2]>NoData)
		{
			num+=1;
		    P[0]=0;
		}
		if(a[1][xrc-2]>NoData)
		{
			num+=1;
			P[1]=0;
		}
		if(a[1][xrc-1]>NoData)
		{
			num+=1;
			P[2]=0;
		}
		if(num==3)
			pafsoilM1[xrc-1] = a[0][xrc-1];
		else
            pafsoilM1[xrc-1]=(a[0][xrc-2]*P[0]+a[1][xrc-2]*P[1]+a[1][xrc-1]*P[2])/(3-num);
	}
   else
	pafsoilM1[xrc-1] = a[0][xrc-1];

	for (int i = 1; i < yrc-1; i++)  
   {  
		//Process the first column of data
		for(int m=0;m<5;m++)    //Template operator
             P[m]=1;
		float b=a[2][0];
  
		if(a[i][0]>NoData)
		{
			int num=0;
	        if(a[i-1][0]>NoData)
			{
				num+=1;
				P[0]=0;   
			}
	        if(a[i-1][1]>NoData)
			{
				num+=1;
				P[1]=0;
			}
	        if(a[i][1]>NoData)
			{
				num+=1;
				P[2]=0;
			}
	        if(a[i+1][0]>NoData)
			{
				num+=1;
				P[3]=0;
			}
	        if(a[i+1][1]>NoData)
			{
				num+=1;
				P[4]=0;
			}
			if(num==5)
				pafsoilM1[xrc*i] = a[i][0];
			else

               pafsoilM1[xrc*i]=(a[i][0]*P[0]+a[i][1]*P[1]+a[i][1]*P[2]+a[i+1][0]*P[3]+a[i+1][1]*P[4])/(5-num);
		}
        else
        pafsoilM1[xrc*i] = a[i][0];
		//deal intermediate data
		for (int j = 1; j < xrc - 1; j++)  
        {  
	       float afWin[9];  
          afWin[0] = a[i-1][j-1];  
          afWin[1] = a[i-1][j];  
          afWin[2] = a[i-1][j+1];  
          afWin[3] = a[i][j-1];  
          afWin[4] = a[i][j];  
          afWin[5] = a[i][j+1];   
          afWin[6] = a[i+1][j-1]; 
          afWin[7] = a[i+1][j];  
          afWin[8] = a[i+1][j+1]; 
		   int count=0;    //Count the number of outliers in the sliding window
           for(int m=0; m<9; m++)
		   
			    if(afWin[m]>NoData)
			       count+=1;
		   
		   if(count==9)    //If all nine pixels are outliers, the central pixel value remains the same, and it is still an outlier
			   pafsoilM1[xrc*i+j]=a[i][j];
			//Otherwise, interpolation is carried out
		   else            
		   {
            if (afWin[0]>NoData)
				afWin[0]=0.0;
			if (afWin[1] >NoData)
				afWin[1]=0.0;
		    if (afWin[2] >NoData)
				afWin[2]=0.0;
	        if (afWin[3] >NoData)
				afWin[3]=0.0;
			if (afWin[5]>NoData)
				afWin[5]=0.0;
			if (afWin[6]>NoData)
				afWin[6]=0.0;
			if (afWin[7]>NoData)
				afWin[7]=0.0;
			if (afWin[8]>NoData)
				afWin[8]=0.0;
		    if (afWin[4]>NoData)
			{
			    afWin[4]=0.0;
               //Water value   
               pafsoilM1[xrc*i+j] =(afWin[0]+afWin[1]+afWin[2]+afWin[3]+afWin[5]+afWin[6]+afWin[7]+afWin[8])/(9.0-count) ;
			}
			else
				pafsoilM1[xrc*i+j]=a[i][j];
		   }
     }
	 //Process the last column of data

		for(int m=0;m<5;m++)    //Template operator
             P[m]=1;
  
		if(a[i][xrc - 1]>NoData)
		{
			int num=0;
	        if(a[i-1][xrc - 1]>NoData)
			{
				num+=1;
				P[0]=0;   
			}
	        if(a[i-1][xrc - 2]>NoData)
			{
				num+=1;
				P[1]=0;	 
			}
	        if(a[i][xrc - 2]>NoData)
			{
				num+=1;
				P[2]=0;	 
			}
	        if(a[i+1][xrc - 2]>NoData)
			{
				num+=1;
				P[3]=0;	
			}
	        if(a[i+1][xrc - 1]>NoData)
			{
				num+=1;
				P[4]=0;	 
			}
			if(num==5)
				pafsoilM1[xrc*i+xrc-1] = a[i][xrc - 1];
			else

               pafsoilM1[xrc*i+xrc-1]=(a[i-1][xrc - 1]*P[0]+a[i-1][xrc - 2]*P[1]+a[i][xrc - 2]*P[2]+a[i+1][xrc - 2]*P[3]+a[i+1][xrc - 1]*P[4])/(5-num);
			//   pafThreeLineWin[j]=pafOutputBuf[j];
			
		}
        else
        pafsoilM1[xrc*i+xrc-1] = a[i][xrc - 1]; 
	   }
		//The last row of data is processed, similar to the first row of data processing
	   if (yrc > 1)  
	   {  

			for(int m=0;m<5;m++)
			P[m]=1;

			//The first pixel is processed
			if(a[yrc-1][0]>NoData) 
			{
				int num=0;
				if(a[yrc-1][1]>NoData)
				{
					num+=1;
					P[0]=0;
				}
				if(a[yrc-2][0]>NoData)
				{
					num+=1;
					P[1]=0;
				}
				if(a[yrc-2][1]>NoData)
				{
					num+=1;
					P[2]=0;
				}
				if(num==3)
					pafsoilM1[xrc*(yrc-1)] = a[yrc-1][0];
				else
					pafsoilM1[xrc*(yrc-1)]=(a[yrc-1][1]*P[0]+a[yrc-2][0]*P[1]+a[yrc-2][1]*P[2])/(3-num);
			}
			else
				pafsoilM1[xrc*(yrc-1)] = a[yrc-1][0];

			//Intermediate pixel processing
			for (int j = 1; j < xrc-1; j++)
			{
				//Parameter operator
				for(int m=0;m<5;m++)
				P[m]=1;
				int num=0;
				if(a[yrc-1][j]>NoData)
				{
					if(a[yrc-1][j-1]>NoData)
					{
						num+=1;
						P[0]=0;
					}
					if(a[yrc-1][j+1]>NoData)
					{
						num+=1;
						P[1]=0;
					}
					if(a[yrc-2][j-1]>NoData)
					{
						num+=1;
						P[2]=0;
					}
					if(a[yrc-2][j]>NoData)
					{
						num+=1;
						P[3]=0;
					}
					if(a[yrc-2][j+1]>NoData)
					{
						num+=1;
						P[4]=0;
					}
					if(num==5)
						pafsoilM1[xrc*(yrc-1)+j] = a[yrc-1][j];
					else
						pafsoilM1[xrc*(yrc-1)+j]=(a[yrc-1][j-1]*P[0]+a[yrc-1][j+1]*P[1]+a[yrc-2][j-1]*P[2]+a[yrc-2][j]*P[3]+a[yrc-2][j+1]*P[4])/(5-num);
				}
				else
					pafsoilM1[xrc*(yrc-1)+j] = a[yrc-1][j];  
			}

			//The last pixel is processed
			for(int m=0;m<3;m++)
				P[m]=1; 
				if(a[yrc-1][xrc-1]>NoData) 
			    {
					int num=0;
					if(a[yrc-1][xrc-2]>NoData)
					{
						num+=1;
						P[0]=0;
					}
					if(a[yrc-2][xrc-2]>NoData)
					{
						num+=1;
						P[1]=0;
					}
					if(a[yrc-2][xrc-1]>NoData)
					{
						num+=1;
						P[2]=0;
					}
					if(num==3)
						pafsoilM1[xrc*yrc-1] = a[yrc-1][xrc-1];
					else
						pafsoilM1[xrc*yrc-1]=(a[yrc-1][xrc-2]*P[0]+a[yrc-2][xrc-2]*P[1]+a[yrc-2][xrc-1]*P[2])/(3-num);
				}
				else
					pafsoilM1[xrc*yrc-1] = a[yrc-1][xrc-1];

				delete []P;
				
				for(int m =0;m<yrc;++m)   
				delete[]   a[m];   
				delete[]   a; 
			}

		}


