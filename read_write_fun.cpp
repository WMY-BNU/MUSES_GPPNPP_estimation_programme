#include "stdafx.h"
#include "IMGALG.h"
#include <gdal_priv.h>
#include<gdal.h>
#include<ogrsf_frmts.h>
#include <iostream>
#include <fstream>
#include <cpl_conv.h>
#include <string.h>
#include "stdlib.h"
#include"ogr_spatialref.h"
#include "npp_fun.h"
#include "init_fun.h"
#include "hdf5.h"
#define FILE "file.h5"
  
using namespace std;
int num=1;
/** 
* Reads and writes data, including read/write HDF and read/write TIF
**/

//string filetype = ".nc";
vector<string> readHDF(const char* filePath,int* index,int size)
{
	vector<string> vecNULL;
	GDALAllRegister();
	GDALDataset *pDataSet = (GDALDataset *) GDALOpen( filePath, GA_ReadOnly );
	if (pDataSet == NULL)
	{
		printf("Cannot open the file %s , please check the file",filePath);
		getchar();
		return vecNULL;
	}

	char ** papszSUBDATASETS = GDALGetMetadata( (GDALDatasetH)pDataSet, "SUBDATASETS");
	vector<string> vSubDataSets;
	vector<string> vSubDataDesc;
	if ( papszSUBDATASETS == NULL )
	{
		string papszMetadata = GDALGetDriverShortName((GDALDriverH)pDataSet);
		vSubDataSets.push_back(papszMetadata);
		vSubDataDesc.push_back(papszMetadata);
	}
	else
	{
		int iCount = CSLCount(papszSUBDATASETS);
		if(  iCount <= 0 )
		{
			GDALClose((GDALDriverH)pDataSet);
			return vecNULL;
		}
		for(int i=0;i<size;i++)
		{
			string tmpstr = string(papszSUBDATASETS[index[i]]);
			tmpstr = tmpstr.substr(tmpstr.find_first_of("=") + 1);
			const char *tmpc_str = tmpstr.c_str();
			string tmpdsc = string(papszSUBDATASETS[index[i]+1]);
			tmpdsc = tmpdsc.substr(tmpdsc.find_first_of("=") + 1);
			GDALDataset  *hTmpDt = (GDALDataset*)GDALOpen(tmpc_str, GA_ReadOnly);

			if(hTmpDt != NULL)
			{
				vSubDataSets.push_back(tmpstr);
				vSubDataDesc.push_back(tmpdsc);

				GDALClose(hTmpDt);
			}
		}

	}
	if(!vSubDataSets.empty())
	{
		GDALClose((GDALDriverH)pDataSet);
		return vSubDataSets;
	}
	else
	{
		GDALClose((GDALDriverH)pDataSet);
		return vecNULL;
	}
}
float* getRasterdata(const char* path,int* xrc,int* yrc, char* tem_str,double* adfGeoTransform)
{
	GDALDataset  *hTmpDt = (GDALDataset*)GDALOpen(path, GA_ReadOnly);
	*xrc=hTmpDt->GetRasterXSize();
	*yrc=hTmpDt->GetRasterYSize();
	float* pafdata=(float*)CPLMalloc(sizeof(float)**xrc**yrc);
	GDALRasterBand* rasterBand=hTmpDt->GetRasterBand(1);
	rasterBand->RasterIO(GF_Read,0,0,*xrc,*yrc,pafdata,*xrc,*yrc,GDT_Float32,0,0);
	hTmpDt->GetGeoTransform(adfGeoTransform);
	const char*projectionref=hTmpDt->GetProjectionRef();
	strcpy(tem_str,projectionref);
	GDALClose(hTmpDt);
	return pafdata;
}
void save2TIFF(const char* savePath,int imgsizeX,int imgsizeY,int imgsizeBand,float* pafData,
	double* adfGeoTransform)
{
		GDALDataset *poDataset;
		//The driver used to create a new file
		GDALDriver *poDriver;
		//Create a new file
		const char *fomat="GTiff";
		poDriver = GetGDALDriverManager()->GetDriverByName(fomat);
		//Get the file format type
		char **papszMetadata = poDriver->GetMetadata();
		poDataset=poDriver->Create(savePath,imgsizeX,imgsizeY,1,GDT_Float32,papszMetadata);
		poDataset->SetGeoTransform(adfGeoTransform);
		//Write the processed image to the dataset
		poDataset->RasterIO(GF_Write, 0,0,imgsizeX,imgsizeY, pafData,imgsizeX,imgsizeY,GDT_Float32,1,0,0,0,0);
		GDALClose(poDataset);
}
//The entire image is imported into the result
void save2TIFF(const char* savePath,int imgsizeX,int imgsizeY,int imgsizeBand,float* pafData,
	double* adfGeoTransform,const char* projectionref)
{
		GDALDataset *poDataset;
		//The driver used to create a new file
		GDALDriver *poDriver;
		//Create a new file
		const char *fomat="GTiff";
		poDriver = GetGDALDriverManager()->GetDriverByName(fomat);
		//Get the file format type
		char **papszMetadata = poDriver->GetMetadata();
		poDataset=poDriver->Create(savePath,imgsizeX,imgsizeY,1,GDT_Float32,papszMetadata);
		poDataset->SetGeoTransform(adfGeoTransform);
		poDataset->SetProjection(projectionref);
		//Write the processed image to the dataset
		poDataset->RasterIO(GF_Write, 0,0,imgsizeX,imgsizeY, pafData,imgsizeX,imgsizeY,GDT_Float32,1,0,0,0,0);
		GDALClose(poDataset);
}
//The final result is saved in TIF format and imported into the result line by line
void save2TIFF(GDALDataset *poDataset, const char* savePath,int imgsizeX,int imgsizeY,int imgsizeBand,GDALDataType dataType,
	double* adfGeoTransform,const char* projectionref)
{
		
		//The driver used to create a new file
		GDALDriver *poDriver;
		//Create a new file
		const char *fomat="GTiff";
		poDriver = GetGDALDriverManager()->GetDriverByName(fomat);
		//Get the file format type
		char **papszMetadata = poDriver->GetMetadata();
		poDataset=poDriver->Create(savePath,imgsizeX,imgsizeY,imgsizeBand,dataType,papszMetadata);
		poDataset->SetGeoTransform(adfGeoTransform);
		poDataset->SetProjection(projectionref);
		//Write the processed image to the dataset
		//poDataset->RasterIO(GF_Write, 0,0,imgsizeX,imgsizeY, pafData,imgsizeX,imgsizeY,GDT_Float32,1,0,0,0,0);
		
}
void makeHDF(char* hdfchar,const char* GPPpath,const char* NPPpath,char* AcquisitionTime,char* ProductionTime,
             char* RawDataName,char* size,char* AlgorithmName,char* GridNum,char* QPProductName,char* NumBand,char* Projection,
             char* ProjectionStr,char* ProjectionPara)
{
	GDALAllRegister();
    //read GPP\NPP tiff file
    GDALDataset  *GPP_DT = (GDALDataset*)GDALOpen(GPPpath, GA_ReadOnly);
    GDALDataset  *NPP_DT = (GDALDataset*)GDALOpen(NPPpath, GA_ReadOnly);
	int xrc=GPP_DT->GetRasterXSize();//The number of rows of the image read
	int yrc=GPP_DT->GetRasterYSize();//The number of columns of the image read
    int* GPPdata=(int*)CPLMalloc(sizeof(int)*xrc*yrc);//Reads GPP data into the cache space
    int* NPPdata=(int*)CPLMalloc(sizeof(int)*xrc*yrc);//Reads NPP data into the cache space
    initPoint(GPPdata,xrc*yrc,0);
    initPoint(NPPdata,xrc*yrc,0);
    GDALRasterBand* rasterBand=GPP_DT->GetRasterBand(1);//Read GPP_DT band data
	GDALDataType dataType =rasterBand->GetRasterDataType();//Test the type of data
	rasterBand->RasterIO(GF_Read,0,0,xrc,yrc,GPPdata,xrc,yrc,dataType,0,0);
	rasterBand=NPP_DT->GetRasterBand(1);
	dataType =rasterBand->GetRasterDataType();
	rasterBand->RasterIO(GF_Read,0,0,xrc,yrc,NPPdata,xrc,yrc,dataType,0,0);

	char *CImageData1=NULL;
	if(ReadTif2Mem(GPPpath,hdfchar,CImageData1)==-1)
		printf("File read error!\n");

	char *CImageData2=NULL;
	if(ReadTif2Mem(NPPpath,hdfchar,CImageData2)==-1)
		printf("File read error!\n");

	delete [] CImageData1;
	CImageData1=NULL;
	delete [] CImageData2;
	CImageData2=NULL;
}
int ReadTif2Mem(const char* c_tif_file,const char* save_file,char *para)
{
	GDALDataset *PDataTIF = (GDALDataset *) GDALOpen( c_tif_file, GA_ReadOnly);

    if(PDataTIF==NULL)
		return -1;
	double adfGeoTransform[6];
	PDataTIF->GetGeoTransform(adfGeoTransform);	//Gets some affine parameters for the angle file
	
	GDALRasterBand *p_ImageBand=PDataTIF->GetRasterBand(1);
	int Hline=p_ImageBand->GetYSize();
	int Hsample=p_ImageBand->GetXSize();
	
	GDALDataType mo_type=p_ImageBand->GetRasterDataType();		//Data type: The test data is of the BYTE type. If the processing data is of a different type, the following request memory will be modified

	if(mo_type==GDT_Byte)
	{
		para=new char[Hline*Hsample];
		memset(para,0,sizeof(char)*Hline*Hsample);
		p_ImageBand->RasterIO(GF_Read,0,0,Hsample,Hline,para,Hsample,Hline,GDT_Byte,0,0);
	}
	else
	{
		return -1;
		printf("Please change the file data type!\n");
	}
	WriteTif2H5(save_file,para,Hline,Hsample);
	GDALClose(PDataTIF);

	return 0;
}
//Write the read TIF memory to the H5 file
void WriteTif2H5(const char* savepath,char*DATA1,int line,int sample) 
{

   hid_t       file_id;  /* identifiers */
   herr_t      status;
   hid_t       group_id, dataset_id, dataspace_id;  /* identifiers */
   hsize_t     dims[2];

   //Create an H5 file

   file_id = H5Fcreate(savepath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   //The group name of the H5 file is DATA
   group_id = H5Gcreate1(file_id, "/DATA", 0);

  // group2_id = H5Gcreate1(file_id, "/MyGroup/Group_A", 0);

   //Open the created H5 file and write data
   file_id = H5Fopen(savepath, H5F_ACC_RDWR, H5P_DEFAULT);

   //The size of the H5 file created
   dims[0] = line;
   dims[1] = sample;

   dataspace_id = H5Screate_simple(2, dims, NULL);

   //The dataset name of the file
    char *F_id="/DATA/dset_";
	char str[256];
	sprintf(str,"%s%d",F_id,num);
	num+=1;

	//Create a dataset file H5T_STD_I8BE byte type, and modify it to the corresponding code if it is written to a different type.
   dataset_id = H5Dcreate1(file_id,str , H5T_STD_I8BE, dataspace_id,
                          H5P_DEFAULT);

   //If the data is written to bytes, the H5T_NATIVE_CHAR is byte, and if the data is written to other types, you need to change the code to the corresponding type.
   H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     DATA1);

   //Close the dataset
    H5Sclose(dataspace_id);

   /* Close the first dataset. */
    H5Dclose(dataset_id);

   /* Close the group. */
   H5Gclose(group_id);

   /* Close the file. */
   H5Fclose(file_id);
}

