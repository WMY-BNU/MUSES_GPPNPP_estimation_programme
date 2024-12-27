#pragma once
#include "IMGALG.h"
#include <gdal_priv.h>
#include <iostream>
#include <fstream>
#include <cpl_conv.h>
#include <string.h>
#include "stdlib.h"
#include"ogr_spatialref.h"

using namespace std;

void initPoint(float* paf,int count);
/**
* @brief                             Calculate the distance between two points
* @param x1,x2,y1,y2                 Coordinates of two points 
*/
float Distance(float x1,float x2,float y1,float y2);
/**
* @brief                             Nearest neighbor method resampling
* @return                            Returns the float value of the sampling point
* @param x_low_index,y_low_index     Index of the upper left corner of the area to be sampled on a 0.25бу image
* @param x_low,y_low                 coordinate of the upper left corner of the area to be sampled on a 0.25бу image
* @param x,y                         The coordinates of the point to be sampled
* @param data                        0.25бу image data
* @param yrc                         The number of rows in a 0.25бу image
* @param flag                        Determine whether the sampling point coincides with the coordinate point of the 0.25бу image, 1 is coincident, and 0 is not coincident
* @param pix_size                    The resolution of the 0.25бу image, here is 0.25 degree
*/
float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	float* data,int yrc,int flag,float pix_size);
//Overloading for int type data
float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	int* data,int yrc,int flag,float pix_size);

//Overloading for DT_16U type data
float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_16U* data,int yrc,int flag,float pix_size);
//Overloading for DT_32U type data
float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_32U* data,int yrc,int flag,float pix_size);

//Overloading for DT_8U type data
float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_8U* data,int yrc,int flag,float pix_size);
/**
* @brief                             Bilinear resampling
* @return                            Returns the float value of the sampling point
* @param x_low_index,y_low_index     Index of the upper left corner of the area to be sampled on a 0.25бу image
* @param x_low,y_low                 Coordinate of the upper left corner of the area to be sampled on a 0.25бу image
* @param x,y                         The coordinates of the point to be sampled
* @param data                        0.25бу image data (float)
* @param yrc                         The number of rows in a 0.25бу image
* @param flag                        Determine whether the sampling point coincides with the coordinate point of the 0.25бу image, 1 is coincident, and 0 is not coincident
* @param pix_size                    The resolution of the 0.25бу image, here is 0.25 degree
* @param nodata                      Invalid value for the image
*/
float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	float* data,int yrc,int flag,float pix_size);
float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	float* data,int yrc,int flag,float pix_size,float nodata);
/**
* @brief                             Bilinear resampling
* @return                            Returns the float value of the sampling point
* @param x_low_index,y_low_index     Index of the upper left corner of the area to be sampled on a 0.25бу image
* @param x_low,y_low                 Coordinate of the upper left corner of the area to be sampled on a 0.25бу image
* @param x,y                         The coordinates of the point to be sampled
* @param data                        0.25бу image data (DT_16U)
* @param yrc                         The number of rows in a 0.25бу image
* @param flag                        Determine whether the sampling point coincides with the coordinate point of the 0.25бу image, 1 is coincident, and 0 is not coincident
* @param pix_size                    The resolution of the 0.25бу image, here is 0.25 degree
*/
float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_16U* data,int yrc,int flag,float pix_size);
float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_32U* data,int yrc,int flag,float pix_size);


/**
* @brief                             Bilinear resampling
* @return                            Returns the float value of the sampling point
* @param x_low_index,y_low_index     Index of the upper left corner of the area to be sampled on a 0.25бу image
* @param x_low,y_low                 Coordinate of the upper left corner of the area to be sampled on a 0.25бу image
* @param x,y                         The coordinates of the point to be sampled
* @param data                        0.25бу image data (DT_8U)
* @param yrc                         The number of rows in a 0.25бу image
* @param flag                        Determine whether the sampling point coincides with the coordinate point of the 0.25бу image, 1 is coincident, and 0 is not coincident
* @param pix_size                    The resolution of the 0.25бу image, here is 0.25 degree
*/
float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_8U* data,int yrc,int flag,float pix_size);


/**
* @brief              Read HDF
* @return             Returns a vector that contains the path where the specified subset of data is located
* @param filePath     HDF file path
* @param index        Subdataset index
* @param size         Number of sub-datasets
*/
vector<string> readHDF(const char* filePath,int* index,int size);
/**
* @brief              Reads a single-band file and returns a pointer of the stored data
* @return             No returns
* @param path         File path
* @param xrc          The number of columns of the data
* @param yrc          Number of rows of data
* @param projectionref          Projection information of the data
* @param adfGeoTransform        Spatial reference of the data
*/
float* getRasterdata(const char* path,int* xrc,int* yrc,char*tem_str,double* adfGeoTransform);
/**
* @brief                       Save data to a TIFF file
* @return                      No returns
* @param savePath              Save the path
* @param imgsizeX              The number of columns where the data is saved
* @param imgsizeY              The number of rows where the data is saved
* @param imgsizeBand           The number of bands where the data is saved
* @param pafData               A pointer that contains data
* @param adfGeoTransform       Spatial coordinate information
*/
void save2TIFF(const char* savePath,int imgsizeX,int imgsizeY,int imgsizeBand,float* pafData,double* adfGeoTransform);

/**
* @brief                        Save data to a TIFF file
* @return                      No returns
* @param savePath              Save the path
* @param imgsizeX              The number of columns where the data is saved
* @param imgsizeY              The number of rows where the data is saved
* @param imgsizeBand           The number of bands where the data is saved
* @param pafData               A pointer that contains data
* @param adfGeoTransform       Spatial coordinate information
* @param projectionref          Projection information of the data
*/
void save2TIFF(const char* savePath,int imgsizeX,int imgsizeY,int imgsizeBand,float* pafData,
	double* adfGeoTransform,const char* projectionref);

/**
* @brief                  Soil moisture interpolation function
* @return                 No returns
* @param pafsoilM1        The float pointer to be interpolated
* @param xrc              Number of rows of the data
* @param yrc              Number of columns of the data
*/
void SoilMoisture(float* pafsoilM1,int xrc,int yrc);

/**
* @brief                  Initialize the pointer
* @return                 No returns
* @param paf              Float pointer
* @param count            Pointer size
* @param value            Initialized values
*/
void initPoint(double* paf,int count,double value);
void initPoint(float* paf,int count,float value);
void initPoint(DT_8U* paf,int count,DT_8U value);
void initPoint(int* paf,int count,int value);
void initPoint(DT_16U* paf,int count,double value);
void initPoint(DT_8U* paf,int count);
void initPoint(DT_32U* paf,int count,double value);
/**
* @brief                  In Windows, get the list of files with the specified suffix in the directory
* @return                 No returns
* @param foldname         The name of the path
* @param filelist         The vector that stores the file path
* @param filetype		  The type of file
*/
void get_filelist(const char *dir,vector<string> &filelist,string filetype);

/**
* @brief                       Resample 0.25бу ERA5 data using latitude and longitude data
* @return                      No returns
* @param pafLong               A pointer that stores longitude
* @param pafLat                A pointer that stores latitude
* @param pix_size              The cell size of the ERA5 data
* @param Lup_long			   Longitude in the upper left corner of the ERA5 data
* @param Lup_lat               Latitude in the upper left corner of the ERA5 data
* @param pafnew                The pointer after resampling
* @param rsMath			       The resampling method, 1 is bilinear and 0 is the nearest neighbor
* @param Lxrc			       The number of columns for latitude and longitude data
* @param pafData			   Pointer to the data to be sampled(float)
* @param xrc			       The number of columns of ERA5 data
* @param nodata			       Invalid value
*/
void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,float* pafData,int xrc);
void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,float* pafData,int xrc,float nodata);

/**
* @brief                       Resample 0.25бу ERA5 data using latitude and longitude data
* @return                      No returns
* @param pafLong               A pointer that stores longitude
* @param pafLat                A pointer that stores latitude
* @param pix_size              The cell size of the ERA5 data
* @param Lup_long			   Longitude in the upper left corner of the ERA5 data
* @param Lup_lat               Latitude in the upper left corner of the ERA5 data
* @param pafnew                The pointer after resampling
* @param rsMath			       The resampling method, 1 is bilinear and 0 is the nearest neighbor
* @param Lxrc			       The number of columns for latitude and longitude data
* @param pafData			   Pointer to the data to be sampled(DT_16U)
* @param xrc			       The number of columns of ERA5 data
*/
void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,DT_16U* pafData,int xrc);

/**
* @brief                       Resample 0.25бу ERA5 data using latitude and longitude data
* @return                      No returns
* @param pafLong               A pointer that stores longitude
* @param pafLat                A pointer that stores latitude
* @param pix_size              The cell size of the ERA5 data
* @param Lup_long			   Longitude in the upper left corner of the ERA5 data
* @param Lup_lat               Latitude in the upper left corner of the ERA5 data
* @param pafnew                The pointer after resampling
* @param rsMath			       The resampling method, 1 is bilinear and 0 is the nearest neighbor
* @param Lxrc			       The number of columns for latitude and longitude data
* @param pafData			   Pointer to the data to be sampled(float)(DT_8U)
* @param xrc			       The number of columns of ERA5 data
*/
void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,DT_8U* pafData,int xrc);


void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,DT_32U* pafData,int xrc);

/**
* @brief                       Calculate the maintenance of respiration for leaves, stems, and roots
* @return                      No returns
* @param pafLAI                A pointer to store LAI data
* @param pafBiomass            A pointer to store biomass
* @param pafrm                 Return the calculated autotrophic respiration
* @param SLA-Q10			   Vegetation parameters
	parameters		unit
	SLA		m2/kgC
	rm1		kgC/day/kg
	rm2		kgC/day/kg
	rm31	kgC/day/kg
	rm32	kgC/day/kg
	Q10
* @param Lati_xrc              Number of data columns
* @param Tb                    base temperature

*/
float coutRM(double pafLAI,unsigned int pafLAImaxy,float pafTcor,
	float SLA,float stema, float stemb, float rm1,float rm2,float rm31,float rm32,
	float frleaf, float crstem, float Q10,	int Lati_xrc,float Tb);
/**
* @brief                       Output log information to a file
* @return                      No returns
* @param logopath              The path to the log file
* @param logopath              The contents of log outputs 
*/
void makeLogo(const char* logopath,const char* content);


/**
* @brief                       Get the start days from which the NPP calculation starts from the file name of the FPAR
* @return                      Returns the start days
* @param fparPath			   The path to the input FPAR file
*/
string GetStartDay(string fparPath);


/**
* @brief                       Converts row and column numbers to geographic coordinates
* @return                      bool returns True if successful, otherwise it returns False
* @param adfGeoTransform       Geographic information of the source data
* @param iCol		           Image row number
* @param iRow                  Image column number
* @param dProjX                The returned X coordinates
* @param dProjY                The returned Y coordinates
*/
bool ImageRowCol2Projection(double *adfGeoTransform, int iCol, int iRow, double &dProjX, double &dProjY);

/**
* @brief                       Convert geographic coordinates to row and column numbers
* @return                      bool returns True if successful, otherwise it returns False
* @param adfGeoTransform       Geographic information of the source data
* @param dProjX                The returned X coordinates
* @param dProjY                The returned Y coordinates
* @param iCol		           Image row number
* @param iRow                  Image column number
*/
bool Projection2ImageRow(double* adfGeoTransform,double dProjX,double dProjY,int&iCol,int&iRow);

/**
* @brief            Generate a latitude and longitude image of the image
* @return           No returns
* @param inPath		The path of the source image
* @param outPath    The path to the temporary directory
*/
void CreateLL(string inPath,string outPath);

/**
* @brief					Calculate LUE based on land use type
* @return					No returns
* @param LANDCOVER_DT		A land-use type image file pointer
* @param j					Line j of the the land use image
* @param Lati_xrc			The number of columns for the land-use image
* @param pafLandcoverM		A pointer to a row of data that stores land-use image
* @param pafLUE				The pointer to the LUE is stored
*/
void computerLUE(GDALDataset *LANDCOVER_DT, int j, int Lati_xrc, DT_8U* pafLandcoverM, float* pafLUE, float* CI, GDALDataset*FPAR_DT, GDALDataset*LAI_DT, DT_8U*pafFparH, DT_16U*pafLaiH);

//**
//* @brief							Calculate the temperature stress factor
//* @return							No returns
//* @param Lati_DT					Pointer to the latitude and longitude image file
//* @param j						Line j of the the land use image
//* @param Lati_xrc					The number of columns for the land-use image
//* @param GLDAS_xrc				The number of columns in the ERA5 image
//* @param DEMM_xrc					Number of columns for 500m DEM
//* @param GLDAS_GeoTransform		ERA5 image (0.25 degrees) spatial coordinate information
//* @param STRMM_DEM_GeoTransform	DEM image (500m) spatial coordinate information
//* @param pafLongH					The longitude pointer of a row of data in a latitude and longitude image (500m).
//* @param pafLatH					The latitude pointer of a row of data in a latitude and longitude image (500m).
//* @param pafTavH					Pointer to a row of average temperature data at 500 m resolution
//* @param pafTavL					Pointer to the entire average temperature image at 0.25 degree resolution
//* @param pafMatH					Pointer to a row of maximum temperature data at 500 m resolution
//* @param pafMatL					Pointer to the entire maximum temperature image at 0.25 degree resolution
//* @param pafToptH					Pointer to a row of optimum temperature data at 500 m resolution
//* @param pafToptL					Pointer to the entire optimum temperature image at 0.25 degree resolution
//* @param pafDemM					Data pointer for the entire frame of the DEM image at 500 meter resolution
//* @param pafDemH					A data pointer for a row of DEM images with a resolution of 500 meters
//* @param pafTE1H					A pointer to a row of data for TE1
//* @param pafK1H					A pointer to a row of data in K1
//* @param pafK2H					A pointer to a row of data in K2
//* @param pafTE2H					A pointer to a row of data in TE2
//*/

///revised by MJ Wang 2019.4.14

void comTempStress(GDALDataset *Lati_DT,GDALDataset *DEMM_DT, int j, int Lati_xrc, int GLDAS_xrc,DT_8U* pafLandcoverH,int DEMM_xrc, double GLDAS_GeoTransform[6],double STRMM_DEM_GeoTransform[6],
				   float* pafLongH, float* pafLatH, float* pafTavH, float* pafTavL,float* pafToptH,
				   float* pafDemM, float* pafDEMH, double* pafK1H, double* pafK2H, double* pafTE2H);


//**
//* @brief							Calculate the PAR
//* @return							No returns
//* @param Lati_DT					Pointer to the latitude and longitude image file
//* @param j						Line j of the the land use image
//* @param Lati_xrc					The number of columns for the land-use image
//* @param GLDAS_xrc				The number of columns in the ERA5 image
//* @param RsPar_xrc				Number of columns for remote sensing inversted PAR image (5km).
//* @param GLDAS_GeoTransform		ERA5 image (0.25 degrees) spatial coordinate information
//* @param RsPar_GeoTransform		Spatial coordinate information of PAR image (5km) retrieved by remote sensing
//* @param pafLongH					The longitude pointer of a row of data in a latitude and longitude image (500m)
//* @param pafLatH					The latitude pointer of a row of data in a latitude and longitude image (500m).
//* @param pafsGldasParL			Pointer to the entire PAR image at 0.25 degree resolution
//* @param pafGldasParH				Pointer to a row of PAR data at 500 m resolution
//* @param pafRsParML				Pointer to the entire remotely sensed PAR data at 5 km resolution
//* @param pafRsParH				Pointer to a row of remotely sensed PAR data at 500 m resolution
//* @param pafParH					Pointer to a row of PAR (500m) data that is finally used for GPP calculations
//*/
void comPar(GDALDataset *Lati_DT, int j, bool IsDataExist, int Lati_xrc, int GLDAS_xrc, int RsPar_xrc, double GLDAS_GeoTransform[6], double RsPar_GeoTransform[6],
	float* pafLongH, float* pafLatH, float* pafsGldasParL, float* pafGldasParH, float* pafRsParML, float* pafRsParH, float* pafParH, float nodata);

//**
//* @brief							Calculate GPP
//* @return							No returns
//* @param FPAR_DT					FPAR image file pointer
//* @param j						Line j of the the land use image
//* @param Lati_xrc					The number of columns for the land-use image
//* @param GLDAS_xrc				The number of columns in the ERA5 image
//* @param pafFparH					Pointer to a row of data in FPAR image (500M).
//* @param pafGppH					Pointer to a row of data in GPP image (500M).
//* @param pafLUEH					Pointer to a row of data in LUE image (500M).
//* @param paff2H					Pointer to a row of data in f2 image (500M).
//* @param pafTE1H					Pointer to a row of data in te1 image (500M).
//* @param pafTE2H					Pointer to a row of data in te2 image (500M).
//* @param pafParH					Pointer to a row of PAR (500m) data that is finally used for GPP calculations
//*/
void comGPP(GDALDataset *FPAR_DT, int j, int Lati_xrc, DT_8U* pafFparH, float FPARnodata, DT_8U* pafLandcoverH, double* pafGppH, float* pafLUEoutp,
			float* pafLUEH, float* paff2H, double* pafTE2H, float* pafParH);///revised by MJ Wang 20190601 in order to output LUE
//**
//* @brief							Calculate GPP
//* @return							No returns
//* @param LAI_DT					Pointer to the latitude and longitude image file
//* @param FPAR_DT					FPAR image file pointer
//* @param j						Line j of the the land use image
//* @param Lati_xrc					The number of columns for the land-use image
//* @param DEMM_xrc					Number of columns for 500m DEM
//* @param STRMM_DEM_GeoTransform	DEM image (500m) spatial coordinate information
//* @param pafLongH					The longitude pointer of a row of data in a latitude and longitude image (500m).
//* @param pafLatH					The latitude pointer of a row of data in a latitude and longitude image (500m).
//* @param pafLaiH					Pointer to a row of data in LAI image (500M)
//* @param pafLaimaxH				Pointer to a row of data in LAImax image (500M)
//* @param pafLandcoverH			Pointer to a row of data in landcover image (500M).
//* @param pafRmH					Pointer to a row of data in maintenance respiration image (500M).
//* @param pafTavH					Pointer to a row of data in average temperature image (500M)
//* @param pafGppH					Pointer to a row of data in GPP image (500M)
//* @param pafNppH					Pointer to a row of data in NPP image (500M)
//*/
///revised by MJ Wang
void comNPP(GDALDataset *LAI_DT, GDALDataset *LAI_MAX_Y,int j, int Lati_xrc, 
	DT_16U* pafLaiH, double* pafLaimaxH, DT_8U* pafLandcoverH, float* pafRmH, float* pafTavH,
			double* pafGppH, double* pafNppH,float LAInodata);

void comET(GDALDataset *Lati_DT, int j, int Lati_xrc, int GLDAS_xrc,  double GLDAS_GeoTransform[6],
				   float* pafLongH, float* pafLatH, float* pafTavH, float* pafTavGldasL,float* pafSWnetH, float* pafSWnetL, float* pafLWnetH, float* pafLWnetL, float* pafDTH, float* pafDTL,float* pafDEMH,
				   GDALDataset *FPAR_DT, GDALDataset *LAI_DT, GDALDataset *LANDCOVER_DT, DT_8U* pafFparH, DT_16U* pafLaiH, DT_8U* pafLandcoverH, double* pafEvaporation, double* pafTransporation, double* pafEvapotransporation, double* pafPEvapotransporation, float* paff2H, float* pafPar);////revised by MJ Wang 2019.4.24
//revised by MJ Wang 2019.4.25
void comCI(GDALDataset *Lati_DT, int j, int d,int Lati_xrc, int GLDAS_xrc, double GLDAS_GeoTransform[6], float* pafLatH, float* pafParH,float* CI);


//The final result is saved in TIF format and imported into the result line by line
//**
//* @brief							GDALDataset of the generated file
//* @return							No returns
//* @param poDataset				A file pointer to the NPP that is generated
//* @param savePath					The path where the final generated NPP file is saved
//* @param imgsizeX					The number of columns in the final generated NPP file
//* @param imgsizeY					The number of rows in the final generated NPP file
//* @param imgsizeBand				The number of bands in the final generated NPP file
//* @param dataType					The data type of the final generated NPP file
//* @param adfGeoTransform			The coordinate information of the final generated NPP file
//* @param projectionref			The projection information of the final generated NPP file
//*/
void save2TIFF(GDALDataset *poDataset, const char* savePath,int imgsizeX,int imgsizeY,int imgsizeBand,GDALDataType dataType,
	double* adfGeoTransform,const char* projectionref);

void	makeHDF(char* hdfchar,const char* gppsavePath,const char* nppsavePath,char* AcquisitionTime,char* ProductionTime,char* RawDataName,
             char* size,char* AlgorithmName,char* GridNum,char* QPProductName,char* NumBand,char* Projection,char* ProjectionStr,char* ProjectionPara);
void Clean_Dir(const char *dir,int depth,string filetype);
int ReadTif2Mem(const char* c_tif_file,const char* save_file,char *para);
void WriteTif2H5(const char* savepath,char*DATA1,int line,int sample); 