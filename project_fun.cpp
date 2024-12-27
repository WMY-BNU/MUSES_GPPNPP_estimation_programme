#include "stdafx.h"
#include "IMGALG.h"
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <string.h>
#include "stdlib.h"
#include"ogr_spatialref.h"
#include "npp_fun.h"
  
using namespace std;

/** 
* Projection transforms
**/

bool Projection2ImageRow(double* adfGeoTransform,double dProjX,double dProjY,int&iCol,int&iRow)
{
	try
	{
		iCol=int((adfGeoTransform[0] -dProjX)/adfGeoTransform[1] );
		iRow=int((adfGeoTransform[3] -dProjY)/adfGeoTransform[5] );
		return true;
	}
		catch(...)
	{
		return false;
	}
}
bool ImageRowCol2Projection(double *adfGeoTransform, int iCol, int iRow, double &dProjX, double &dProjY)
{
	//adfGeoTransform[6]  The array adfGeoTransform holds some of the parameters in the affine transformation, which are explained below
	//adfGeoTransform[0]  X coordinates in the upper left corner 
	//adfGeoTransform[1]  East-west resolution
	//adfGeoTransform[2]  Rotation angle, 0 indicates the image "North up"
	//adfGeoTransform[3]  Y coordinates in the upper left corner 
	//adfGeoTransform[4]  Rotation angle, 0 indicates the image "North up"
	//adfGeoTransform[5]  North-south resolution

	try
	{
		dProjX = adfGeoTransform[0] + adfGeoTransform[1] * iCol + adfGeoTransform[2] * iRow;
		dProjY = adfGeoTransform[3] + adfGeoTransform[4] * iCol + adfGeoTransform[5] * iRow;
		return true;
	}
	catch(...)
	{
		return false;
	}
}
