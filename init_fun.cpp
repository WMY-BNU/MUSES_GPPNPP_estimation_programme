#include "stdafx.h"
#include "IMGALG.h"
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <string.h>
#include "stdlib.h"
#include "npp_fun.h"

using namespace std;

/** 
* Initialize the pointer function
**/

void initPoint(double* paf,int count,double value)
{
	for(int i=0;i<count;i++)
		paf[i]=value;
}
void initPoint(DT_16U* paf,int count,double value)
{
	for(int i=0;i<count;i++)
		paf[i]=value;
}
void initPoint(DT_8U* paf,int count,DT_8U value)
{
	for(int i=0;i<count;i++)
		paf[i]=value;
}
void initPoint(DT_8U* paf,int count)
{
	for(int i=0;i<count;i++)
		paf[i]=(DT_8U)0;
}
void initPoint(int* paf,int count,int value)
{
	for(int i=0;i<count;i++)
		paf[i]=value;
}
void initPoint(float* paf,int count)
{
	for(int i=0;i<count;i++)
		paf[i]=0;
}
void initPoint(float* paf,int count,float value)
{
	for(int i=0;i<count;i++)
		paf[i]=value;
}
void initPoint(DT_32U* paf,int count,double value)
{
	for(int i=0;i<count;i++)
		paf[i]=value;
}
