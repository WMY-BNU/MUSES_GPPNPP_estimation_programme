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
* Interpolation and resampling of data
**/

float Distance(float x1,float x2,float y1,float y2)
{
	return M_SQUARE(x1-x2)+M_SQUARE(y1-y2);
}

float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	float* data,int yrc,int flag,float pix_size)
{
	//Handle cases where the interpolation point falls on the boundary
	if(flag==1)
		return data[x_low_index+y_low_index*yrc];
	float* p_max=(float*)malloc(sizeof(float)*4);
	p_max[0]=Distance(x_low,x,y_low,y);
	p_max[1]=Distance(x_low+pix_size,x,y_low,y);
	p_max[2]=Distance(x_low,x,y_low+pix_size,y);
	p_max[3]=Distance(x_low+pix_size,x,y_low+pix_size,y);
	float min=p_max[0];
	int min_index=0;
	for(int i=0;i<4;i++)
	{
		if (p_max[i] < min)
		{
			min=p_max[i];
			min_index=i;
		}
	}
	switch(min_index)
	{
	case 0:
		{return data[x_low_index+y_low_index*yrc];
		break;}
	case 1:
		{return data[x_low_index+1+y_low_index*yrc];
		break;}
	case 2:
		{return data[x_low_index+(y_low_index+1)*yrc];
		break;}
	case 3:
		{return data[x_low_index+1+(y_low_index+1)*yrc];
		break;}
	default:
		return 0;
	}
	delete(p_max);
	return 0;
}

float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_16U* data,int yrc,int flag,float pix_size)
{
	//Handle cases where the interpolation point falls on the boundary
	if(flag==1)
		return data[x_low_index+y_low_index*yrc];
	int* p_max=(int*)malloc(sizeof(int)*4);
	p_max[0]=Distance(x_low,x,y_low,y);
	p_max[1]=Distance(x_low+pix_size,x,y_low,y);
	p_max[2]=Distance(x_low,x,y_low+pix_size,y);
	p_max[3]=Distance(x_low+pix_size,x,y_low+pix_size,y);
	int min=p_max[0];
	int min_index=0;
	for(int i=0;i<4;i++)
	{
		if (p_max[i] < min)
		{
			min=p_max[i];
			min_index=i;
		}
	}
	switch(min_index)
	{
	case 0:
		{return data[x_low_index+y_low_index*yrc];
		break;}
	case 1:
		{return data[x_low_index+1+y_low_index*yrc];
		break;}
	case 2:
		{return data[x_low_index+(y_low_index+1)*yrc];
		break;}
	case 3:
		{return data[x_low_index+1+(y_low_index+1)*yrc];
		break;}
	default:
		return 0;
	}
	delete(p_max);
	return 0;
}

float Nearest(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_8U* data,int yrc,int flag,float pix_size)
{
	//Handle cases where the interpolation point falls on the boundary
	if(flag==1)
		return data[x_low_index+y_low_index*yrc];
	int* p_max=(int*)malloc(sizeof(int)*4);
	p_max[0]=Distance(x_low,x,y_low,y);
	p_max[1]=Distance(x_low+pix_size,x,y_low,y);
	p_max[2]=Distance(x_low,x,y_low+pix_size,y);
	p_max[3]=Distance(x_low+pix_size,x,y_low+pix_size,y);
	int min=p_max[0];
	int min_index=0;
	for(int i=0;i<4;i++)
	{
		if (p_max[i] < min)
		{
			min=p_max[i];
			min_index=i;
		}
	}
	switch(min_index)
	{
	case 0:
		{return data[x_low_index+y_low_index*yrc];
		break;}
	case 1:
		{return data[x_low_index+1+y_low_index*yrc];
		break;}
	case 2:
		{return data[x_low_index+(y_low_index+1)*yrc];
		break;}
	case 3:
		{return data[x_low_index+1+(y_low_index+1)*yrc];
		break;}
	default:
		return 0;
	}
	delete(p_max);
	return 0;
}
float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	float* data,int yrc,int flag,float pix_size,float nodata)
{
	//Handle cases where the interpolation point falls on the boundary
	if(flag==1)
		return data[x_low_index+y_low_index*yrc];
	if(data[x_low_index+y_low_index*yrc]==nodata||data[x_low_index+1+y_low_index*yrc==nodata]||
		data[x_low_index+(y_low_index+1)*yrc]==nodata||data[x_low_index+1+(y_low_index+1)*yrc]==nodata)
		return nodata;
	return data[x_low_index+y_low_index*yrc]*(x_low+pix_size-x)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+1+y_low_index*yrc]*(x-x_low)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+(y_low_index+1)*yrc]*(x_low+pix_size-x)*(y_low-y)/(pix_size*pix_size)+
		data[x_low_index+1+(y_low_index+1)*yrc]*(x-x_low)*(y_low-y)/(pix_size*pix_size);
}
float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	float* data,int yrc,int flag,float pix_size)
{
	//Handle cases where the interpolation point falls on the boundary
	if(flag==1)
		return data[x_low_index+y_low_index*yrc];
	return data[x_low_index+y_low_index*yrc]*(x_low+pix_size-x)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+1+y_low_index*yrc]*(x-x_low)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+(y_low_index+1)*yrc]*(x_low+pix_size-x)*(y_low-y)/(pix_size*pix_size)+
		data[x_low_index+1+(y_low_index+1)*yrc]*(x-x_low)*(y_low-y)/(pix_size*pix_size);
}

float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_16U* data,int yrc,int flag,float pix_size)
{
	//Handle cases where the interpolation point falls on the boundary
	if(flag==1)
		return data[x_low_index+y_low_index*yrc];
	return data[x_low_index+y_low_index*yrc]*(x_low+pix_size-x)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+1+y_low_index*yrc]*(x-x_low)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+(y_low_index+1)*yrc]*(x_low+pix_size-x)*(y_low-y)/(pix_size*pix_size)+
		data[x_low_index+1+(y_low_index+1)*yrc]*(x-x_low)*(y_low-y)/(pix_size*pix_size);
}

float MyBilinear(int x_low_index,int y_low_index,float x_low,float y_low,float x,float y,
	DT_8U* data,int yrc,int flag,float pix_size)
{
	//Handle cases where the interpolation point falls on the boundary
	if(flag==1)
		return data[x_low_index+y_low_index*yrc];
	return data[x_low_index+y_low_index*yrc]*(x_low+pix_size-x)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+1+y_low_index*yrc]*(x-x_low)*(y-y_low+pix_size)/(pix_size*pix_size)+
		data[x_low_index+(y_low_index+1)*yrc]*(x_low+pix_size-x)*(y_low-y)/(pix_size*pix_size)+
		data[x_low_index+1+(y_low_index+1)*yrc]*(x-x_low)*(y_low-y)/(pix_size*pix_size);
}


void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,float* pafData,int xrc,float nodata)
{

		//Interpolate a row of data
		for(int i=0;i<Lxrc;i++)
		{
			float long_dis=pafLong[i]-Lup_long;
			float lat_dis=Lup_lat-pafLat[i];
			//Obtain an index of the upper-left coordinates of the area located in the 0.25бу image
			int x_Low=int(long_dis/pix_size);
			int y_Low=int(lat_dis/pix_size);

			int flag=0;
			if(long_dis/pix_size==int(long_dis/pix_size)&&lat_dis/pix_size==int(lat_dis/pix_size))
				flag=1;
			if(rsMath==1)
			{
			//Interpolation is performed using the bilinear method
			pafnew[i]=MyBilinear(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
				pafLong[i],pafLat[i],pafData,xrc,flag,pix_size,nodata);
			}
			else
			{
			//Interpolation is performed using the nearest neighbor method
			pafnew[i]=Nearest(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
				pafLong[i],pafLat[i],pafData,xrc,flag,pix_size);
			}
		}
}
void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,float* pafData,int xrc)
{

		//Interpolate a row of data
		for(int i=0;i<Lxrc;i++)
		{
			float long_dis=pafLong[i]-Lup_long;
			float lat_dis=Lup_lat-pafLat[i];
			//Obtain an index of the upper-left coordinates of the area located in the 0.25бу image
			int x_Low=int(long_dis/pix_size);
			int y_Low=int(lat_dis/pix_size);

			//cout<<x_Low<<endl;
			//cout<<y_Low<<endl;

			int flag=0;
			if(long_dis/pix_size==int(long_dis/pix_size)&&lat_dis/pix_size==int(lat_dis/pix_size))
				flag=1;
			if(rsMath==1)
			{
			//Interpolation is performed using the bilinear method
			pafnew[i]=MyBilinear(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
				pafLong[i],pafLat[i],pafData,xrc,flag,pix_size);
			}
			else
			{
			//Interpolation is performed using the nearest neighbor method
			pafnew[i]=Nearest(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
				pafLong[i],pafLat[i],pafData,xrc,flag,pix_size);
			}
		}
}

void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,DT_16U* pafData,int xrc)
{

		//Interpolate a row of data
		for(int i=0;i<Lxrc;i++)
		{
			float long_dis=pafLong[i]-Lup_long;
			float lat_dis=Lup_lat-pafLat[i];
			//Obtain an index of the upper-left coordinates of the area located in the 0.25бу image
			int x_Low=int(long_dis/pix_size);
			int y_Low=int(lat_dis/pix_size);

			int flag=0;
			if(long_dis/pix_size==int(long_dis/pix_size)&&lat_dis/pix_size==int(lat_dis/pix_size))
				flag=1;
			if(rsMath==1)
			{
				//Interpolation is performed using the bilinear method
				pafnew[i]=MyBilinear(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
					pafLong[i],pafLat[i],pafData,xrc,flag,pix_size);
			}
			else
			{
				//Interpolation is performed using the nearest neighbor method
				pafnew[i]=Nearest(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
					pafLong[i],pafLat[i],pafData,xrc,flag,pix_size);
			}
		}
}

void resampleL2H(float* pafLong,float* pafLat,float pix_size,float Lup_long,float Lup_lat,float* pafnew,
	int rsMath,int Lxrc,DT_8U* pafData,int xrc)
{

		//Interpolate a row of data
		for(int i=0;i<Lxrc;i++)
		{
			float long_dis=pafLong[i]-Lup_long;
			float lat_dis=Lup_lat-pafLat[i];
			//Obtain an index of the upper-left coordinates of the area located in the 0.25бу image
			int x_Low=int(long_dis/pix_size);
			int y_Low=int(lat_dis/pix_size);

			int flag=0;
			if(long_dis/pix_size==int(long_dis/pix_size)&&lat_dis/pix_size==int(lat_dis/pix_size))
				flag=1;
			if(rsMath==1)
			{
				//Interpolation is performed using the bilinear method
				pafnew[i]=MyBilinear(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
					pafLong[i],pafLat[i],pafData,xrc,flag,pix_size);
			}
			else
			{
				//Interpolation is performed using the nearest neighbor method
				pafnew[i]=Nearest(x_Low,y_Low,Lup_long+x_Low*pix_size,Lup_lat-y_Low*pix_size,
					pafLong[i],pafLat[i],pafData,xrc,flag,pix_size);
			}
		}
}



