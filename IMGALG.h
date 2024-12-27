/***************************************************************************
*
* Time: 2009-09-21
* Project: Remote sensing platform
* Purpose: Core library files
* Author:  Li Minlu
* Copyright (c) 2009, liminlu0314@gmail.com
* Describe:Provides definitions of commonly used data types, etc
*
****************************************************************************/
#pragma once
//#ifndef IMGALG_H
#define IMGALG_H

/**
* \file ImgCore.h
* @brief Core type definitions
*
* Export interface (in C) and define the core type
*/

/**
* Ignore warnings on the MS Widows platform lml 2010-10-19
* warning C4100: ¡°*¡±: Unreferenced parameters
* warning C4190: ¡°identifier1¡±There is a specified C link, but a UDT "identifier2" that is not compatible with C is returned
* warning C4251: The class "type" needs to be used by the client of the class "type2" using the DLL interface
* warning C4275: The non-DLL interface class key "identifier" is used as the base for the DLL interface class key "identifier".
* warning C4305: Truncated from "type1" to "type2".
* warning C4309: Truncation constants
* warning C4819: The file contains characters that cannot be represented on the current code page (936). Please save the file in Unicode format to prevent data loss
* warning C4996: Non-standard extensions are used: Enumerations are used in qualifiers
*/
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
#pragma warning(disable: 4100 4190 4251 4275 4305 4309 4819 4996 )
#endif

/**
* @brief Whether to use the Vld memory leak monitoring tool
*/
#if _USEVLD
#if _DEBUG	//Detect memory leaks in debug mode
#include "vld.h"
#endif
#endif

/**
* @brief SPECIFIES WHETHER TO USE THE LOG TOOL TO WRITE LOGS
*/
#if _USELOG
#define USE_LOG4CPP
#endif

#include <float.h>
#include <algorithm>
#include <deque>
#include <fstream>
#include <limits>
#include <map>
#include <stack>
#include <string>
#include <vector>


/**
* @brief Export symbol definitions
*/
#ifdef IMGALG_EXPORTS
#define IMGALG_API __declspec(dllexport)
#else
#define IMGALG_API __declspec(dllimport)
#endif

/**
* @brief Define NULL
*/
#ifndef NULL
#  define NULL  0
#endif

/**
* @brief Define FALSE
*/
#ifndef FALSE
#  define FALSE 0
#endif
#ifndef proLEN
#	define proLEN 200
#endif
/**
* @brief Define TRUE
*/
#ifndef TRUE
#  define TRUE  1
#endif

#ifndef MAX
/*! Find the maximum */
#  define MIN(a, b)      ((a<b) ? a : b)
/*! Find the minimum */
#  define MAX(a, b)      ((a>b) ? a : b)
#endif

/**
* @brief Define ABS and find the absolute value
*/
#ifndef ABS
#  define ABS(x)        ((x<0) ? (-1*(x)) : x)
#endif

/**
* @brief Define PI=3.141592653... and degree and radian conversion
*/
#ifndef M_PI
/*! Define Pi PI */
# define M_PI  3.1415926535897932384626433832795
/*! Radians to degrees */
# define DEG_PER_RAD      ((double)(180.0/M_PI))
/*! degrees to Radians */
# define RAD_PER_DEG      ((double)(M_PI/180.0))
#endif

/**
* @brief Define the square
*/
#ifndef M_SQUARE
# define M_SQUARE(x)  (x)*(x)
#endif
#ifndef NoData
# define NoData 1000
#endif
/**
* @brief Define the cube
*/
#ifndef M_CUBE
# define M_CUBE(x)  (x)*(x)*(x)
#endif

///*! Determine whether the floating-point number is NaN */
//inline bool isnan(const float& v)  { return _isnan(v) ? true : false; }
///*! Determine whether the double number is NaN */
//inline bool isnan(const double& v) { return _isnan(v) ? true : false; }
///*! Obtain the NaN value of double */
//inline double nan() { return numeric_limits<double>::quiet_NaN(); }

/**
* @brief Extremum of float type
*/
#ifndef FLT_EQUALS
/*! Whether floating-point numbers are equal */
#define FLT_EQUALS(x, y)  (fabs((double)x-y)<FLT_EPSILON)
/*! Whether floats are equal (specifies the comparison threshold) */
#define FLT_EQUALS_N(x, y, z)  (fabs((double)x-y)<z)
#endif

#ifndef FLT_ZERO
/*! Whether the floating-point number is 0 or not */
#define FLT_ZERO(x)  (fabs(x)<FLT_EPSILON)
#endif

/**
* @brief Release the array
*/
#define RELEASE(x)	if(x!=NULL) {delete []x; x = NULL;}

/**
* @brief Set the global region as the default region of the operating system
*/
#define SET_LOCAL	{ locale::global(locale("")); setlocale(LC_ALL,"Chinese-simplified"); }
/**
* @brief Reset the global locale
*/
#define REVERT_LOCAL	locale::global(locale("C"))


#ifndef EQUAL
#if defined(WIN32) || defined(WIN32CE)
/*! Compare whether the strings are equal */
#  define EQUALN(a, b, n)           (_strnicmp(a, b, n) == 0)
/*! Compare whether the strings are equal */
#  define EQUAL(a, b)              (_stricmp(a, b) == 0)
#else
/*! Compare whether the strings are equal */
#  define EQUALN(a, b, n)           (strncasecmp(a, b, n) == 0)
/*! Compare whether the strings are equal */
#  define EQUAL(a, b)              (strcasecmp(a, b) == 0)
#endif
#endif
//#endif
/*! byte */
typedef unsigned char	byte;
/*! 8U */
typedef unsigned char	DT_8U;
/*! 16U */
typedef unsigned short	DT_16U;
/*! 16S */
typedef short			DT_16S;
/*! 32U */
typedef unsigned int	DT_32U;
/*! 32S */
typedef int				DT_32S;
/*! 32F */
typedef float			DT_32F;
/*! 64F */
typedef double			DT_64F;

/*! Successfully executed */
const int RE_SUCCESS		= 0;
/*! The file does not exist */
const int RE_FILENOTEXIST	= 1;
/*! The file format is not supported */
const int RE_FILENOTSUPPORT	= 2;
/*! The image data type is incorrect */
const int RE_FILETYPEERROR	= 3;
/*! Failed to create image */
const int RE_CREATEFAILED	= 4;
/*! The input parameter is incorrect */
const int RE_PARAMERROR		= 5;
/*! Other errors */
const int RE_FAILED			= 6;
/*! There is no public area for the image */
const int RE_NOSAMEEXTENT	= 7;
/*! The user cancels the action */
const int RE_USERCANCEL		= 8;
/*! The file is already in use */
const int RE_FILEISUESED	= 9;
/*! Unsupported pixel depth */
const int RE_DEPTHNOTSUPPORT	= 10;
/*! The number of bands does not meet the requirements */
const int RE_BANDCOUNTERROR		= 11;
/*! The file does not have a projection */
const int RE_NOPROJECTION		= 12;
/*! The projections are inconsistent */
const int RE_PROJECTIONDIFF		= 13;
/***************************************************************************
*
* Time: 2013-05-29
* Project: NPP estimates
* Purpose: Core library files
* Author:  Hao Lvyuan
* Copyright (c) 2013, Rwatermoon@gmail.com
* Describe:Related macro definitions
*
****************************************************************************/
/*! Light use efficiency kg/MJ*/
const float LUE[]={0,1.008,1.159,1.103,
	1.044,1.116,0.888,0.774,0.800,0.768,
	0.680,0.680,0.680,0.680,0.680,0,0};
/*! Temperature correction parameters */
const float AIR_COR =0.0065;
const float Topt=20;
//Kelvin temperature
#ifndef K_value
#define K_value 273.16
#endif