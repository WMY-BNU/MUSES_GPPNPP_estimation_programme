#include "stdafx.h"
#include "IMGALG.h"
#include <gdal_priv.h>
#include <iostream>
#include <fstream>
#include <cpl_conv.h>
#include <string.h>
#include "stdlib.h"
#include "npp_fun.h"


using namespace std;

/** 
* A recording function that outputs temporary text, which mainly records the running status of the program
**/

void makeLogo(const char* logopath,const char* content){
	ofstream in;
	time_t timep;
	time (&timep);
	char *pszCurrTime = (char*)malloc(sizeof(char)*20);
	memset(pszCurrTime, 0, sizeof(char)*20);
	strftime(pszCurrTime, 20 , "%Y/%m/%d %H:%M:%S", localtime(&timep));
	in.open(logopath,ios::app);
	in<<pszCurrTime<<"\t"<<content<<"\n";
	in.close();

}