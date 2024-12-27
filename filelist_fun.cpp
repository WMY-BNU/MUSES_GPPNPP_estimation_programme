#include "stdafx.h"
#include "IMGALG.h"
#include <gdal_priv.h>
#include <iostream>
#include <fstream>
#include<algorithm>
#include <cpl_conv.h>
#include <string.h>
#include "stdlib.h"
#include "npp_fun.h"
#include"ogr_spatialref.h"
#include <windows.h>  
using namespace std;

/** 
* Get all file names of a certain file type in a specified folder under Windows system,
* and save in the vector of the string
*/ 


void get_filelist(const char *foldname,vector<string> &filelist,string filetype)
{
	//string filetype = ".nc";
    HANDLE hFind;
    WIN32_FIND_DATA fileData;
    string line;
    char fn[MAX_PATH];
    char tmpfn[MAX_PATH];
    strcpy(fn,foldname);

    //Need to process the string of the folder name
    if(fn[strlen(fn) -1] != '\\' )
    {
        strcat(fn, "\\");
    }

    //Note the order, at this point the fn has joined"\\"
    strcpy(tmpfn,fn);
    //If you don't add *, you'll get an error!
    strcat(fn, "*");

    hFind = FindFirstFile(fn, &fileData);
    FindNextFile(hFind, &fileData);
    while(FindNextFile(hFind, &fileData)){
        //If scaned a folder
        if (fileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
        {
            if(fileData.cFileName[0]!='.')
            {
                char szFile[MAX_PATH];
                strcpy(szFile,tmpfn);
                strcat(szFile,fileData.cFileName);
                get_filelist(szFile,filelist,filetype);
            }
        }
        //If scaned a file
        else
        {
            line = (string)tmpfn;
            line+= fileData.cFileName;
            if (line.find(filetype,0)!=string::npos)
            {
                filelist.push_back(line);
            }
            else continue;
        }
	}
	sort(filelist.begin(),filelist.end());
}

//Delete all temporary files in the folder
void Clean_Dir(const char *dir,int depth,string filetype)
{

}