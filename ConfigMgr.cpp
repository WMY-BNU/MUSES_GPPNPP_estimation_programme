// ConfigMgr.cpp: implementation of the CConfigMgr class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "ConfigMgr.h"
#include "Markup.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CConfigMgr::CConfigMgr()
{

}

CConfigMgr::~CConfigMgr()
{

}

bool CConfigMgr::Parse(const char *pszFile)
{
	CMarkup xml;
	if (!xml.Load(pszFile))
	{
		return false;
	}

	if (xml.FindElem("parameters"))
	{
		xml.IntoElem();
		string strCode;
		while (xml.FindElem("param"))
		{
			string strName = xml.GetAttrib("name");
			string strData = xml.GetData();
			m_mapData[strName] = strData;
		}
		xml.OutOfElem();
	}
	return true;
}

string CConfigMgr::GetPathBy(const char *pszKey)
{
	map<string, string>::iterator iterFind = m_mapData.find(pszKey);
	if (iterFind != m_mapData.end())
	{
		return iterFind->second;
	}
	return "";
}