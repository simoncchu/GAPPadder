//
//  ScaffoldUtils.h
//  
//
//  Created by Chong Chu on 03/24/15.
//
//  Process Scaffold information 
//
#ifndef ___ScaffoldUtils_h__
#define ___ScaffoldUtils_h__

#include <map>
#include <string>
using namespace std;

typedef string ScaffoldType;
typedef map<ScaffoldType,int> MapContigInfo;
typedef map<ScaffoldType, MapContigInfo> MapScaffoldInfo;

//load scaffold information 
void loadScaffoldInfo(const char* sfilename, MapScaffoldInfo& mapScaffoldInfo, int miniSupportPairCutOff);

#endif