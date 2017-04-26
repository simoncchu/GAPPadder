//
//  ScaffoldUtils.cpp
//  
//
//  Created by Chong Chu on 3/24/15.
//
//

#include "ScaffoldUtils.h"
#include<fstream>

using namespace std;

//load scaffold information
void loadScaffoldInfo(const char* sfilename, MapScaffoldInfo& mapScaffoldInfo, int miniSupportPairCutOff)
{
	ifstream finScaffold;
	finScaffold.open(sfilename);

	//id_1 contig_1 length_1 direction_1 id_2 contig_2 length_2 direction_2 Num-of-supported-pairs Min-Distance Max-Distance Mean-Distance
	int id1, id2, len1, len2, numOfSupport;
	double minDist, maxDist, meanDist;
	string contigName1,contigName2, dirct1, dirct2;
	while(finScaffold>>id1>>contigName1>>len1>>dirct1>>id2>>contigName2>>len2>>dirct2>>numOfSupport>>minDist>>maxDist>>meanDist)
	{
		//if(minDist> 0) continue;//distance is larger than 0, means no overlap 
		if(numOfSupport < miniSupportPairCutOff)
		{
			continue;
		}

		if(dirct1=="-")
			contigName1+="_R";
		if(dirct2=="-")
			contigName2+="_R";

		//insert [contigName1,[contigName2,1]]
		if(mapScaffoldInfo.find(contigName1)==mapScaffoldInfo.end())
		{//not found
			MapContigInfo mci;
			mci[contigName2]=1;
			mapScaffoldInfo[contigName1]=mci;
		}
		else
		{//found
			mapScaffoldInfo[contigName1][contigName2]=1;
		}

		//insert [contigName2,[contigName1,1]]
		if(mapScaffoldInfo.find(contigName2)==mapScaffoldInfo.end())
		{//not found
			MapContigInfo mci;
			mci[contigName1]=1;
			mapScaffoldInfo[contigName2]=mci;
		}
		else
		{//found
			mapScaffoldInfo[contigName2][contigName1]=1;
		}
		
	}

	finScaffold.close();
}