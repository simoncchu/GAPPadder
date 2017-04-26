#include<fstream>
#include"fai_parser.h"
#include<iostream>
#include"public_func.h"

using namespace std;

FaiParser::FaiParser()
{
	this->path="";
}
	
FaiParser::FaiParser(std::string path)
{
	this->path=path;
}

void FaiParser::setPath(std::string path)
{
	this->path=path;
}

/*
Description
	Parse .fasta.fai file, and returen chromosome start position and number of chars each line.
*/
void FaiParser::parseFai()
{
	ifstream fin;
	fin.open(path.c_str(), ifstream::in);

	string chrom_name, slength, sstart_pos, ssize_symbol, ssize_ascii;
	int cnt=0;
	vchroms.clear();
	
	while(!fin.eof())
	{
		string aline="";
		getline(fin,aline);
		ChromInfo ci;
		ci.id=cnt;
		std::stringstream ss(aline);

		std::getline(ss, chrom_name, '\t');
		ci.cname=chrom_name;
		std::getline(ss, slength, '\t');
		ci.length=PubFuncs::cvtStr2Int(slength);
		std::getline(ss, sstart_pos, '\t');
		ci.startpos=PubFuncs::cvtStr2Int(sstart_pos);
		std::getline(ss, ssize_symbol, '\t');
		ci.size_chars=PubFuncs::cvtStr2Int(ssize_symbol);
		std::getline(ss, ssize_ascii, '\t');
		ci.size_ascii=PubFuncs::cvtStr2Int(ssize_ascii);
		
		vchroms.push_back(ci);
		cnt++;
	}

	fin.close();
}

/*
Description:
	Parse .fasta.fai file, and returen chromosome length.
Input:
	chrom name
Output:
	length of chrom
*/
int FaiParser::getChromLen(std::string chrom)
{	
	int size=vchroms.size();
	for(int i=0;i<size;i++)
	{
		if(vchroms[i].cname==chrom)
		{
			return vchroms[i].length;
		}
	}
	return 0;
}