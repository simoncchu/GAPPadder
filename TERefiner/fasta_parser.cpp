#include"fasta_parser.h"
#include<fstream>
#include<iostream>
#include <sstream>
#include "public_func.h"
#include "public_parameters.h"


using namespace std;


FastaParser::FastaParser()
{
	this->vid_name.clear();
}

FastaParser::FastaParser(std::string path)
{
	this->path=path;
	this->vid_name.clear();
}

void FastaParser::setPath(std::string path)
{
	this->path=path; 
}

/*
Description: return the specified sub-sequence.
*/
string FastaParser::parseFasta(std::string chrom, long pos, long length)
{	
//cout<<chrom<<" "<<pos<<" "<<length<<endl;
	long chrom_start=0;//start position of chromosome, including all ascii chars, like '\0''\n'.
	int size_each_line=0;//number of symbols, like without  '\0\n'.
	int size_ascii=0;//number of all ascii chars, like with  '\0\n'.
	long chrom_len=0;
	parseFai(chrom, chrom_start, size_each_line, size_ascii, chrom_len);//get chrom start position from fasta.fai file.
//cout<<"chrom info:"<< chrom_start<<" "<<size_each_line<<" "<<size_ascii<<" "<<chrom_len<<endl;
	if(chrom_start==0)
	{
		return "";
	}

	if(length==-1) //get the whole length;
		length=chrom_len;

	ifstream fin;
	fin.open(path.c_str(), ifstream::in);//open fasta file 

//cout<<"Fasta is open"<<endl;

	long start_part_len=size_each_line-size_ascii;
	long lines=0;
	long left=pos%size_each_line;
	if(pos>size_each_line)
	{
		fin.seekg(chrom_start);//move to chrom_start position 
		string start_part;
		fin>>start_part;
		start_part_len=start_part.length();

		lines=(pos-start_part_len)/size_each_line;
		left=(pos-start_part_len)%size_each_line;

	}
	long newpos=(start_part_len+(size_ascii-size_each_line)) + lines*size_ascii + left;
	
//cout<<"File is seek to "<<newpos+chrom_start<<endl;

	fin.seekg(newpos+chrom_start);

	string sub_ref="";
	long cnt=0;
	while(true)
	{
		string temp;
		fin>>temp;
		cnt+=temp.length();
		if(cnt>=length)
		{
			sub_ref+=temp.substr(0,(size_each_line-(cnt-length)));
			break;
		}
		else
		{
			sub_ref+=temp;
		}
	}
	fin.close();

	return sub_ref;
}

void FastaParser::mapChromIDName()//map chrom ID with Name, by using XX.fasta.fai file. 
{
	ifstream fin;
	std::string fai_path=path+".fai";
	fin.open(fai_path.c_str(), ifstream::in);

	int id=0;
	ofstream fout;
	fout.open(CHROM_ID_NAME.c_str(),ofstream::out);
	std::string chrom_name;
	while(fin>>chrom_name)
	{
		char buffer[MAX_LEN_EACH_LINE_FAI];
		fin.getline(buffer,MAX_LEN_EACH_LINE_FAI);
		
		fout<<id<<" "<<chrom_name<<endl;
		id++;
	}
	fout.close();
	fin.close();
}

void FastaParser::loadChromIDName()
{
	ifstream fin;
	fin.open(CHROM_ID_NAME.c_str(), ifstream::in);
	int id;
	string name;
	vid_name.clear();
	while(fin>>id>>name)
	{
		this->vid_name.push_back(std::make_pair(id,name));
	}
	fin.close();
}

//------------private functions------------------------------------------------------
/*
Description: 
	Parse .fasta.fai file, and returen chromosome start position and number of chars each line.
*/
void FastaParser::parseFai(std::string chrom, long& start, int& size_aline, int& size_all_aline, long& chrom_len)
{
	ifstream fin;
	std::string fai_path=path+".fai";
	fin.open(fai_path.c_str(), ifstream::in);
	
	string chrom_name, slength, sstart_pos, ssize_symbol, ssize_ascii;
	while(!fin.eof())
	{
		string aline="";
		getline(fin,aline);
		std::stringstream ss(aline);
		std::getline(ss, chrom_name, '\t');
		if(chrom_name==chrom)
		{
			std::getline(ss, slength, '\t');
			chrom_len=PubFuncs::cvtStr2Int(slength);
			std::getline(ss, sstart_pos, '\t');
			start=PubFuncs::cvtStr2Int(sstart_pos);
			std::getline(ss, ssize_symbol, '\t');
			size_aline=PubFuncs::cvtStr2Int(ssize_symbol);
			std::getline(ss, ssize_ascii, '\t');
			size_all_aline=PubFuncs::cvtStr2Int(ssize_ascii);
			break;
		}
	}

	fin.close();
}