#include"Alignment.h"
#include"public_parameters.h"
#include"public_func.h"
#include<string>
#include<iostream>

using namespace std;


Alignment::Alignment(BamAlignmentRecord*& bar)
{
	this->bar=bar;
	cigar="";
}

void Alignment::setBar(BamAlignmentRecord*& bar)
{
	this->bar=bar;
	cigar="";
}

/*
Description:
	return alignment type according to flag. 
*/
int Alignment::getReadType()
{
	int flag=this->bar->flag;
	if((flag&0x4)==0)
	{//read mapped
		if((flag&0x800)==0)
		{//not supplementary map
			//check whether soft-clip
			if(isClipped()==false)
			{//fully mapped 
				if(isFullyMapped()==true)
					return READ_TYPE_FULLMAP;
				else 
					return READ_TYPE_OTHER;
			}
			else
			{//clipped
				return READ_TYPE_CLIP;
			}
		}
		else
		{//supplementary map, like hard-clip
			//check whether hard-clip;
			if(this->isHardClipped() == true)
				return READ_TYPE_CLIP;

			return READ_TYPE_OTHER;
		}
	}
	else
	{//read unmapped 
		//check whether mate is mapped
		if((flag&0x8)==0)
		{//mate is mapped 
			return READ_PAIR_MAP_TYPE_01;
		}
		else
		{//mate is unmapped
			return READ_PAIR_MAP_TYPE_00;
		}
	}
}

int Alignment::getPEReadType()
{
	int flag=this->bar->flag;
	if((flag&0x4)==0)
	{//read mapped
		if((flag&0x8)==0)
		{//mate is mapped 
			return READ_PAIR_MAP_TYPE_11;
		}
		else
		{//mate is unmapped
			return READ_PAIR_MAP_TYPE_10;
		}
	}
	else
	{//read unmapped
		if((flag&0x8)==0)
		{//mate is mapped 
			return READ_PAIR_MAP_TYPE_01;
		}
		else
		{//mate is unmapped 
			return READ_PAIR_MAP_TYPE_00;
		}
	}
}

bool Alignment::isFullyMapped()
{
	int size=bar->cigar.size();
	if(size==1 && bar->cigar[0].second==READ_LENGTH && bar->cigar[0].first=="M")
	{
		return true;
	}
	else
	{
		int cnt=0;
		for(int i=0;i<size;i++)
		{
			if(bar->cigar[i].first=="M")
			{
				cnt+=bar->cigar[i].second;
			}
		}
		double percent=(double)cnt/(double)READ_LENGTH; 
		if( percent > READ_FULL_MAPPED_CUTOFF)
			return true;
		else
			return false;
	}
	return false;
}



/*
Description:
	check whether alignment is clipped. 
Output:
	0, fully mapped.
*/
bool Alignment::isClipped()
{
	int size=bar->cigar.size();
	int cnt=0;
	if(size==1 && bar->cigar[0].second==READ_LENGTH && bar->cigar[0].first=="M")
	{//perfect fully mapped 
		return false;
	}
	else 
	{//clipped or other
		if(size>1 && (bar->cigar[0].first=="S" || bar->cigar[0].first=="H"))
			return true;
		else if(size>1 && (bar->cigar[size-1].first=="S" || bar->cigar[size-1].first=="H"))
			return true;
	}
}


bool Alignment::isHardClipped()
{
	int size=bar->cigar.size();
	if(size>1 && (bar->cigar[0].first=="H" || bar->cigar[size-1].first=="H"))
	{
		return true;
	}
	return false;
}

/*
Output:
	1, left soft-clip
	2, right soft-clip
	3, both side soft-clip 
	4, left hard-clip
	5, right hard-clip
	6, both hard-clip.
*/
int Alignment::getClipType(int& pos1, int& len1, int& pos2, int& len2)
{
	pos1=-1;pos2=-1;len1=-1;len2=-1;

	int size=bar->cigar.size();

	if(size>1 && bar->cigar[0].first=="S" && bar->cigar[size-1].first=="S")
	{//both soft-clipped 
		pos1=0;
		len1=bar->cigar[0].second;
		pos2=READ_LENGTH-bar->cigar[size-1].second;
		len2=bar->cigar[size-1].second;
		return READ_TYPE_BOTH_SOFTCLIP;
	}
	else if(size>1 && bar->cigar[0].first=="S")
	{//left soft-clip
		pos1=0;
		len1=bar->cigar[0].second;
		return READ_TYPE_LEFT_SOFTCLIP;
	}
	else if(size>1 && bar->cigar[size-1].first=="S")
	{//right soft-clip
		pos2=READ_LENGTH - bar->cigar[size-1].second - 1;
		len2=bar->cigar[size-1].second;
		return READ_TYPE_RIGHT_SOFTCLIP;
	}
	else if(size>1 && bar->cigar[0].first=="H" && bar->cigar[size-1].first=="H")
	{
		pos1=0;
		len1=bar->cigar[0].second;
		pos2=READ_LENGTH - bar->cigar[size-1].second;
		len2=bar->cigar[size-1].second;
		return READ_TYPE_BOTH_HARDCLIP;
	}
	else if(size>1 && bar->cigar[0].first=="H")
	{
		pos1=0;
		len1=bar->cigar[0].second;
		return READ_TYPE_LEFT_HARDCLIP;
	}
	else if(size>1 && bar->cigar[size-1].first=="H")
	{
		pos2=READ_LENGTH - bar->cigar[size-1].second - 1;
		len2=bar->cigar[size-1].second;
		return READ_TYPE_RIGHT_HARDCLIP;
	}
	else
	{
		return READ_TYPE_OTHER;
	}
}

std::string Alignment::getCigar()
{
	this->cigar="";
	int size=bar->cigar.size();
	for(int i=0;i<size;i++)
	{
		cigar+=PubFuncs::cvtInt2Str(bar->cigar[i].second);
		cigar+=bar->cigar[i].first;
	}
	return this->cigar;
}

//convert string to cigar vector
void Alignment::str2Cigar(std::string cigar, std::vector<std::pair<std::string, int> >& vcigar)
{
	int len=cigar.length();
	string temp="";
	for(int i=0;i<len;i++)
	{
		if(cigar[i]>='0' && cigar[i]<='9')
		{
			temp+=cigar[i];
		}
		else
		{
			int value=PubFuncs::cvtStr2Int(temp);
			temp="";
			string opt=""; opt+=cigar[i];
			vcigar.push_back(std::make_pair(opt,value));
		}
	}
	
}

bool Alignment::isDuplicate()//whether optimal duplicate 
{
	
	if((this->bar->flag & 0x400) == 0)
	{//not duplicate
		return false;
	}
	else 
		return true;
}

bool Alignment::isPrimaryAlign()//whether is primary alignment
{
	if((this->bar->flag & 0x100) == 0)
	{//primary alignment
		return true;
	}
	else
		return false;
}

bool Alignment::passQualityCK()//whether fails quality check 
{
	if((this->bar->flag & 0x200) == 0)
	{//pass quality control
		return true;
	}
	else
		return false;
}

bool Alignment::isMateMapped()
{
	if((this->bar->flag & 0x8) == 0)
	{//mate is mapped 
		return true;
	}
	else
	{
		return false;
	}
}

bool Alignment::isFirstInPair()//for pair-end reads, it is the first in pair.
{
	if((this->bar->flag & 0x40) == 0)
	{//not first in pair
		return false;
	}
	else
		return true;
}

bool Alignment::isSecondInPair()
{
	if((this->bar->flag & 0x80) == 0)
	{//not first in pair
		return false;
	}
	else
		return true;
}

bool Alignment::isSecondaryAlignment()//whether it is secondary alignment
{
	if((this->bar->flag & 0x100) == 0)
	{//is not secondary alignment
		return false;
	}
	else
		return true;
}

bool Alignment::isReverseStrand()
{
	if((this->bar->flag & 0x10) == 0)
	{
		return false;
	}
	else
		return true;
}

bool Alignment::isMateReverseStrand()
{
	if((this->bar->flag & 0x20) == 0)
	{
		return false;
	}
	else
		return true;
}

//////added 11/17/14///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//For some sam/bam files, there are varied read length. So we need to give the specific length and then check read alignment type 
/*
Description:
	return alignment type according to flag.
*/
int Alignment::getReadType(int rlength)
{
	int flag=this->bar->flag;
	cout<<"flag is "<<flag<<endl;///////////////////////////////////////////////////////////////////////////////////
	if((flag&0x4)==0)
	{//read mapped
		if((flag&0x800)==0)
		{//not supplementary map
			//check whether soft-clip
			if(isClipped(rlength)==false)
			{//fully mapped 
				if(isFullyMapped(rlength)==true)
					return READ_TYPE_FULLMAP;
				else 
					return READ_TYPE_OTHER;
			}
			else
			{//clipped
				return READ_TYPE_CLIP;
			}
		}
		else
		{//supplementary map, like hard-clip 
			//check whether hard-clip;
			if(this->isHardClipped() == true)
				return READ_TYPE_CLIP;

			return READ_TYPE_OTHER;
		}
	}
	else
	{//read unmapped 
		//check whether mate is mapped
		if((flag&0x8)==0)
		{//mate is mapped 
			return READ_PAIR_MAP_TYPE_01;
		}
		else
		{//mate is unmapped
			return READ_PAIR_MAP_TYPE_00;
		}
	}
}


bool Alignment::isFullyMapped(int rlength)
{
	int size=bar->cigar.size();
	if(size==1 && bar->cigar[0].second<=rlength && bar->cigar[0].first=="M")
	{
		return true;
	}
	else
	{
		int cnt=0;
		int total_len=0;
		for(int i=0;i<size;i++)
		{
			if(bar->cigar[i].first=="M")
			{
				cnt += bar->cigar[i].second;
				total_len += bar->cigar[i].second;
			}
			if(bar->cigar[i].first=="S" || bar->cigar[i].first=="H" || bar->cigar[i].first=="I")
				total_len += bar->cigar[i].second;
		}
		double percent=(double)cnt/(double)total_len; 

		if( percent > READ_FULL_MAPPED_CUTOFF)
			return true;
		else
			return false;
	}
	return false;
}

bool Alignment::isPerfectMapped(int rlength)//perfect mapped
{
	int size=bar->cigar.size();
	if((size==1) && (bar->cigar[0].second == rlength) && (bar->cigar[0].first=="M"))
	{
		return true;
	}
	else
		return false;
}

bool Alignment::isClipped(int rlength)
{
	int size=bar->cigar.size();
	int cnt=0;
	if(size==1 && bar->cigar[0].second==rlength && bar->cigar[0].first=="M")
	{//fully mapped 
		return false;
	}
	else 
	{//clipped or other
		if(size>1 && (bar->cigar[0].first=="S" || bar->cigar[0].first=="H"))
			return true;
		else if(size>1 && (bar->cigar[size-1].first=="S" || bar->cigar[size-1].first=="H"))
			return true;
	}
}
