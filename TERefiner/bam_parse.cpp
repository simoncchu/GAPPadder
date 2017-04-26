#include"bam_parse.h"

#include<iostream>
#include <algorithm>
#include<cassert>

using namespace std;

BamParse::BamParse(string filename)
{
	this->filename=filename;
	this->filename_index=filename+".bai";
	min_read_length=500;
	max_read_length=0;
	insert_size=0;
	std_variation=0.0;
}

BamParse::BamParse(vector<string>*& listname)
{
	this->listfname=listname;
	min_read_length=500;
	max_read_length=0;
	insert_size=0;
	std_variation=0.0;
}

BamParse::~BamParse()
{
	int size=0;
	size=bam_head_records.size();
	for(int i=0;i<size;i++)
	{
		if(bam_head_records[i]!=NULL)
		{
			delete bam_head_records[i];
			bam_head_records[i]=NULL;
		}
	}
	vector<BamHeaderRecord*>().swap(bam_head_records);

	size=bam_aln_records.size();
	for(int i=0;i<size;i++)
	{
		if(bam_aln_records[i]!=NULL)
		{
			delete bam_aln_records[i];
			bam_aln_records[i]=NULL;
		}
	}
	vector<BamAlignmentRecord*>().swap(bam_aln_records);
}

void BamParse::parseCigar(vector<CigarOp>& cgdata, BamAlignmentRecord*& bar)//parse cigar
{
	pair<string,int> ptemp;
	for(int i=0;i<cgdata.size();i++)
	{
		ptemp.first=cgdata[i].Type;
		
		ptemp.second=cgdata[i].Length;

		bar->cigar.push_back(ptemp);
	}
}

void BamParse::parseHeader()
{
	
}

bool BamParse::loadIndex()
{
	bool bsignal = reader.OpenIndex(this->filename_index);
	if(bsignal==false) 
	{
		cerr << "Index file not found, now create it!!!"<<endl;
		bool bci=reader.CreateIndex();
		if(bci==false)
		{
			cerr << "Index file cannot be created!!!"<<endl;
			return false;
		}
		else
		{
			bool blct=reader.LocateIndex();
			if(blct==false)
			{
				cerr << "Index file cannot be located!!!"<<endl;
				return false;
			}
		}

	}
	return true;
}

/*
Description:
	Load all the bam into memory at one time if no parameters set, otherwise load the needed part of the bam.
	Save the parsed info into vector. 
*/
bool BamParse::parseAlignment(int chrom1, int chrom1_begin, int chrom2, int chrom2_end)
{
	//BamReader reader;
    if ( !reader.Open(filename) ) {
        cerr << "Bamtools ERROR: could not open input BAM file: " << filename << endl;
        return false;
    }
		
	//check whether need to set a region.
	if(chrom1>-1 && chrom1_begin>-1 && chrom2>-1 && chrom2_end>-1)
	{
		//this->loadIndex(reader);
		/*if(reader.HasIndex()==false)  //////////////////////////why return false, even have loaded the index???????????
		{
			cerr << "No index loaded!!"<<endl;
			return false;
		}*/

		BamRegion br(chrom1,chrom1_begin,chrom2,chrom2_end);
		bool is_set=reader.SetRegion(br);
		//if(is_set==false)
		//{
		//	cerr << "Cannot set the region!!"<<endl;
		//	return false;//cannot set the region.  
		//}
	}

	//process input data 
    BamAlignment al;
	while ( reader.GetNextAlignment(al) )
	{
		if(al.Position<0) continue;

		BamAlignmentRecord* bar=new BamAlignmentRecord();
		setAlignmentRecord(al,bar);
		setReadGroup(al,bar);//set read group information 
		bam_aln_records.push_back(bar);
	}

	reader.Close();
	return true;
}

bool BamParse::openReader()
{
	//BamReader reader;
    if ( !reader.Open(filename) ) {
        cerr << "Bamtools ERROR: could not open input BAM file: " << filename << endl;
        return false;
    }
}

bool BamParse::dumpAlignments(int chrom1, int chrom1_begin, int chrom2, int chrom2_end)
{
	bool rtn = true;
	//check whether need to set a region. 
	if(chrom1>-1 && chrom1_begin>-1 && chrom2>-1 && chrom2_end>-1)
	{
		//this->loadIndex(reader); 
		/*if(reader.HasIndex()==false)  //////////////////////////why return false, even have loaded the index???????????
		{
			cerr << "No index loaded!!"<<endl;
			return false;
		}*/

		BamRegion br(chrom1,chrom1_begin,chrom2,chrom2_end);
		bool is_set=reader.SetRegion(br);
		//if(is_set==false)
		//{
		//	cerr << "Cannot set the region!!"<<endl;
		//	return false;//cannot set the region.  
		//}
	}

	//process input data
    BamAlignment al;
	while ( reader.GetNextAlignment(al) )
	{
		if(al.Position<0) continue;

		BamAlignmentRecord* bar=new BamAlignmentRecord();
		setAlignmentRecord(al,bar);
		setReadGroup(al,bar);//set read group information 
		bam_aln_records.push_back(bar);
	}
	return rtn;
}

void BamParse::closeReader()
{
	reader.Close();
}


/*
Description:
	parse bam chrom by chrom.
Input:
	int chromID, start from 0;
	int length, chrom length.
*/
bool BamParse::parseAlignmentByChrom(int chromID,int length)
{
	bool bsignal=parseAlignment(chromID,0,chromID,length);//read the whole chromosome into memory 
	if(bsignal==false)
		std::cout<<"Cannot parse bam file!! "<<endl;//////////////////////////////////////////////////////////
	return bsignal;
}

void BamParse::setAlignmentRecord(BamAlignment& al,BamAlignmentRecord*& bar)
{
	// mandatory fields 
	//parse QNAME 
	bar->qName=al.Name;
	//parse FLAG
	bar->flag=al.AlignmentFlag;
	//parse REF SEQ ID 
	bar->rID=al.RefID;// start from 0, for human_g1k_v37, 0 represent chr1, and 23 represent chrx, and so on.
	//parse 1-based mappping position, POS 
	bar->pos=al.Position;
	//parse MAPQ
	bar->mapQ=al.MapQuality;
	//parse bin  
	bar->bin=al.Bin;
	//parse CIGAR
	parseCigar(al.CigarData, bar);
	//parse RNEXT ID
	bar->rNextID=al.MateRefID;
	//parse PNEXT, 1-based
	bar->pNext=al.MatePosition;
	//parse TLEN
	bar->tLen=al.Length;
	//parse SEQ
	//bar->seq="*";
	bar->seq=al.QueryBases;
	//parse QUAL
	bar->qual=al.Qualities;
		
	/*option fields*/
	//string stroptfileds="";
	//for(int i=10;i<vemp.size();i++)
	//	stroptfileds+=vemp[i];
}

void BamParse::setReadGroup(BamAlignment& al,BamAlignmentRecord*& bar)
{
	std::string rg;
	if ( al.GetTag("RG", rg) ) 
	{
		bar->tags=rg;
	} 
	else 
	{
		bar->tags="";
	}
}

bool cmpAlignmentRecord(BamAlignmentRecord* bar1, BamAlignmentRecord* bar2)
{
	if(bar1->rID==bar2->rID)
	{
		if(bar1->pos==bar2->pos)
			return bar1->qName < bar2->qName;
		else 
			return bar1->pos < bar2->pos;
	}
	else 
		return bar1->rID < bar2->rID;
}

/*
Description:
Sort all the alignment, in the priority of Chrom first, then read Name.
*/
void BamParse::sortAlignmentByChrom()
{
	sort(bam_aln_records.begin(),bam_aln_records.end(),cmpAlignmentRecord);
}


void BamParse::clearAll()
{
	int size=0;
	size=bam_head_records.size();
	for(int i=0;i<size;i++)
	{
		if(bam_head_records[i]!=NULL)
		{
			delete bam_head_records[i];
			bam_head_records[i]=NULL;
		}
	}
	vector<BamHeaderRecord*>().swap(bam_head_records);

	size=bam_aln_records.size();
	for(int i=0;i<size;i++)
	{
		if(bam_aln_records[i]!=NULL)
		{
			delete bam_aln_records[i];
			bam_aln_records[i]=NULL;
		}
	}
	vector<BamAlignmentRecord*>().swap(bam_aln_records);
}

int BamParse::getRefInfo()
{
	//BamReader reader; mapping quality
    if ( !reader.Open(filename) ) {
        cerr << "Bamtools ERROR: could not open input BAM file: " << filename << endl;
        return false;
    }
	const BamTools::RefVector& refs1 = reader.GetReferenceData();
	refs=refs1;

	reader.Close();
	return refs.size();
}

void BamParse::getChromNameLength(int refid, string& name, int& length)
{
	if(refid==-1) return; //unmapped

	assert(refid < refs.size()); // just a sanity check before accessing the vector 
	const string& refname = refs[refid].RefName;
	name=refname;
	const int len=refs[refid].RefLength;
	length=len;
}

//void BamParse::openWriter(string filename)
//{
//	sh=reader.GetHeader();
//	if(writer.IsOpen()==false)
//	{
//		bool bopen=writer.Open(filename,sh,refs);
//		if(bopen==false)
//			cerr<<"Cannot write to "<<filename<<endl;
//	}
//}
//	
//void BamParse::saveAlignment(BamAlignment al)
//{	
//	bool bw=writer.SaveAlignment(al);
//	if(bw==false)
//		cout<<"Save alignment error!!"<<endl;
//}
//
//void BamParse::closeWriter()
//{
//	writer.Close();
//}

