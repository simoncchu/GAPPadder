#include<iostream>
#include<string>
#include<fstream>
#include"insert_size.h"
#include"bam_parse.h"

using namespace std;

InsertSize::InsertSize()
{
	sfbam="";
	sfref="";
	mean=0.0;
	std_dev=0.0;
}

//
void InsertSize::estimate4Contigs(int min_contig_len)//estimate insert size for contigs (align reads to contigs) 
{
	//read in bam file
	BamParse bp(sfbam);
    bp.loadIndex();
	bp.openReader();

	//read in fai file
	string fref_fai=sfref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, itemp1, itemp2;
	long long  ctg_start;
	int istart=0, iend;
	int ctg_id=0;

	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp1>>itemp2)
	{//for each contig 
		if(ctg_lenth<min_contig_len)
			continue;
		
		
	}

	bp.closeReader();//close BamReader 
	fin_fai.close();
}