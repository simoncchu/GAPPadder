#include"Coverage.h"
#include"bam_parse.h"
#include"Alignment.h"
#include"public_parameters.h"
#include<iostream> 
#include<fstream>

Coverage::Coverage()
{
	this->region_lenth=0;
	this->read_len=102;
}

double Coverage::calcRegionCoverage(BamParse& bp, int& covered_length)
{
	covered_length = 0;
	double cov=0.0;
	if(this->region_lenth<=0) return cov;

	//cout<<"region lenth "<<region_lenth;///////////////////////////////////////////////////////////////////////
	int* cov_bases=new int[region_lenth+1];
	for(int i=0;i<=region_lenth;i++) cov_bases[i]=0;
	int vbamsize=bp.bam_aln_records.size();//number of reads
	
	for(int i=0;i<vbamsize;i++)
	{		
		int map_pos=bp.bam_aln_records[i]->pos;
		Alignment alnmt;
		alnmt.setBar(bp.bam_aln_records[i]);
		
		//first check whether is a qualified read.
		if((alnmt.isDuplicate()==true) || (alnmt.isPrimaryAlign() == false) || (alnmt.passQualityCK()==false))
		{
			continue;
		}
		
		//then check reads type 
		int aln_type=alnmt.getReadType();
		if(aln_type==READ_TYPE_FULLMAP)
		{
			for(int i=0;i<READ_LENGTH;i++)
			{
				int itmp=map_pos+i;
				//if(itmp >= this->region_lenth) cout<<"full map bypass "<<map_pos<<" "<<alnmt.getCigar()<<" "<<itmp<<" "<<region_lenth<<endl;/////////////////////////////////
				if(itmp>region_lenth)
				{
					break;
				}
				cov_bases[itmp]++;
			}
		}
		else if(aln_type==READ_TYPE_CLIP)
		{
			int cigar_size=bp.bam_aln_records[i]->cigar.size();
			//first check whether a qualified read 
			int len_read=0;
			int len_region=0;
			for(int k=0;k<cigar_size;k++)
			{
				string scg=bp.bam_aln_records[i]->cigar[k].first;
				int icg=bp.bam_aln_records[i]->cigar[k].second;
				
				if(scg=="M" || scg=="S" || scg=="H" || scg=="I")
					len_read+=icg;

				if(scg=="M" || scg!="D")
					len_region+=icg;
			}
			len_region+=map_pos;

			//if(len_read > READ_LENGTH) continue;////////////////////////////for now bypass this kind of reads 
			//if(len_region>region_lenth) continue;
			
			int j=0;
			bool bbrk=false;
			for(int k=0;k<cigar_size;k++)
			{
				string scg=bp.bam_aln_records[i]->cigar[k].first;
				int icg=bp.bam_aln_records[i]->cigar[k].second;

				if(scg=="S" || scg=="H" || scg=="I")
					continue;
				else if(scg=="M")
				{
					for(int ii=0;ii<icg;ii++)
					{
						int itmp=map_pos+j;
						if(itmp > region_lenth) 
						{
							//cout<<"clip overflow "<<bp.bam_aln_records[i]->qName<<" "<<map_pos<<" "<<alnmt.getCigar()<<" "<<itmp<<" "<<region_lenth<<" "<<len_read<<endl;
							bbrk=true;
							break;
						}
						//if(itmp>=region_lenth) cout<<itmp<<endl;////////////////////////////////////////////////////////////////////////////////////////
						cov_bases[itmp]++;
						j++;
					}
				}
				else if(scg=="D")
				{
					j+=icg;
				}

				if(bbrk==true)
					break;
			}
			
		}
		else
		{
			//cout<<READ_LENGTH<<" "<<bp.bam_aln_records[i]->qName<<" "<<alnmt.getCigar()<<" "<<bp.bam_aln_records[i]->mapQ<<endl; ////////////////////////////////////////////////////////////
		}
	}

	
	ofstream fout_base_cov;
	fout_base_cov.open("base_coverage.txt",ofstream::app);

	//calc average coverage 
	long long total_cov=0;
	for(int i=0;i<region_lenth;i++)
	{
//cout<<i<<" "<<cov_bases[i]<<endl;
fout_base_cov<<region_lenth<<":"<<i<<" "<<cov_bases[i]<<endl;

		total_cov+=cov_bases[i];
		if(cov_bases[i]>0) covered_length+=1;
	}
	cov=(double)total_cov/(double)region_lenth;
	
	fout_base_cov.close();

	if(cov_bases!=NULL)
	{
		delete[] cov_bases;
		cov_bases=NULL;
	}

	//std::cout<<"delete "<<endl;/////////////////////////////////////////////////////////////////////////////////////////////
	return cov;
}

//calc coverage with cutoff that used to filter out reads that only have small part is aligned
double Coverage::calcRegionCoverageWithCutOff(BamParse& bp, double cutoff)
{
	double cov=0.0;
	if(this->region_lenth<=0) return cov;

	//cout<<"region lenth "<<region_lenth;/////////////////////////////////////////////////////////////////////
	//int* cov_bases=new int[region_lenth+1];
	//for(int i=0;i<=region_lenth;i++) cov_bases[i]=0; 
	int vbamsize=bp.bam_aln_records.size();//number of reads 
	long cnt_total=0;
	for(int i=0;i<vbamsize;i++)
	{		
		int map_pos=bp.bam_aln_records[i]->pos;
		Alignment alnmt;
		alnmt.setBar(bp.bam_aln_records[i]);
		
		//first check whether is a qualified read.
		if((alnmt.isDuplicate()==true) || (alnmt.isPrimaryAlign() == false) || (alnmt.passQualityCK()==false))
		{
			cout<<READ_LENGTH<<" "<<bp.bam_aln_records[i]->qName<<" "<<alnmt.getCigar()<<" "<<bp.bam_aln_records[i]->mapQ<<endl;//////////////////////////
			continue;
		}

		int cigar_size=bp.bam_aln_records[i]->cigar.size();
		int cnt_read=0;
		int read_length=0;
		for(int k=0;k<cigar_size;k++)
		{
			string scg=bp.bam_aln_records[i]->cigar[k].first;
			int icg=bp.bam_aln_records[i]->cigar[k].second;
			read_length+=icg;
			if(scg=="M")
				cnt_read+=icg;
		}

		double ratio=(double)cnt_read/(double)READ_LENGTH;
		if(ratio>=cutoff)
			cnt_total+=cnt_read;
	}
	cov=(double)cnt_total/(double)region_lenth;
	return cov;
}

void Coverage::setRegionLenth(int len)
{
	this->region_lenth=len;
}

void Coverage::setReadLenth(int read_lenth)
{
	this->read_len=read_lenth;
}