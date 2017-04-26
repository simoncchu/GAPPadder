#include"refiner.h"
#include<fstream>
#include<string>
#include<vector>
#include<utility>
#include<iostream>
#include<algorithm>
#include <map>
#include"bam_parse.h"
#include"Alignment.h"
#include"public_parameters.h"
#include"public_func.h"
#include"contigs.h"
#include"fai_parser.h"
#include"Coverage.h"

#include "api/BamReader.h"
#include"api/BamWriter.h"
#include"api/SamHeader.h"
#include "api/BamAux.h"
using namespace BamTools;


using namespace std;

const int LENN=200;//length of the "N" between each two contigs. 
const int ISREVERSE=1;
const int ISFORWARD=0;

Refiner::Refiner(std::string fbam, std::string ffai, std::string fout_folder)
{
	this->fbam=fbam;
	this->ffai=ffai;
	this->fout_folder=fout_folder;
}


void Refiner::refineByReads(string fbam, string fref, string fout)
{
	//read in bam file 
	BamParse bp(fbam);
    bp.loadIndex();

	//read in fai file 
	string ffai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(ffai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, ctg_start, itemp;
	int ctg_id=0;

	vector<int> vid_rm; vid_rm.clear();
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig 
		int cnt_clip=0;
		int cnt_fullmap=0;

		bp.clearAll();//clear all the records
		bool bsignal=bp.parseAlignment(ctg_id,0,ctg_id,ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}

		//ofstream fout_contig;
		//string fctg_overall=fout_folder+sctg;
		//fout_contig.open(fctg_overall.c_str());

		/*ofstream fout_clip_reads;
		string fname_clip=fout_folder+PubFuncs::cvtInt2Str(ctg_id) + "_" + sctg + "_clipped_reads.txt";
		fout_clip_reads.open(fname_clip.c_str());*/

		/*ofstream fout_fullmap_reads;
		string fname_fullmap=fout_folder+PubFuncs::cvtInt2Str(ctg_id) + "_" + sctg + "_fully_mapped_reads.txt";
		fout_fullmap_reads.open(fname_fullmap.c_str());*/

		//read in reads
		int vbamsize=bp.bam_aln_records.size();//number of reads
		for(int i=0;i<vbamsize;i++)
		{
			Alignment alnmt;
			alnmt.setBar(bp.bam_aln_records[i]);
			
			int map_pos,rid,rnext_id,pnext,flag;
			map_pos=bp.bam_aln_records[i]->pos;//mapping position
			string seq=bp.bam_aln_records[i]->seq;//seqence of reads
			string qname=bp.bam_aln_records[i]->qName;//reads id
			rid=bp.bam_aln_records[i]->rID;
			rnext_id=bp.bam_aln_records[i]->rNextID;//ref id of mate-read 
			pnext=bp.bam_aln_records[i]->pNext;//mate read mapping position.
			flag=bp.bam_aln_records[i]->flag;
			std::string cigar=alnmt.getCigar();

			int aln_type=alnmt.getReadType();
			if(aln_type==READ_TYPE_CLIP)
			{
				cnt_clip++;

				int pos1,len1,pos2,len2, clip_pos1, clip_pos2;
				int clip_type=alnmt.getClipType(pos1,len1,pos2,len2);
				if(clip_type==READ_TYPE_LEFT_SOFTCLIP || clip_type==READ_TYPE_LEFT_HARDCLIP)
				{
					clip_pos1=map_pos;
					clip_pos2=-1;
				}
				else if(clip_type==READ_TYPE_RIGHT_SOFTCLIP || clip_type==READ_TYPE_RIGHT_HARDCLIP)
				{
					int len_mapped = READ_LENGTH - len2;
					clip_pos2 = map_pos + len_mapped;
					clip_pos1=-1;
				}
				else if(clip_type==READ_TYPE_BOTH_SOFTCLIP || clip_type==READ_TYPE_BOTH_HARDCLIP)
				{
					clip_pos1=map_pos;
					clip_pos2=map_pos + (READ_LENGTH - len1 - len2); 
				}

				/*fout_clip_reads<<clip_type<<" "<<qname<<" "<<flag<<" "<<sctg<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext \
						<<" "<<len1<<" "<<len2<<" "<<clip_pos1<<" "<<clip_pos2<<endl;*/
			}//end of if
			else if(aln_type==READ_TYPE_FULLMAP)
			{
				cnt_fullmap++;
				/*fout_fullmap_reads<<qname<<" "<<flag<<" "<<sctg<<" "<<map_pos<<" "<<cigar<<" "<<seq<<" "<<pnext<<endl;*/
			}
			else
			{
				//other type of reads
			}

		}//end of for 

		//fout_contig<<ctg_id<<" "<<sctg<<":"<<endl;
		//fout_contig<<"No. of clipped reads are "<<cnt_clip<<", and No. of fully mapped reads are "<<cnt_fullmap<<endl;
		//fout_contig.close(); 

		double cf_ratio=(double)cnt_fullmap/((double)cnt_clip+(double)cnt_fullmap);
		//cout<<cf_ratio<<endl;
		if(cf_ratio < this->cf_cutoff)
		{
			//cout<<"Ratio:"<<cf_ratio<<"    ID:"<<sctg<<endl;//output contig name 
			vid_rm.push_back(ctg_id);
		}

		/*fout_clip_reads.close();
		fout_fullmap_reads.close();*/
		ctg_id++;
	}//end of while
	
	//remove low-quality ones 
	rmCotigs(fref,vid_rm,fout);

	fin_fai.close();
}

//Concatenate all the contigs togethor and see the concatenated one as a whole chromosome.
void Refiner::refineByReadsCombinedContigs()
{
	//read in bam file
	BamParse bp(this->fbam);
    bp.loadIndex();

	//read in fai file 
	ifstream fin_fai;
	fin_fai.open(ffai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, ctg_start, itemp;
	int istart=0, iend;
	int ctg_id=0;
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig 
		int cnt_clip=0;
		int cnt_fullmap=0;	
		iend=istart + ctg_lenth - 1;

		bp.clearAll();//clear all the records 
		bool bsignal=bp.parseAlignment(0,istart,0,iend);//read in region [start, end] of chromosome id (start from 0) 
		if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}

		ofstream fout_contig;
		string fctg_overall=fout_folder+sctg;
		fout_contig.open(fctg_overall.c_str());

		ofstream fout_clip_reads;
		string fname_clip=fout_folder+PubFuncs::cvtInt2Str(ctg_id) + "_" + sctg + "_clipped_reads.txt";
		fout_clip_reads.open(fname_clip.c_str());

		ofstream fout_fullmap_reads;
		string fname_fullmap=fout_folder+PubFuncs::cvtInt2Str(ctg_id) + "_" + sctg + "_fully_mapped_reads.txt";
		fout_fullmap_reads.open(fname_fullmap.c_str());

		//read in reads 
		int vbamsize=bp.bam_aln_records.size();//number of reads 
		for(int i=0;i<vbamsize;i++)
		{
			Alignment alnmt;
			alnmt.setBar(bp.bam_aln_records[i]);
			
			int map_pos,rid,rnext_id,pnext,flag;
			map_pos=bp.bam_aln_records[i]->pos;//mapping position
			string seq=bp.bam_aln_records[i]->seq;//seqence of reads
			string qname=bp.bam_aln_records[i]->qName;//reads id
			rid=bp.bam_aln_records[i]->rID;
			rnext_id=bp.bam_aln_records[i]->rNextID;//ref id of mate-read 
			pnext=bp.bam_aln_records[i]->pNext;//mate read mapping position.
			flag=bp.bam_aln_records[i]->flag;
			std::string cigar=alnmt.getCigar();

			int aln_type=alnmt.getReadType();
			if(aln_type==READ_TYPE_CLIP)
			{
				cnt_clip++;

				int pos1,len1,pos2,len2, clip_pos1, clip_pos2;
				int clip_type=alnmt.getClipType(pos1,len1,pos2,len2);
				if(clip_type==READ_TYPE_LEFT_SOFTCLIP || clip_type==READ_TYPE_LEFT_HARDCLIP)
				{
					clip_pos1=map_pos;
					clip_pos2=-1;
				}
				else if(clip_type==READ_TYPE_RIGHT_SOFTCLIP || clip_type==READ_TYPE_RIGHT_HARDCLIP)
				{
					int len_mapped = READ_LENGTH - len2;
					clip_pos2 = map_pos + len_mapped;
					clip_pos1=-1;
				}
				else if(clip_type==READ_TYPE_BOTH_SOFTCLIP || clip_type==READ_TYPE_BOTH_HARDCLIP)
				{
					clip_pos1=map_pos;
					clip_pos2=map_pos + (READ_LENGTH - len1 - len2); 
				}

				fout_clip_reads<<clip_type<<" "<<qname<<" "<<flag<<" "<<sctg<<" "<<(map_pos-istart)<<" "<<cigar<<" "<<seq<<" "<<pnext \
						<<" "<<len1<<" "<<len2<<" "<<(clip_pos1-istart)<<" "<<(clip_pos2-istart)<<endl;
			}//end of if
			else if(aln_type==READ_TYPE_FULLMAP)
			{
				cnt_fullmap++;
				fout_fullmap_reads<<qname<<" "<<flag<<" "<<sctg<<" "<<(map_pos-istart)<<" "<<cigar<<" "<<seq<<" "<<pnext<<endl;
			}
			else
			{
				//other type of reads 
			}

		}//end of for

		fout_contig<<ctg_id<<" "<<sctg<<":"<<endl;
		fout_contig<<"No. of clipped reads are "<<cnt_clip<<", and No. of fully mapped reads are "<<cnt_fullmap<<endl;
		fout_contig.close();

		double cf_ratio=(double)cnt_fullmap/((double)cnt_clip+(double)cnt_fullmap);
		//cout<<cf_ratio<<endl;
		if(cf_ratio < this->cf_cutoff)
		{
			cout<<"Ratio:"<<cf_ratio<<"    ID:"<<sctg<<endl;//output contig name
		}

		fout_clip_reads.close();
		fout_fullmap_reads.close();

		istart+=ctg_lenth;
		istart+=LENN;
		ctg_id++;
	}//end of while

	fin_fai.close();
}

void Refiner::setCutOff(double cf_cutoff)
{
	this->cf_cutoff=cf_cutoff;
}

void Refiner::setThreshold(int threshold)
{
	this->threshold=threshold;
}

void Refiner::setContigFile(std::string fcontig)
{
	this->fcontig=fcontig;
}

void Refiner::setReadLength(int readlen)
{
	this->read_length=readlen;
}


void Refiner::removeRepeatsOfTwoContigSets(std::string fbam, std::string fref, std::string bam_fasta, std::string fnew_contig_fa)
{
	//read in bam file
	BamParse bp(fbam);
    bp.loadIndex();
	
	string sctg;
	int ctg_lenth, ctg_start, itemp;
	int istart=0, iend;
	//first build a map
	string sbam_fa_fai=bam_fasta+".fai";
	ifstream fbam_fai;
	fbam_fai.open(sbam_fa_fai.c_str());
	map<string,pair<int,int> > bam_mfa;
	int i=0;
	while(fbam_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{
		bam_mfa[sctg]=std::make_pair(i,ctg_lenth);
		i++;
	}
	fbam_fai.close();
	
	//count reads contig by contig 
	//read in fai file 
	string fref_fai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());
	vector<string> vid_dup;
	vid_dup.clear();
	int ctg_id=0;
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig
		bp.clearAll();//clear all the records 
		bool bsignal=bp.parseAlignment(ctg_id,0,ctg_id,ctg_lenth);//read in region [start, end] of chromosome id (start from 0) 
		if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}

		//read in reads 
		int vbamsize=bp.bam_aln_records.size();//number of reads 
		for(int i=0;i<vbamsize;i++)
		{
			Alignment alnmt;
			alnmt.setBar(bp.bam_aln_records[i]);
			string qname=bp.bam_aln_records[i]->qName;//reads id
			
			/*int aln_type=alnmt.getReadType(bam_mfa[qname].second);
			if(aln_type==READ_TYPE_FULLMAP)
			{
				vid_dup.push_back(qname);
				//cout<<qname<<endl; 
			}*/

			bool is_fully_map = alnmt.isFullyMapped(bam_mfa[qname].second);
			if(is_fully_map==true)
			{
				vid_dup.push_back(qname);
				//cout<<qname<<endl; 
			}

		}
		ctg_id++;
	}
	fin_fai.close();

	//remove the duplicate ones from bam_fasta file
	//find the index of those need to be removed
	vector<int> v_index_dup;
	v_index_dup.clear();
	for(int i=0;i<vid_dup.size();i++)
	{
		v_index_dup.push_back(bam_mfa[vid_dup[i]].first);
	}
	sort(v_index_dup.begin(),v_index_dup.end());//sort
	vector<int> v_index_dup2;
	v_index_dup2.clear();
	int vn=v_index_dup.size();
	if(vn>0)
	{
		int ivlast=v_index_dup[0];
		for(int i=1;i<vn;i++)
		{
			if(ivlast!=v_index_dup[i])
				v_index_dup2.push_back(ivlast);
			ivlast=v_index_dup[i];
		}
		v_index_dup2.push_back(ivlast);
	}
	//remove according to index 
	rmCotigs(bam_fasta, v_index_dup2, fnew_contig_fa);
}

/*
Input:
	bam_fastq, original file;
	v_index_dup, the ids that need to be removed from the original file;
Output:
	fnew_contig_fa, the generated new file
*/
void Refiner::rmCotigs(std::string bam_fasta, std::vector<int>& v_index_dup, std::string fnew_contig_fa)
{
	ifstream fin_bam_fa;
	ofstream fout_new_contig;

	if(v_index_dup.size()==0)
	{	
		if(bam_fasta==fnew_contig_fa) //may have problem here, because for example "./fa.test"!="fa.test", but they are same file.//////////////////////////////////////
			return;
		else
		{
			fin_bam_fa.open(bam_fasta.c_str());
			fout_new_contig.open(fnew_contig_fa.c_str());
			fout_new_contig<<fin_bam_fa.rdbuf();//copy the file 
			fin_bam_fa.close();
			fout_new_contig.close();
			return;
		}
	}
	else
	{
		fin_bam_fa.open(bam_fasta.c_str());
		string fa_id="", fa_seq="";
		vector<std::pair<string,string> > vcontigs;
		vcontigs.clear();

		string sline="";
		bool bfirst=true;
		while(std::getline(fin_bam_fa,sline))
		{
			if(sline.length()>0 && sline[0]=='>')
			{
				if(bfirst==true)
				{
					bfirst=false;
					fa_id=sline;
					continue;
				}
				vcontigs.push_back(std::make_pair(fa_id, fa_seq));
				fa_id=sline;
				fa_seq="";
			}
			else
			{
				fa_seq+=sline;
				fa_seq+="\n";
			}
		}
		vcontigs.push_back(std::make_pair(fa_id, fa_seq));//the last one 
		fin_bam_fa.close();

		fout_new_contig.open(fnew_contig_fa.c_str());
		int vsize=vcontigs.size();
		bool* bsignal=new bool[vsize];
		for(int i=0;i<vsize;i++)
			bsignal[i]=true;
		for(int i=0;i<v_index_dup.size();i++)
			bsignal[v_index_dup[i]]=false;

		for(int i=0;i<vsize;i++)
		{
			if(bsignal[i]==true)
			{
				fout_new_contig<<vcontigs[i].first<<endl; 
				fout_new_contig<<vcontigs[i].second;
			}
			
		}
		delete[] bsignal;
		fout_new_contig.close();
	}
}

bool cmp_vfa(const std::pair<string,int>& p1, const pair<string,int>& p2)
{
	if(p1.first==p2.first)
		return p1.second<p2.second;
	else
		return p1.first<p2.first;
}

void Refiner::removeRepeatsOfOneContigSet(std::string fbam, std::string fref, std::string fnew)
{
	//read in bam file
	BamParse bp(fbam);
    bp.loadIndex();
	
	//read in fai file 
	string fref_fai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, ctg_start, itemp;
	int istart=0, iend;
	int ctg_id=0;
	
	vector<std::pair<string,int> > vfa;
	vfa.clear();
	map<string,pair<int,int> > mfa;
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig 
		mfa[sctg]=std::make_pair(ctg_id,ctg_lenth);
		vfa.push_back(std::make_pair(sctg,ctg_id));
		ctg_id++;
	}
	fin_fai.close();

	//remove duplicates contained in another 
	fin_fai.open(fref_fai.c_str());
	vector<int> vid_rm;
	ctg_id=0;
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig
		bp.clearAll();//clear all the records 
		bool bsignal=bp.parseAlignment(ctg_id,0,ctg_id,ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}

		//read in reads 
		int vbamsize=bp.bam_aln_records.size();//number of reads
		
		//cout<<"bam file size: "<<vbamsize<<endl;//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		for(int i=0;i<vbamsize;i++)
		{
			Alignment alnmt;
			alnmt.setBar(bp.bam_aln_records[i]);
			string qname=bp.bam_aln_records[i]->qName;//reads id
			
			int rid=bp.bam_aln_records[i]->rID;//ref id
			string rname=vfa[rid].first;//bp.bam_aln_records[i]->rName;//ref name
			
			bool is_fully_map=alnmt.isFullyMapped(mfa[qname].second);
			if((is_fully_map == true) && (qname!=rname))
			{
				//cout<<qname<<" "<<rname<<" "<<mfa[qname].first<<" "<<mfa[rname].first<<endl;//report pairs ////////////////////////////////////////
				//for further work, for each pair, remove the short one...
				int iq=mfa[qname].second;
				int ir= mfa[rname].second;
				if(iq==ir)
				{
					//cout<<"equal......"<<endl;
					if(qname<rname)
					{
						//cout<<"equal remove "<<qname<<endl;///////////////////////////////////////////////////
						vid_rm.push_back(mfa[qname].first);
					}
				}
				else
				{
					//cout<<"insert some..."<<qname<<" "<<rname<<" "<<alnmt.getCigar()<<endl;//////////////////////////////////////////////
					vid_rm.push_back(mfa[qname].first);
				}
			}

		}
		ctg_id++;
	}
	fin_fai.close();

	//cout<<"remove duplicate and output"<<endl;/////////////////////////////////////////////////////////////////////////////////
	//remove duplicate and output 
	sort(vid_rm.begin(),vid_rm.end());
	vector<int> vid_rm2;
	vid_rm2.clear();
	int vn=vid_rm.size();
	if(vn>0)
	{
		int ivlast=vid_rm[0];
		for(int i=1;i<vn;i++)
		{
			if(ivlast!=vid_rm[i])
				vid_rm2.push_back(ivlast);
			ivlast=vid_rm[i];
		}
		vid_rm2.push_back(ivlast);
	}

	rmCotigs(fref, vid_rm2, fnew);
}

void Refiner::removeContainedContigs(std::string fbam, std::string fref, std::string fnew)
{
	//read in bam file
	BamParse bp(fbam);
	bp.loadIndex();

	//read in fai file 
	string fref_fai = fref + ".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, ctg_start, itemp;
	int istart = 0, iend;
	int ctg_id = 0;

	vector<std::pair<string, int> > vfa;
	vfa.clear();
	map<string, pair<int, int> > mfa;
	while (fin_fai >> sctg >> ctg_lenth >> ctg_start >> itemp >> itemp)
	{//for each contig 
		mfa[sctg] = std::make_pair(ctg_id, ctg_lenth);
		vfa.push_back(std::make_pair(sctg, ctg_id));
		ctg_id++;
	}
	fin_fai.close();


	vector<int> vid_rm;
	map<int, int> mrm;//<id,cnt> 
	bp.clearAll();//clear all the records 	
	bool bsignal = bp.parseAlignment(-1, -1, -1, -1);
	if (bsignal == false)
	{
		std::cout << "Cannot parse bam file " << this->fbam << std::endl;
		return;
	}

	//read in reads 
	int vbamsize = bp.bam_aln_records.size();//number of reads 
//cout<<"Number of alignments: "<<vbamsize<<endl;/////////////////////////////////////////////////////////////////////////////////////////////////////	
	for (int i = 0; i < vbamsize; i++)
	{
		Alignment alnmt;
		alnmt.setBar(bp.bam_aln_records[i]);
		string qname = bp.bam_aln_records[i]->qName;//reads id

		int rid = bp.bam_aln_records[i]->rID;//ref id
		string rname = vfa[rid].first;//bp.bam_aln_records[i]->rName;//ref name 

		if (qname == rname) continue;

		bool is_fully_map=alnmt.isFullyMapped(mfa[qname].second);
		if (is_fully_map == true)
		{
//cout << mfa[qname].first << endl;////////////////////////////////////////////////////////////////////////////
			if (mrm.count(mfa[qname].first)>0)
			{//already exist 
				mrm[mfa[qname].first]++;
			}
			else
			{
				mrm[mfa[qname].first] = 0;
				vid_rm.push_back(mfa[qname].first);
			}
		}

	}
	rmCotigs(fref, vid_rm, fnew);
}

//remove the duplicate ones with/withou contained ones
void Refiner::removeDupRepeatsOfOneContigSet(std::string fbam, std::string fref, std::string fnew, double cutoff_ratio, bool brm_cntn)
{
	//read in bam file
	BamParse bp(fbam);
    bp.loadIndex();
	
	//read in fai file 
	string fref_fai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, ctg_start, itemp;
	int istart=0, iend;
	int ctg_id=0;
	
	vector<std::pair<string,int> > vfa;
	vfa.clear();
	map<string,pair<int,int> > mfa;
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig 
		mfa[sctg]=std::make_pair(ctg_id,ctg_lenth);
		vfa.push_back(std::make_pair(sctg,ctg_id));
		ctg_id++;
	}
	fin_fai.close();

	
	vector<int> vid_rm;
	map<int,int> mrm;//<id,cnt> 
	bp.clearAll();//clear all the records 	
	bool bsignal=bp.parseAlignment(-1,-1,-1,-1);
	if(bsignal==false)
	{
		std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
		return;
	}

	//read in reads 
	int vbamsize=bp.bam_aln_records.size();//number of reads 
	//cout<<"Number of alignments: "<<vbamsize<<endl;/////////////////////////////////////////////////////////////////////////////////////////////////////	
	for(int i=0;i<vbamsize;i++)
	{
		Alignment alnmt;
		alnmt.setBar(bp.bam_aln_records[i]);
		string qname=bp.bam_aln_records[i]->qName;//reads id	
		int rid=bp.bam_aln_records[i]->rID;//ref id
		string rname=vfa[rid].first;//bp.bam_aln_records[i]->rName;//ref name 
		if(brm_cntn==false)
		{//remove duplicate 
			bool is_fully_map=alnmt.isFullyMapped(mfa[qname].second);
			if((is_fully_map == true) && (qname > rname))
			{
				//cout<<qname<<" "<<rname<<" "<<mfa[qname].first<<" "<<mfa[rname].first<<endl;//report pairs ////////////////////////////////////////
				//for further work, for each pair, remove the short one...
				int iq=mfa[qname].second; //length of the reads
				int ir= mfa[rname].second;//length of the ref 
				if(iq==ir) 
				{
					if(mrm.count(mfa[qname].first) > 0)
					{//already exist
						mrm[mfa[qname].first]++;
					}
					else
					{
						mrm[mfa[qname].first]=0;
						vid_rm.push_back(mfa[qname].first);
					}
				}
				else
				{
					//only when length are similar, they will be removed 
					int idiff=0;
					int imin=0;
					if(iq>ir)
					{
						idiff=iq-ir;
						imin=ir;
					}
					else
					{
						idiff=ir-iq;
						imin=iq;
					}

					if(((double)idiff/(double)imin) <= (1.0-cutoff_ratio))
					{
						if(mrm.count(mfa[qname].first)>0)
						{//already exist
							mrm[mfa[qname].first]++;
						}
						else
						{
							mrm[mfa[qname].first]=0;
							vid_rm.push_back(mfa[qname].first);
						}							 
					}
				}
			}
		}
		else 
		{//remove contained
			if(qname==rname) continue;

			bool is_perfect_map=alnmt.isPerfectMapped(mfa[qname].second);
			//bool is_perfect_map=alnmt.isFullyMapped(mfa[qname].second);
			if(is_perfect_map==true)
			{
				if(mrm.count(mfa[qname].first)>0)
				{//already exist 
					mrm[mfa[qname].first]++;
				}
				else
				{
					mrm[mfa[qname].first]=0;
					vid_rm.push_back(mfa[qname].first);
				}
			}
		}

	}
	//cout<<"remove duplicate and output"<<endl;/////////////////////////////////////////////////////////////////////////////////
	//remove duplicate and output
	
	/*sort(vid_rm.begin(),vid_rm.end());
	vector<int> vid_rm2;
	vid_rm2.clear();
	int vn=vid_rm.size();
	if(vn>0)
	{
		int ivlast=vid_rm[0];
		for(int i=1;i<vn;i++)
		{
			if(ivlast!=vid_rm[i])
				vid_rm2.push_back(ivlast);
			ivlast=vid_rm[i];
		}
		vid_rm2.push_back(ivlast);
	}*/
	rmCotigs(fref, vid_rm, fnew);
}

/* 
Max non-overlapping coverage of a benchmark repeat. The higher, the better.
#Total fraction of benchmark repeats that are covered by some contigs. 
Total coverage of contigs mapped to a benchmark repeat. If higher than 1, then smaller is better; if much less than 1, then larger is better.
#Longest single contig that is mapped to benchmark (excluding indel/clipping).
Total fraction of bases of repeat contigs that are not mapped to a benchmark. 
#Average coverage of benchmark contigs where there are mapped contigs (i.e. excluding positions where there is no mapped repeat contigs). 
*/
bool cmpSgmt(std::pair<int,int>p1, std::pair<int,int>p2)
{
	if(p1.second < p2.second)
		return true;
	else
		return false;
}

void Refiner::getLongestUncover(std::vector<std::pair<int,int> >& vsgmts, int& luncover)
{
	luncover=vsgmts[0].second-vsgmts[0].first+1;
	int cur_end=vsgmts[0].second;
	for(int i=1;i<vsgmts.size();i++)
	{
		if(vsgmts[i].first<cur_end)
			continue;
		cur_end=vsgmts[i].second;
		luncover+=(vsgmts[i].second-vsgmts[i].first+1);
	}
}

void Refiner::evaluateWithBenchmark(std::string fbam, std::string fref, double cutoff_map_fraction)
{
//cout<<fbam<<" "<<fref<<endl;
	//read in bam file 
	BamParse bp(fbam);
    bp.loadIndex();
	
	//read in fai file 
	string fref_fai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str(),ifstream::in);

	//count reads contig by contig 
	string sctg, slength;
	int ctg_lenth, ctg_start, itemp;
	int istart=0, iend;
	int ctg_id=0;
	
	vector<std::pair<string,int> > vfa;
	vfa.clear();
	map<string,pair<int,int> > mfa;
	while(!fin_fai.eof())
	{
		string aline="";
		getline(fin_fai,aline);
//cout<<aline<<endl;
		std::stringstream ss(aline);	
		std::getline(ss, sctg, '\t');
		std::getline(ss, slength, '\t');
		ctg_lenth=PubFuncs::cvtStr2Int(slength);
//cout<<sctg<<" "<<ctg_lenth<<" "<<ctg_start<<endl;
		mfa[sctg]=std::make_pair(ctg_id,ctg_lenth);
		vfa.push_back(std::make_pair(sctg,ctg_id));
		ctg_id++;
	}

//	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
//	{//for each contig 
//cout<<sctg<<" "<<ctg_lenth<<" "<<ctg_start<<" "<<itemp<<endl;
//		mfa[sctg]=std::make_pair(ctg_id,ctg_lenth);
//		vfa.push_back(std::make_pair(sctg,ctg_id));
//		ctg_id++;
//	}
	fin_fai.close();
	
	ofstream fout_cover_ratio;
	string fout_cr_name=fbam+".statistic.txt";
	fout_cover_ratio.open(fout_cr_name.c_str());

	ofstream fout_cov_ratio_table;
	string fout_cr_name_table=fbam+".statistic.table.txt";
	fout_cov_ratio_table.open(fout_cr_name_table.c_str());
	fout_cov_ratio_table<<"RefName\tLength\ttotal-covered-region\ttotal_covered_region_ratio\tLongest-signle-cover\tLongest-signle-cover-ratio\t";
	//fout_cov_ratio_table<<"Longest-uncover\tLongest-uncover-ratio\ttotal-mapped-bases\tAverage-coverage-of-mapped-region"<<endl;
	fout_cov_ratio_table<<"total-mapped-bases\tAverage-coverage-of-mapped-region"<<endl;

	fin_fai.open(fref_fai.c_str());
	int total_covered=0;
	int total_ref=0;
	ctg_id=0;
	
	while(!fin_fai.eof()) //while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig
		string aline="";
		getline(fin_fai,aline);
		std::stringstream ss(aline);
		std::getline(ss, sctg, '\t');
		std::getline(ss, slength, '\t');
		ctg_lenth=PubFuncs::cvtStr2Int(slength);
		string sstart_pos="";
		std::getline(ss, sstart_pos, '\t');
		ctg_start=PubFuncs::cvtStr2Int(sstart_pos);

//cout<<sctg<<" "<<ctg_lenth<<" "<<ctg_start<<endl;
		bp.clearAll();//clear all the records 
		
		bool bsignal=bp.parseAlignment(ctg_id,0,ctg_id,ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}
		
		//read in reads 
		int vbamsize=bp.bam_aln_records.size();//number of reads
		if(vbamsize<=0) 
		{
			ctg_id++;
			continue;
		}

		//get the longest alignment 
		int max_match=0;
		int rid;
		int* cov_bases=new int[ctg_lenth+1];//base coverage of fa 
		for(int i=0;i<=ctg_lenth;i++)
			cov_bases[i]=0;

		vector<std::pair<int,int> > vsgmts;
		vsgmts.clear();
//cout<<"vbamsize: "<<vbamsize<<endl;
		for(int i=0;i<vbamsize;i++)
		{
			Alignment alnmt;
			alnmt.setBar(bp.bam_aln_records[i]);
			string qname=bp.bam_aln_records[i]->qName;//reads id 

			if((alnmt.getReadType()==READ_PAIR_MAP_TYPE_01) || (alnmt.getReadType()==READ_PAIR_MAP_TYPE_00))//read unmapped
			{
				//get the length of the unmapped contig 
				continue;
			}

			int map_pos=bp.bam_aln_records[i]->pos;
			rid=bp.bam_aln_records[i]->rID;//ref id 

			int mlen=0;
			int cigar_size=bp.bam_aln_records[i]->cigar.size();
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
					mlen+=icg;

					for(int ii=0;ii<icg;ii++)
					{
						int itmp=map_pos+j;
						if(itmp > ctg_lenth) 
						{
							bbrk=true;
							break;
						}
						
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
			
			vsgmts.push_back(std::make_pair(map_pos, map_pos+mlen-1));

			if(((double)mlen/(double)ctg_lenth) < cutoff_map_fraction )
			{
				continue;
			}

			if(mlen>max_match)
				max_match=mlen;
		}//end of bam alignment loop
		
		sort(vsgmts.begin(),vsgmts.end(),cmpSgmt);
		int longest_uncover=0;
		getLongestUncover(vsgmts,longest_uncover);

		int total_covered_region=0;//total length of the covered region
		int total_mapped_bases=0;
		for(int r=0;r<ctg_lenth;r++)
		{
			if(cov_bases[r]>0)
			{
				total_covered_region++;
				total_mapped_bases+=cov_bases[r];
			}
		}
		ctg_id++;
		
		if(max_match<=0) continue;
		string rname=vfa[rid].first;//ref name
		int ir= mfa[rname].second;//length of the ref
		fout_cover_ratio<<"------------------------------------------------------------------------------------------"<<endl;
		fout_cover_ratio<<rname<<" "<<ir<<endl;
		fout_cover_ratio<<"Total covered: "<<total_covered_region<<" "<<(double)total_covered_region/(double)ir<<endl;
		fout_cover_ratio<<"Longest signle cover: "<<max_match<<" "<<(double)max_match/(double)ir<<endl;
		fout_cover_ratio<<"Total length of all uncovered sgmts:"<<longest_uncover<<" "<<(double)longest_uncover/(double)ir<<endl;
		fout_cover_ratio<<"Average coverage of mapped region: "<<total_mapped_bases<<" "<<(double)total_mapped_bases/(double)total_covered_region<<endl;
		fout_cover_ratio<<"------------------------------------------------------------------------------------------"<<endl;
		
		//fout_cov_ratio_table<<rname<<"\t"<<ir<<"\t"<<total_covered_region<<"\t"<<(double)total_covered_region/(double)ir<<"\t"<<max_match<<"\t"<<(double)max_match/(double)ir<<"\t"\
		//					<<longest_uncover<<"\t"<<(double)longest_uncover/(double)ir<<"\t"<<total_mapped_bases<<"\t"<<(double)total_mapped_bases/(double)total_covered_region<<endl;
		fout_cov_ratio_table<<rname<<"\t"<<ir<<"\t"<<total_covered_region<<"\t"<<(double)total_covered_region/(double)ir<<"\t"<<max_match<<"\t"<<(double)max_match/(double)ir<<"\t"\
							<<total_mapped_bases<<"\t"<<(double)total_mapped_bases/(double)total_covered_region<<endl;

		total_covered+=max_match;
		total_ref+=ir;
		
		if(cov_bases!=NULL)
		{
			delete[] cov_bases;
			cov_bases=NULL;
		}
	}

	//fout_cover_ratio<<total_ref<<" "<<total_covered<<" "<<(double)total_covered/(double)total_ref<<endl;
	fout_cover_ratio.close();
	fout_cov_ratio_table.close();
	fin_fai.close();
	
}

void Refiner::gnrtUniqueFa(std::string fa_old, std::string fa_new)
{
	//read in fai file 
	string fref_fai=fa_old+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	vector<std::pair<string,int> > vfa;
	string sctg;
	int ctg_lenth, ctg_start, itemp;
	int istart=0, iend;
	int ctg_id=0;
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
	{//for each contig 
		vfa.push_back(std::make_pair(sctg,ctg_id));
		ctg_id++;
	}
	fin_fai.close();

	vector<int> vid_rm;
	vector<int> vid_rm2;
	vid_rm2.clear();
	//remove those duplicate themselves 
	if(vfa.size()>0)
	{
		sort(vfa.begin(),vfa.end(),cmp_vfa);
		string vfa_slast=vfa[0].first;
		for(int i=1;i<vfa.size();i++)
		{
			if(vfa[i].first==vfa_slast)
			{
				vid_rm.push_back(vfa[i].second);
			}
			vfa_slast=vfa[i].first;
		}

		//remove duplicate and output 
		sort(vid_rm.begin(),vid_rm.end());
		int vn=vid_rm.size();
		if(vn>0)
		{
			int ivlast=vid_rm[0];
			for(int i=1;i<vn;i++)
			{
				if(ivlast!=vid_rm[i])
					vid_rm2.push_back(ivlast);
				ivlast=vid_rm[i];
			}
			vid_rm2.push_back(ivlast);
		}

		///////////////////////////////////////
		//ofstream fout_test1;//////////////////////////////////////////////////////////////////////////////////////
		//fout_test1.open("test_test1.txt"); 
		//for(int i=0;i<vid_rm.size();i++)
		//{
		//	fout_test1<<vid_rm[i]<<endl;
		//}
		//fout_test1.close();///////////////////////////////////////////////////////////////////////////////////////
		//
		//ofstream fout_test2;//////////////////////////////////////////////////////////////////////////////////////
		//fout_test2.open("test_test2.txt");
		//for(int i=0;i<vid_rm2.size();i++)
		//{
		//	fout_test2<<vid_rm2[i]<<endl;
		//}
		//fout_test2.close();///////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////

		rmCotigs(fa_old, vid_rm2, fa_new);
	}
	
}

bool cmpContigs(pair<pair<Contigs,Contigs>,double> pcd1, pair<pair<Contigs,Contigs>,double> pcd2)
{
	pair<Contigs,Contigs> pc1=pcd1.first;
	pair<Contigs,Contigs> pc2=pcd2.first;
	if(pc1.first.getContigName() == pc2.first.getContigName())
	{
		if(pc1.second.getContigName() == pc2.second.getContigName())
		{
			if(pc1.first.getOrientation() == pc2.first.getOrientation())
				return pc1.second.getOrientation() < pc2.second.getOrientation();
			else
				return pc1.first.getOrientation() < pc2.first.getOrientation();
		}
		else
			return pc1.second.getContigName() < pc2.second.getContigName();
	}
	else
	{
		return pc1.first.getContigName() < pc2.first.getContigName();
	}
}

void Refiner::cntContigLinkage(std::string fbam, std::string fref, std::string fcontig_info)// cnt number of paired-end reads that link each two contigs
{
	//read in bam file 
	BamParse bp(fbam);
    bp.loadIndex();
	bp.openReader();

	//read in fai file 
	string fref_fai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, itemp1, itemp2;
	long long  ctg_start;
	int istart=0, iend;
	int ctg_id=0;

	//vector<pair<Contigs,Contigs> > vcontig_pairs; 
	vector<pair<pair<Contigs,Contigs>,double> > vcontig_pairs;
	vcontig_pairs.clear();
	//vector<double> vcp_distance;//contig pair distance, relative to vcontig_pairs 
	//vcp_distance.clear();

	int cnt_full_map_pairs=0;
	FaiParser faipsr(fref_fai);
	faipsr.parseFai();

	vector<pair<string,int> > vchroms;
	vchroms.clear();
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp1>>itemp2)
	{//for each contig 
		
		vchroms.push_back(std::make_pair(sctg,ctg_lenth));
	}
	fin_fai.close();

	//cout<<"No. of contigs is: "<<vchroms.size()<<endl;

	vector<double> vcoverage;
	vcoverage.clear();
	vector<int> vcovered_length;
	vcovered_length.clear();

	ofstream fout_connect;
	string fname_temp=fcontig_info+"_temp_ContigsConnections.txt";
	fout_connect.open(fname_temp.c_str());
	
	for(int k=0;k<vchroms.size();k++)
	{//for each contig 
		sctg=vchroms[k].first;
		ctg_lenth=vchroms[k].second;

		//cout<<sctg<<" "<<ctg_lenth<<endl;/////////////////////////////////////////////////////////////////////////////////////////////
		//ofstream fout_contig;//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//string sctg_path="temp/"+sctg+".txt";
		//fout_contig.open(sctg_path.c_str());///////////////////////////////////////////////////////////////////////////////////////////

		bp.clearAll();//clear all the records
		bp.dumpAlignments(ctg_id, 0, ctg_id, ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		//bool bsignal=bp.dumpAlignments(ctg_id,0,ctg_id,ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		/*if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}*/
		//cout<<"dumped..."<<endl;//////////////////////////////////////////////////////////////////////////

		Coverage coverage;
		coverage.setRegionLenth(ctg_lenth);
		coverage.setReadLenth(this->read_length);
		int covered_length;
		double dcov=coverage.calcRegionCoverage(bp,covered_length);
		vcoverage.push_back(dcov);
		vcovered_length.push_back(covered_length);

		//load in reads
		int vbamsize=bp.bam_aln_records.size();//number of reads 
		//cout<<"No. of reads is "<< vbamsize<<endl;///////////////////////////////////////////////////////
		for(int i=0;i<vbamsize;i++)
		{	
			
			Alignment alnmt;
			alnmt.setBar(bp.bam_aln_records[i]);
			int refid=bp.bam_aln_records[i]->rID;//start from 0
			int mate_refid=bp.bam_aln_records[i]->rNextID;
			int pos=bp.bam_aln_records[i]->pos;
			int mpos=bp.bam_aln_records[i]->pNext;
			
			int lcontig_len=vchroms[refid].second;
			int rcontig_len=vchroms[mate_refid].second;
			double outer_is=MEAN_INSERT_SIZE + 3*SD_INSERT_SIZE;//outer insert size
			
			//check whether both are aligned
			int pe_type=alnmt.getPEReadType();
			if(pe_type != READ_PAIR_MAP_TYPE_11) continue; //make sure both reads are aligned.
			else
			{
				//need to filter out those,  first in pair (also mate in other contig B), but insert size in smaller than contig length.
				double l_inner_dist=(double)(lcontig_len-pos);
				double r_inner_dist=(double)mpos;
				double max_allowed_pdist=outer_is-READ_LENGTH;//maximum allowed paired read distance, if larger, than see as unqualified pair
				if((l_inner_dist > max_allowed_pdist) || (r_inner_dist > max_allowed_pdist) ) continue;

				cnt_full_map_pairs++;
			}	
			
			string rname = faipsr.vchroms[refid].cname;//faipsr.getChromName(refid);
			string mrname= faipsr.vchroms[mate_refid].cname;//faipsr.getChromName(mate_refid);
			
			if(rname!=mrname)
			{	
				
				Contigs lctig, rctig;
				if(alnmt.isFirstInPair())
				{//first in pair 
					lctig.setContigName(rname);
					rctig.setContigName(mrname);
					lctig.setContigID(refid);
					rctig.setContigID(mate_refid);

					if(alnmt.isReverseStrand())
						lctig.setOrientation(ISREVERSE);
					else
						lctig.setOrientation(ISFORWARD);

					if(alnmt.isMateReverseStrand())
						rctig.setOrientation(ISREVERSE);
					else
						rctig.setOrientation(ISFORWARD);

					//calc contig distance according to insert size 
					double cp_dist;
					this->calcContigDistance(cp_dist,pos,mpos,vchroms[refid].second, vchroms[mate_refid].second);
					vcontig_pairs.push_back(std::make_pair(std::make_pair(lctig,rctig),cp_dist));			
				}
			}
			
		}//end of for 
		ctg_id++; 
		//output contig connections related to current contig 
		getUniqueContigPairs(fout_connect, vcontig_pairs, vchroms);
		
		std::vector<std::pair<std::pair<Contigs,Contigs>,double> >().swap(vcontig_pairs); //release vcontig_pairs 
		
	}//end of for
	

	bp.closeReader();//close BamReader 
	fout_connect.close();

	filterByCoverage(fname_temp, vcoverage, fcontig_info);

	//output contig information 
	ofstream fout_contig_info;
	string fctg_cov=fcontig_info+"_after_filter_cov_info.txt";
	fout_contig_info.open(fctg_cov.c_str());
	for(int i=0;i<vchroms.size();i++)
	{
		fout_contig_info<<vchroms[i].first<<" "<<vchroms[i].second<<" "<<vcoverage[i]<<" "<<vcovered_length[i]<<endl;
	}
	fout_contig_info.close();
}

void Refiner::calcCoverage(std::string fref, std::string fbam, std::string sfcov)
{
	//read in bam file 
	BamParse bp(fbam);
    bp.loadIndex();
	bp.openReader();

	//read in fai file 
	string fref_fai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, itemp1, itemp2;
	long long  ctg_start;
	int istart=0, iend;
	int ctg_id=0;

	int cnt_full_map_pairs=0;
	FaiParser faipsr(fref_fai);
	faipsr.parseFai();

	vector<pair<string,int> > vchroms;
	vchroms.clear();
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp1>>itemp2)
	{//for each contig 
		vchroms.push_back(std::make_pair(sctg,ctg_lenth));
	}
	fin_fai.close();

	//cout<<"Number of contigs is: "<<vchroms.size()<<endl;////////////////////////////////////////////////////////////////////////////////////

	vector<double> vcoverage;
	vector<int> vcovered_length;
	vcoverage.clear();
	vcovered_length.clear();
	for(int k=0;k<vchroms.size();k++)
	{//for each contig 
		sctg=vchroms[k].first;
		ctg_lenth=vchroms[k].second;

		bp.clearAll();//clear all the records
		bool bsignal=bp.dumpAlignments(ctg_id,0,ctg_id,ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}
		//cout<<"dumped # of reads is: "<<bp.bam_aln_records.size()<<endl;////////////////////////////////////////////////////////////////////////

		Coverage coverage;
		coverage.setRegionLenth(ctg_lenth);
		coverage.setReadLenth(this->read_length);
		int covered_len=0;
		double dcov=coverage.calcRegionCoverage(bp,covered_len);
		vcoverage.push_back(dcov);
		vcovered_length.push_back(covered_len);
		ctg_id++; 
	}
	bp.closeReader();//close BamReader 

	//output contig information  
	ofstream fout_contig_info;
	string fctg_cov=sfcov+"_before_filter_cov_info.txt";
	fout_contig_info.open(fctg_cov.c_str());
	for(int i=0;i<vchroms.size();i++)
	{
		fout_contig_info<<vchroms[i].first<<" "<<vchroms[i].second<<" "<<vcoverage[i]<<" "<<vcovered_length[i]<<endl;
	}
	fout_contig_info.close();
	
}

//calc
void Refiner::calcCoveageWithCutoff(std::string fref, std::string fbam, double read_cutoff, std::string sfcov)
{
	//read in bam file
	BamParse bp(fbam);
    bp.loadIndex();
	bp.openReader();

	//read in fai file 
	string fref_fai=fref+".fai";
	ifstream fin_fai;
	fin_fai.open(fref_fai.c_str());

	//count reads contig by contig 
	string sctg;
	int ctg_lenth, itemp1, itemp2;
	long long  ctg_start;
	int istart=0, iend;
	int ctg_id=0;

	int cnt_full_map_pairs=0;
	FaiParser faipsr(fref_fai);
	faipsr.parseFai();

	vector<pair<string,int> > vchroms;
	vchroms.clear();
	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp1>>itemp2)
	{//for each contig 
		vchroms.push_back(std::make_pair(sctg,ctg_lenth));
	}
	fin_fai.close();

	vector<double> vcoverage;
	vcoverage.clear();
	for(int k=0;k<vchroms.size();k++)
	{//for each contig 
		sctg=vchroms[k].first;
		ctg_lenth=vchroms[k].second;

		bp.clearAll();//clear all the records
		bp.dumpAlignments(ctg_id, 0, ctg_id, ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		//bool bsignal=bp.dumpAlignments(ctg_id,0,ctg_id,ctg_lenth);//read in region [start, end] of chromosome id (start from 0)
		/*if(bsignal==false)
		{
			std::cout<<"Cannot parse bam file "<<this->fbam<<std::endl;
			return;
		}*/

		Coverage coverage;
		coverage.setRegionLenth(ctg_lenth);
		coverage.setReadLenth(this->read_length);
		double dcov=coverage.calcRegionCoverageWithCutOff(bp,read_cutoff); //if the ratio of mapped part is smaller than cutoff, then read is filted out
		vcoverage.push_back(dcov);

		ctg_id++; 
	}
	bp.closeReader();//close BamReader

	ofstream fout_contig_info;
	string fctg_cov=sfcov+"_cov_info_with_cutoff.txt";
	fout_contig_info.open(fctg_cov.c_str());
	double total_bases=0;
	for(int i=0;i<vchroms.size();i++)
	{
		fout_contig_info<<vchroms[i].first<<" "<<vchroms[i].second<<" "<<vcoverage[i]<<endl;
		total_bases+=(double)vchroms[i].second * (double)vcoverage[i];
	}
	fout_contig_info<<"Total bases is: "<<total_bases<<endl;
	fout_contig_info.close();

	
}

//private functions---------------------------------------------------------------------------------------------------------------------
/*
Output:
	double dist, the distance between two contigs.
*/
void Refiner::calcContigDistance(double& gap, int pos, int mpos, int lcontig_len, int rcontig_len)
{
	int ldist=lcontig_len-pos;
	int rdist=mpos+READ_LENGTH; //not precisely maybe clipped 
	gap=MEAN_INSERT_SIZE-(double)ldist-(double)rdist;
}


void Refiner::getUniqueContigPairs(ofstream& fout, std::vector<std::pair<std::pair<Contigs,Contigs>,double> >& vcontig_pairs, vector<pair<string,int> >& vchroms)
{
	sort(vcontig_pairs.begin(),vcontig_pairs.end(),cmpContigs);
	
	int cnt_same=1;
	string llast="", rlast="";
	int lid_last, rid_last;
	string sl, sr;
	int ldir_last=-1, rdir_last=-1;
	double dist_last=0.0;
	int vcp_size=vcontig_pairs.size();
	double min_dist=MEAN_INSERT_SIZE, max_dist=-100000000000, total_dist=0;

	for(int i=0;i<vcp_size;i++)
	{
		
		string lcur=vcontig_pairs[i].first.first.getContigName();
		string rcur=vcontig_pairs[i].first.second.getContigName();
		int ldir_cur=vcontig_pairs[i].first.first.getOrientation();
		int rdir_cur=vcontig_pairs[i].first.second.getOrientation();
		int lid_cur=vcontig_pairs[i].first.first.getContigID();
		int rid_cur=vcontig_pairs[i].first.second.getContigID();
		//int llen_cur=vcontig_pairs[i].first.first.. 
		double dist_cur=vcontig_pairs[i].second;
		total_dist+=dist_cur;
		if(dist_cur>max_dist) max_dist=dist_cur;
		if(dist_cur<min_dist) min_dist=dist_cur;

		if(llast=="" || rlast=="")
		{
			llast = lcur;
			rlast = rcur;
			ldir_last=ldir_cur;
			rdir_last=rdir_cur;
			lid_last=lid_cur;
			rid_last=rid_cur;
		}
		else
		{
			if(lcur==llast && rcur==rlast && ldir_cur==ldir_last && rdir_cur==rdir_last)
				cnt_same++;
			else
			{
				if(ldir_last==ISFORWARD) sl="+"; else sl="-";	
				if(rdir_last==ISREVERSE) sr="+"; else sr="-";

				if(cnt_same > this->threshold)
					fout<<lid_last<<" "<<llast<<" "<<vchroms[lid_last].second<<" "<<sl<<" "<<rid_last<<" "<<rlast<<" "<<vchroms[rid_last].second<<" "<<sr<<" "<<cnt_same \
					<<" "<<min_dist<<" "<<max_dist<<" "<<total_dist/(double)cnt_same<<endl;

				min_dist=MEAN_INSERT_SIZE; max_dist=-100000000000; total_dist=0;
				cnt_same=1;
				llast=lcur;
				rlast=rcur;
				ldir_last=ldir_cur;
				rdir_last=rdir_cur;
				lid_last=lid_cur;
				rid_last=rid_cur;
			}
		}
	}
	
	//output the last one 
	if(ldir_last==ISFORWARD) sl="+"; else sl="-";			
	if(rdir_last==ISREVERSE) sr="+"; else sr="-";
	if(cnt_same > this->threshold )
		fout<<lid_last<<" "<<llast<<" "<<vchroms[lid_last].second<<" "<<sl<<" "<<rid_last<<" "<<rlast<<" "<<vchroms[rid_last].second<<" "<<sr<<" "<<cnt_same \
					<<" "<<min_dist<<" "<<max_dist<<" "<<total_dist/(double)cnt_same<<endl;

}

void Refiner::filterByCoverage(std::string fname, std::vector<double>& vcoverage, std::string fnew)
{
	ifstream fin;
	fin.open(fname.c_str());
	int lcontig_id, rcontig_id;
	int lcontig_len, rcontig_len;
	string lcontig, rcontig, lcontig_ori, rcontig_ori; 
	int cnt_pairs;
	double max_dist, min_dist, mean_dist;

	ofstream fout;
	fout.open(fnew.c_str());
	while(fin>>lcontig_id>>lcontig>>lcontig_len>>lcontig_ori>>rcontig_id>>rcontig>>rcontig_len>>rcontig_ori>>cnt_pairs>>min_dist>>max_dist>>mean_dist)
	{
		double large_cov = vcoverage[lcontig_id] > vcoverage[rcontig_id] ? vcoverage[lcontig_id] : vcoverage[rcontig_id];
		double small_cov = vcoverage[lcontig_id] > vcoverage[rcontig_id] ? vcoverage[rcontig_id] : vcoverage[lcontig_id];
		
		if((large_cov>0.0) && (((large_cov-small_cov)/large_cov) <= this->cf_cutoff))
		{
			fout<<lcontig_id<<" "<<lcontig<<" "<<lcontig_len<<" "<<lcontig_ori<<" "<<rcontig_id<<" "<<rcontig<<" "<<rcontig_len<<" "<<rcontig_ori<<" "<<cnt_pairs \
				<<" "<<min_dist<<" "<<max_dist<<" "<<mean_dist<<endl;
		}
	}
	fout.close();
	fin.close();
}

//void Refiner::cntContigLinkage_test(std::string fbam, std::string fref)
//{
//	//read in bam file
//	//BamParse bp(fbam);
//    //bp.loadIndex();
//
//	
//	//read in fai file
//	string fref_fai=fref+".fai";
//	ifstream fin_fai;
//	fin_fai.open(fref_fai.c_str());
//
//	//count reads contig by contig 
//	string sctg;
//	int ctg_lenth, ctg_start, itemp;
//	int istart=0, iend;
//	int ctg_id=0;
//
//	vector<pair<Contigs,Contigs> > vcontig_pairs;
//	vcontig_pairs.clear();
//
//	int cnt_full_map_pairs=0;
//	FaiParser faipsr(fref_fai);
//	faipsr.parseFai();
//
//	vector<pair<string,int> > vchroms;
//	vchroms.clear();
//	while(fin_fai>>sctg>>ctg_lenth>>ctg_start>>itemp>>itemp)
//	{//for each contig 
//		vchroms.push_back(std::make_pair(sctg,ctg_lenth));
//	}
//	fin_fai.close();
//
//	cout<<"No. of contigs is: "<<vchroms.size()<<endl; 
//
//	string fbam_index=fbam+".bai";
//	BamReader reader;
//	bool bsignal = reader.OpenIndex(fbam_index);
//	if ( !reader.Open(fbam) ) {
//		cerr << "Bamtools ERROR: could not open input BAM file: " << fbam << endl;
//		return ;
//	}
//
//	for(int i=0;i<501;i++)
//	{
//		int temp_len;
//		temp_len=vchroms[i].second;
//
//		BamRegion br(0,0,i,temp_len);
//		bool is_set=reader.SetRegion(br);
//	
//		int read_cnt=0;
//		BamAlignment al;
//		while ( reader.GetNextAlignment(al) )
//		{
//			if(al.Position<0) continue;
//		
//			read_cnt++;
//		}
//		cout<<"No. of reads is "<<read_cnt<<endl;
//		
//	}
//	reader.Close();
//}
