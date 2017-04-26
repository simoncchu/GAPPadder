#include<iostream>
#include<fstream>
#include<string>

#include"./algorithms/local_alignment.h"
#include"./algorithms/stdaln.h"
#include"scaffolding.h"
#include"fai_parser.h"
#include"public_func.h"

using namespace std;

const int OVERLAP=7;

Scaffolding::Scaffolding()
{

}


void Scaffolding::filterByDistance()
{
	
}

void Scaffolding::scaffold(std::string sfcontig, std::string sfconnect, std::string sfout)
{
	loadContigs(sfcontig);
	mergeContigs(sfconnect);
	std::string sfcontig_cnct_merge=sfconnect+"_merged.txt";
	constructConnectedContigs(sfcontig_cnct_merge, sfout);
}

//check contigs connections, and combine overlapped 
void Scaffolding::mergeContigs(std::string fcontig_connection)
{
	ifstream fin_connct;
	fin_connct.open(fcontig_connection.c_str());

	int lcontig_id, rcontig_id;
	int lcontig_len, rcontig_len;
	string lcontig, rcontig, lcontig_ori, rcontig_ori; 
	int cnt_pairs;
	double max_dist, min_dist, mean_dist;
	int overlap=0;

	std::string sfcontig_cnct_merge=fcontig_connection+"_merged.txt";
	std::string sfcontig_discard=fcontig_connection+"_discarded.txt";
	ofstream fout_merged;
	fout_merged.open(sfcontig_cnct_merge.c_str());
	ofstream fout_discarded;
	fout_discarded.open(sfcontig_discard.c_str());

	//string scontig_id;
	while(fin_connct>>lcontig_id>>lcontig>>lcontig_len>>lcontig_ori>>rcontig_id>>rcontig>>rcontig_len>>rcontig_ori>>cnt_pairs>>min_dist\
		>>max_dist>>mean_dist)
	{	
		overlap=0;
		bool bmerged=false;
		//cout<<lcontig_id<<" "<<lcontig<<" "<<rcontig_id<<" "<<rcontig<<endl;//////////////////////////////////////////////////////////////////////////////////////////////
		
		if(mean_dist>0) 
		{
			//output 
			fout_merged<<lcontig_id<<" "<<lcontig<<" "<<lcontig_len<<" "<<lcontig_ori<<" "<<rcontig_id<<" "<<rcontig<<" "<<rcontig_len\
				       <<" "<<rcontig_ori<<" "<<cnt_pairs<<" "<<min_dist\
					   <<" "<<max_dist<<" "<<mean_dist<<" "<<bmerged<<" "<<overlap<<endl;
			continue; //only check those whose gap is negative 
		}

		string lc,rc;
		if(lcontig_ori=="+")
			lc=mcontigs[lcontig].first;
		else 
			lc=mcontigs[lcontig].second;

		if(rcontig_ori=="+")
			rc=mcontigs[rcontig].first;
		else
			rc=mcontigs[rcontig].second;
		
		int min_contig;
		if(lcontig_len>rcontig_len) min_contig=rcontig_len;
		else min_contig=lcontig_len;
		int min_gap=(int)(-1*min_dist);
		if(min_contig<min_gap)
			min_gap=min_contig;
		
		std::string prefix_lsubstr;
		if(min_gap>=lcontig_len)
			prefix_lsubstr="";
		else
			prefix_lsubstr=lc.substr(0,lcontig_len-min_gap);
		std::string lsubstr=lc.substr(lcontig_len-min_gap,min_gap);
		
		std::string rsubstr=rc.substr(0,min_gap);
		std::string suffix_rsubstr;
		if(min_gap>=rcontig_len)
			std::string suffix_rsubstr="";
		else
			suffix_rsubstr=rc.substr(min_gap,rcontig_len-min_gap);

		LocalAlignment la;
		int loptm_start=-1,loptm_end=-1,roptm_start=-1, roptm_end=-1;		
		la.optAlign(lsubstr, rsubstr, loptm_start, loptm_end, roptm_start, roptm_end); //1-based positions 
		
		//check whether they have overlap ....
		overlap=loptm_end-loptm_start+1;
		//also check whether prefix and suffix 
		if((overlap>OVERLAP) && (loptm_end==lsubstr.length()) && (roptm_start==1))
		{			
			//merge 
			bmerged=true;
			std::string smerged=lc+ rsubstr.substr(overlap,(min_gap-overlap)) + suffix_rsubstr;
			string newid=lcontig+lcontig_ori+rcontig+rcontig_ori;
			mmerged_nodes[newid]=overlap;
			//output 
			fout_merged<<lcontig_id<<" "<<lcontig<<" "<<lcontig_len<<" "<<lcontig_ori<<" "<<rcontig_id<<" "<<rcontig<<" "<<rcontig_len\
				       <<" "<<rcontig_ori<<" "<<cnt_pairs<<" "<<min_dist\
					   <<" "<<max_dist<<" "<<mean_dist<<" "<<bmerged<<" "<<overlap<<endl;
		}
		else
		{//they have negative distance, but no apparent overlap, then discard the connection
			fout_discarded<<lcontig_id<<" "<<lcontig<<" "<<lcontig_len<<" "<<lcontig_ori<<" "<<rcontig_id<<" "<<rcontig<<" "<<rcontig_len\
				       <<" "<<rcontig_ori<<" "<<cnt_pairs<<" "<<min_dist\
					   <<" "<<max_dist<<" "<<mean_dist<<endl;
		}
	}
	fin_connct.close();
	fout_merged.close();
	fout_discarded.close();
	
}

/*
Description:
	Load the sequences into the map<id,<seq, reverse_supplementary_seq> > 
*/
void Scaffolding::loadContigs(std::string sfcontigs)
{
	ifstream fin_contig;
	fin_contig.open(sfcontigs.c_str());

	std::string fa_id, fa_seq, fa_rc_seq;
	string sline="";
	bool bfirst=true;
	while(std::getline(fin_contig,sline))
	{
		if(sline.length()>0 && sline[0]=='>')
		{
			if(bfirst==true)
			{
				bfirst=false;
				fa_id=sline.substr(1,sline.length()-1);
				continue;
			}
			getReverseSupplementary(fa_seq, fa_rc_seq);
			mcontigs[fa_id]=std::make_pair(fa_seq,fa_rc_seq);
			fa_id=sline.substr(1,sline.length()-1);
			fa_seq="";
		}
		else
		{
			fa_seq+=sline;
		}
	}
	mcontigs[fa_id]=std::make_pair(fa_seq,fa_rc_seq);//the last one 

	//cout<<"Load contigs.fa finished..."<<endl;
	fin_contig.close();
}

//construct the connected contigs.
void Scaffolding::constructConnectedContigs(std::string sfcontig_cnct, std::string sfout)
{
	//build the graph-------------------------------------------------------------
	int lcontig_id, rcontig_id;
	int lcontig_len, rcontig_len;
	string lcontig, rcontig, lcontig_ori, rcontig_ori; 
	int cnt_pairs;
	double max_dist, min_dist, mean_dist;
	bool bmerged;
	int overlap;

	ifstream fin_connct;
	fin_connct.open(sfcontig_cnct.c_str());

	map<string, int> mgraph;
	mgraph.clear();
	int index_graph=0;
	vector<ContigNode*> vpcn;//vector of contig node pointer 
	while(fin_connct>>lcontig_id>>lcontig>>lcontig_len>>lcontig_ori>>rcontig_id>>rcontig>>rcontig_len\
				   >>rcontig_ori>>cnt_pairs>>min_dist>>max_dist>>mean_dist>>bmerged>>overlap)
	{
		string sunique=lcontig+lcontig_ori;
		map<string,int>::iterator it=mgraph.find(sunique);
		int itemp=-1;
		if(it==mgraph.end())
		{//not found, creat a new one. 
			ContigNode* pnew=new ContigNode();
			pnew->pmate=NULL;
			if(lcontig_ori=="+")
				pnew->brc=true;
			else 
				pnew->brc=false;
			pnew->support_PE=-1;
			pnew->sname=lcontig;
			pnew->dist=mean_dist;
			pnew->overlap=overlap;
			vpcn.push_back(pnew);
			mgraph[sunique]=index_graph;
			itemp=index_graph;
			index_graph++;
		}
		else
		{//already exist 
			itemp=it->second;
		}

		ContigNode* pnode=new ContigNode();
		pnode->sname=rcontig;
		if(rcontig_ori=="+")
			pnode->brc=true;
		else
			pnode->brc=false;
		pnode->support_PE=cnt_pairs;
		pnode->dist=mean_dist;
		pnode->overlap=overlap;
		pnode->pmate=vpcn[itemp]->pmate;
		vpcn[itemp]->pmate=pnode;
	}

	//traverse the graph and output the connected contigs-------------------------------------------
	//For now, only output the two-connected situations---------------------------------------------
	ofstream fcnct_contigs;
	fcnct_contigs.open(sfout.c_str());
	int vsize=vpcn.size();
	for(int i=0;i<vsize;i++)
	{
		int itotal_pe=0;
		int itotal_nodes=0;
		ContigNode* pcur=vpcn[i];
		while(pcur!=NULL)
		{
			if(pcur->support_PE==-1) 
			{
				pcur=pcur->pmate;
				continue;
			}
			else
			{
				itotal_pe+=pcur->support_PE;
				itotal_nodes++;
				pcur=pcur->pmate;
			}
			
		}
		if(itotal_nodes==0)
		{
			cout<<"The constructed graph is wrong!!!!"<<endl;
			break;
		}
		int ave_pe=itotal_pe/itotal_nodes;
		
		//output all those supported PE larger or equal to the average one
		pcur=vpcn[i];
		while(pcur!=NULL)
		{		
			if(pcur->support_PE >= ave_pe)
			{
				string temp_id="";
				string temp_ori1="+", temp_ori2="+";
				if(vpcn[i]->brc==false) temp_ori1="-"; 
				if(pcur->brc==false) temp_ori2="-";
				string sdis="";
				if(pcur->overlap==0)
					sdis=PubFuncs::cvtInt2Str(int(pcur->dist));
				else
					sdis=PubFuncs::cvtInt2Str(-1*pcur->overlap);

				temp_id=">"+vpcn[i]->sname+"$" + temp_ori1 + "$" + pcur->sname + "$" + temp_ori2+"$"+sdis;
				

				string lseq, rseq;
				if(vpcn[i]->brc==true)
					lseq=mcontigs[vpcn[i]->sname].first;
				else
					lseq=mcontigs[vpcn[i]->sname].second;

				if(pcur->brc==true)
					rseq=mcontigs[pcur->sname].first;
				else
					rseq=mcontigs[pcur->sname].second;
				

				//check is overlap or not 
				string final_seq="";
				if(pcur->overlap==0)
				{//non-overlap
					//check distance 
					string seperator="";
					for(int j=0;j<(int)pcur->dist;j++)
						seperator+="N";

					final_seq=lseq+seperator+rseq;
				}
				else
				{//overlap
					if(pcur->overlap>0)
						final_seq=lseq+rseq.substr(pcur->overlap, rseq.length() - pcur->overlap);	
					else
						cout<<"Wrong overlap value, it's negative!!!"<<endl;
				}

				string sfa_fmt=""; 
				this->cvtStr2Fa(final_seq,60,sfa_fmt);//output the connected contigs
				fcnct_contigs<<temp_id<<endl;
				fcnct_contigs<<sfa_fmt;
			}	
			pcur=pcur->pmate;
		}//end of while
	}//end of for
	fcnct_contigs.close();

	//release the nodes------------------------------------------------------------------------------
	ContigNode* prls;
	for(int i=0;i<vsize;i++)
	{
		prls=vpcn[i];
		while(prls!=NULL)
		{
			ContigNode* pcur=prls;
			pcur=prls->pmate;
			delete prls;
			prls=pcur;
		}
	}
	
	fin_connct.close();
}


void Scaffolding::cvtStr2Fa(std::string& seq, int size_aline, std::string& sout)
{
	sout="";
	int start=0;
	while(start<seq.length())
	{
		sout+=seq.substr(start,size_aline);
		sout+="\n";
		start+=size_aline;
	}
}

/*
Get reverse supplementary. 
*/
void Scaffolding::getReverseSupplementary(std::string& sin, std::string& sout)
{
	sout="";
	int slen=sin.length();
	for(int i=slen-1;i>=0;i--)
	{
		sout+=getSupplementary(sin[i]);
	}
}

char Scaffolding::getSupplementary(char ci)
{
	if(ci=='A' || ci=='a')
		return 'T';
	else if(ci=='T' || ci=='t')
		return 'A';
	else if(ci=='C' || ci=='c')
		return 'G';
	else if(ci=='G' || ci=='g')
		return 'C';
	else
		return 'N';
}