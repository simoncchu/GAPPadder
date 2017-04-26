#include<iostream>
#include<string>
#include<vector>
#include<fstream>

#include"RepeatsClassifier.h"
#include"public_func.h"
#include"StrOperation.h"
#include"./algorithms/local_alignment.h"

using namespace std;

void RepeatsClassifier::classifyRepeats(std::string sfcandidates)
{
	//ifstream fin_cand;
	//fin_cand.open(sfcandidates.c_str());

	//vector<vector<string> > vcandidate;
	//vcandidate.clear();
	//string sline;
	//while(std::getline(fin_cand,sline))
	//{
	//	vector<string> vgroups;
	//	vgroups=PubFuncs::split(sline,' ');
	//	int iround=vgroups.size();
	//	if(iround<=0) 
	//	{
	//		cout<<"No candiddate contigs to support tandame repeats!!!";
	//		continue;
	//	}
	//	for(int i=0;i<iround;i++)
	//	{
	//		vector<string> vcontig=PubFuncs::split(vgroups[i],',');
	//		vcandidate.push_back(vcontig);
	//	}

	//	
	//	
	//}

	//

	//fin_cand.close(); 
}


void RepeatsClassifier::validateRepeats(string& seq1, string& seq2)
{
	//align seq1 and seq2 
	LocalAlignment la;
	int optm_start_ref=-1, optm_end_ref=-1, optm_start_sgmt=-1, optm_end_sgmt=-1;
	int scnd_start_ref=-1, scnd_end_ref=-1, scnd_start_sgmt=-1, scnd_end_sgmt=-1;
	//la.optAlign(seq1, seq2, optm_start_ref, optm_end_ref, optm_start_sgmt, optm_end_sgmt);
	//In case of TR rotation, we still need an alignment for the rest part 
	//first concatenate the un-optimal parts 
	la.optAlignWithRestSecondOpt(seq1, seq2, optm_start_ref, optm_end_ref, optm_start_sgmt, optm_end_sgmt,\
		scnd_start_ref, scnd_end_ref, scnd_start_sgmt, scnd_end_sgmt);

	int optm_lenth=0;
	if(optm_end_ref>0)
		optm_lenth=optm_end_ref-optm_start_ref+1;
	else 
		optm_lenth=0;

	int scnd_lenth=0;
	if(scnd_end_ref>0)
		scnd_lenth=scnd_end_ref-scnd_start_ref+1;
	else
		scnd_lenth=0;

	int iff=optm_lenth + scnd_lenth;//forward forward

	//align seq1 and reverse supplementary of seq2 
	StrOperation so;
	string seq2_rs;
	so.getReverseSupplementary(seq2,seq2_rs);
	int optm_start_ref_rs=-1, optm_end_ref_rs=-1;
	int scnd_start_ref_rs=-1, scnd_end_ref_rs=-1;
	la.optAlignWithRestSecondOpt(seq1, seq2_rs, optm_start_ref_rs, optm_end_ref_rs, optm_start_sgmt, optm_end_sgmt,\
	scnd_start_ref_rs, scnd_end_ref_rs, scnd_start_sgmt, scnd_end_sgmt);

	int optm_rs_lenth=0;
	if(optm_end_ref_rs>0)
		optm_rs_lenth=optm_end_ref_rs-optm_start_ref_rs+1;
	else 
		optm_rs_lenth=0;

	int scnd_rs_lenth=0;
	if(scnd_end_ref_rs>0)
		scnd_rs_lenth=scnd_end_ref_rs-scnd_start_ref_rs+1;
	else
		scnd_rs_lenth=0;
	int ifr = optm_rs_lenth + scnd_rs_lenth; //forward reverse 

	//compare which one has a better results 
	if( iff > ifr)
		cout<<iff<<endl;
		//cout<<optm_start_ref<<" "<<optm_end_ref<<endl; 
	else
		cout<<ifr<<endl;
		//cout<<optm_start_ref_rs<<" "<<optm_end_ref_rs<<endl;
}

