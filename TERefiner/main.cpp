/*To do list:
1. classify the tandem repeats and sparsed repeats.
2. check the overlap, if the distance is negative. And concatenate if necessary. 
3. Many TE share same contig, like those L1TA series. So greedy doesn't work for this.
4. Is the coverage is reliable? Because some contigs are short????????????????????????????????????
5. When remove repeats, should have some flexible, not strictly fully map.
6. Error: the output contigs.fa file have one empty line output. Should be the problem of RemoveRepeatsOfTwoContigSets???????????????????? 
*/

#include<iostream>
#include<string>
#include<stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include"refiner.h"
#include"RepeatsClassifier.h"
#include"scaffolding.h"
#include"public_parameters.h"
#include"./algorithms/local_alignment.h"

using namespace std;

int READ_LENGTH;//extern value, declared in public_parameters.h 
double MEAN_INSERT_SIZE;//mean insert size 
double SD_INSERT_SIZE;//standard derivation of insert size 
double READ_FULL_MAPPED_CUTOFF;//larger than this cutoff, than a read will be considered as fully mapped

int main(int argc, char* argv[])
{
	bool ba=false;
	bool bm=false;
	bool bc=false;
	bool bt=false;
	bool bo=false;
	bool bl=false;
	bool bs=false;
	bool bu=false;
	bool bg=false;
	bool bb=false;
	bool bp=false;
	bool be=false;
	bool bk = false;
	bool flag_l=false, flag_c=false, flag_t=false, flag_s=false, flag_b=false, flag_r=false, flag_q=false, flag_o=false, \
		 flag_m=false, flag_d=false, flag_v=false, flag_g=false;
	double cutoff_ratio=0.85;
	//double is_mean, is_derivation;//mean insert size, and standard derivation.
	int threshold=5;
	int read_length=106;
	bool brm_contained=false;
	READ_LENGTH=106;
	SD_INSERT_SIZE=50;
	READ_FULL_MAPPED_CUTOFF=0.85;
	string fbam, fref, fref2, ffastq, fout, fbam_cov;

	int c;
	while ((c = getopt(argc, argv, "AMCTOLSUGBPKEl:c:t:s:b:r:q:o:m:d:v:g")) != -1)
	switch (c)
    {
	case 'A'://
		ba=true;
		break;
	case 'M'://local aligment of two string
		bm=true;
		break;
	case 'C': //remove according to ratio of fully mapped reads
        bc=true;
		break;
    case 'T': //Remove repeats in two contig sets
		bt=true;
		break;
	case 'O'://remove repeats in one congit set
		bo=true;
		break;
	case 'P':
		bp=true;
		break;
	case 'K'://remove contained ones
		bk = true;
		break;
	case 'L'://contigs linkages
		bl=true;
		break;
	case 'S'://Scaffolding 
		bs=true;
		break;
	case 'U'://generate unique fa 
		bu=true;
		break;
	case 'G'://calc coverage of contigs with given bam and cutoff 
		bg=true;
		break;
	case 'B': //calc coverage of contigs with covered length 
		bb=true;
		break;
	case 'E'://evaluate ratio of contigs cover benchmark
		be=true;
		break;
	case 'l'://read length
		READ_LENGTH=atoi(optarg);
		flag_l=true;
		break;
	case 'c'://cutoff ratio 
		cutoff_ratio=atof(optarg);
		flag_c=true;
		break;
	case 't':
		threshold=atoi(optarg);
		flag_t=true;
		break;
	case 's':
		fref2=(string)optarg;
		flag_s=true;
		break;
	case 'b'://bam file path
		fbam=(string)optarg;
		flag_b=true;
		break;
	case 'r'://reference file path
		fref=(string)optarg;
		flag_r=true;
		break;
	case 'q'://fastq file path 
		ffastq=(string)optarg;
		flag_q=true;
		break;
	case 'o'://output file path
		fout=(string)optarg;
		flag_o=true;
		break;
	case 'm'://mean insert size
		MEAN_INSERT_SIZE=atof(optarg);
		flag_m=true;
		break;
	case 'd'://standard derivation of insert size 
		SD_INSERT_SIZE=atof(optarg);
		flag_d=true;
		break;
	case 'v'://bam file for calc coverage
		fbam_cov=(string)optarg;
		flag_v=true;
		break;
	case 'g'://remove perfect contained
		flag_g=true;
		brm_contained=true;
		break;
	default:
		cout<<"Wrong parameters"<<endl;
	}

	if(bc && flag_c&& flag_b && flag_r&& flag_o && flag_l)
	{
		Refiner rfnr;
		rfnr.setCutOff(cutoff_ratio);
		rfnr.refineByReads(fbam,fref,fout);
	}
	else if(bt && flag_r && flag_b && flag_c && flag_o && flag_s)
	{//remove repeats of two sets of contigs
		Refiner rfnr;
		READ_FULL_MAPPED_CUTOFF=cutoff_ratio;
		rfnr.removeRepeatsOfTwoContigSets(fbam,fref,fref2,fout);
	}
	else if(bo && flag_b && flag_r && flag_c && flag_o)
	{//remove repeats of one set of contigs 
		Refiner rfnr;
		READ_FULL_MAPPED_CUTOFF=cutoff_ratio;
		rfnr.removeRepeatsOfOneContigSet(fbam,fref,fout);
	}
	else if(bp && flag_b && flag_r && flag_c && flag_o)
	{//remove duplicates(with or without contained ones)
		Refiner rfnr;
		READ_FULL_MAPPED_CUTOFF=cutoff_ratio;
		rfnr.removeDupRepeatsOfOneContigSet(fbam,fref,fout,cutoff_ratio,brm_contained);
	}
	else if (bk && flag_b && flag_r && flag_c && flag_o)
	{
		Refiner rfnr;
		READ_FULL_MAPPED_CUTOFF = cutoff_ratio;
		rfnr.removeContainedContigs(fbam, fref, fout);
	}
	else if(bl && flag_b && flag_r && flag_c && flag_t && flag_l && flag_o && flag_m && flag_d && flag_v)
	{
		Refiner rfnr;
		rfnr.setCutOff(cutoff_ratio);//coverage cutoff ratio
		rfnr.setThreshold(threshold);//minimum number of supported read pairs

		string sprefix=fbam.substr(0,fbam.size()-9);
		string fcov_bam = sprefix + ".sam_for_coverage.sorted.bam";
		cout<<fcov_bam<<endl;
		rfnr.calcCoveageWithCutoff(fref,fcov_bam,0,fout);
		rfnr.cntContigLinkage(fbam,fref,fout);
	}
	else if(bs && flag_r && flag_s && flag_o)
	{
		Scaffolding scfd;
		scfd.scaffold(fref,fref2,fout);
	}
	else if(bu && flag_r && flag_o)
	{
		Refiner rfnr;
		rfnr.gnrtUniqueFa(fref,fout);
	}
	else if(ba && flag_s && flag_r)
	{
		RepeatsClassifier rc;
		rc.validateRepeats(fref,fref2);
	}
	else if(bm &&  flag_s && flag_r)
	{
		LocalAlignment la;
		int opt_start_ref=-1, opt_end_ref=-1, opt_start_sgmt=-1, opt_end_sgmt=-1;
		la.optAlign(fref,fref2,opt_start_ref, opt_end_ref, opt_start_sgmt, opt_end_sgmt);
		cout<<opt_start_ref<<" "<<opt_end_ref<<" "<<opt_start_sgmt<<" "<<opt_end_sgmt<<endl;
	}
	else if(bg && flag_b && flag_r && flag_c && flag_o && flag_l)
	{
		Refiner rfnr;
		rfnr.calcCoveageWithCutoff(fref,fbam,cutoff_ratio,fout);
	}
	else if(bb && flag_b && flag_r && flag_o && flag_l)
	{
		Refiner rfnr;
		rfnr.calcCoverage(fref,fbam,fout);
	}
	else if(be && flag_b && flag_r && flag_c)
	{
		Refiner rfnr;
		rfnr.evaluateWithBenchmark(fbam, fref,cutoff_ratio);
	}
	else
	{
		cout<<"Please check parameters setting!"<<endl;
	}

	return 0;
}
