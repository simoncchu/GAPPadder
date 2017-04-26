#ifndef _H_RENIFER_
#define _H_RENIFER_

#include<string>
#include<vector>
#include<utility>
#include<fstream>
#include"contigs.h"

class Refiner
{
public:
	Refiner(){};
	Refiner(std::string fbam, std::string ffai, std::string fout_folder);

public:
	void setCutOff(double cf_cutoff);
	void setContigFile(std::string fcontig);
	void setThreshold(int threshold);
	void setReadLength(int readlen);

public:
	void refineByReadsCombinedContigs();
	void refineByReads(std::string fbam, std::string fref, std::string fout);

	void removeRepeatsOfTwoContigSets(std::string fbam, std::string fref, std::string bam_fasta, std::string fnew_contig_fa);
	void removeRepeatsOfOneContigSet(std::string fbam, std::string fref, std::string fnew);
	void removeDupRepeatsOfOneContigSet(std::string fbam, std::string fref, std::string fnew, double cutoff_ratio, bool brm_cntn);//only remove the duplicates 
	void removeContainedContigs(std::string fbam, std::string fref, std::string fnew);
	void cntContigLinkage(std::string fbam, std::string fref, std::string fcontig_info);//cnt number of paired-end reads that link each two contigs
	void cntContigLinkage_test(std::string fbam, std::string fref);
	void gnrtUniqueFa(std::string fa_old, std::string fa_new);

	void calcCoverage(std::string fref, std::string fbam, std::string sfcov);
	void calcCoveageWithCutoff(std::string fref, std::string fbam, double cutoff, std::string sfcov);

	void evaluateWithBenchmark(std::string fbam, std::string fref, double cutoff_map_fraction);
	
	void filteMergedContigsByClip(std::string fbam, std::string fref, double cutoff_clip_fraction, double coverage);


private:
	void calcContigDistance(double& dist, int pos, int mpos, int lcontig_len, int rcontig_len);
	void getUniqueContigPairs(std::ofstream& fout, std::vector<std::pair<std::pair<Contigs,Contigs>,double> >& vcontig_pairs, \
		std::vector<std::pair<std::string,int> >& vchroms);
	void filterByCoverage(std::string fname, std::vector<double>& vcoverage, std::string fnew);
	void rmCotigs(std::string foriginal, std::vector<int>& vid_rm, std::string foutput);
	void getLongestUncover(std::vector<std::pair<int,int> >& vsgmts, int& luncover);

private:
	std::string fcontig;//contig file 
	double cf_cutoff;
	int threshold;
	int read_length;
	std::string fbam;//bam file path
	std::string ffai;//fai file path
	std::string fout_folder;//output folder
};

#endif