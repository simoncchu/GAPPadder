#ifndef _H_EGAREVOC_
#define _H_EGAREVOC_

#include"bam_parse.h"

class Coverage
{
public:
	Coverage();
	void setRegionLenth(int len);
	void setReadLenth(int read_lenth);

public:
	double calcRegionCoverage(BamParse& bp, int& covered_length);
	double calcRegionCoverageWithCutOff(BamParse& bp, double cutoff);//calc coverage with cutoff that used to filter out reads that only have small part is aligned

private:
	int region_lenth;
	int read_len;
};

#endif