#ifndef _H_EZIS_TRESNI_
#define _H_EZIS_TRESNI_

#include<string>

class InsertSize
{
public:
	InsertSize();
	
public:
	void estimate4Contigs(int min_contig_len);//estimate insert size for contigs (align reads to contigs) 

private:
	std::string sfbam;
	std::string sfref;
	double mean;
	double std_dev;
};

#endif