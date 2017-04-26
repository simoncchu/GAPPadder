#ifndef _H_GNIDLOFFACS_
#define _H_GNIDLOFFACS_

#include<string>
#include<vector>
#include<map>

class ContigNode
{
public:
	ContigNode(){pmate=NULL; overlap=-1; support_PE=-1; dist=0;};
public:
	std::string sname;//contig_name
	bool brc;//is reverse complement
	double cov;//coverage
	int support_PE;//number of supported read pairs
	int overlap;
	double dist;
	ContigNode* pmate;//connected nodes
};


class Scaffolding
{
public:
	Scaffolding();
	void setContigFile(std::string sfcontigs);

public:
	void scaffold(std::string sfcontig, std::string sfconnect, std::string fout);
	void mergeContigs(std::string fcontig_connection);

	void loadContigs(std::string sfcontigs);
	void constructConnectedContigs(std::string sfcontig_cnct, std::string sfout);

private:
	void cvtStr2Fa(std::string& seq, int size_aline, std::string& sout);
	void filterByDistance();
	void getReverseSupplementary(std::string& sin, std::string& sout);
	char getSupplementary(char cin);
	
private:
	std::map<std::string,std::pair<std::string, std::string> > mcontigs; //<ID,<seq, reverse_supplementary> > 
	std::map<std::string, int> mmerged_nodes;//<ID_ori_ID_ori, overlap> 
	std::string sfbam;//bam file
	std::string sfcontigs;//contigs.fa file 
};

#endif


