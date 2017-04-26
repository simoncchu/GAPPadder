#ifndef _H_GITNOC_
#define _H_GITNOC_

#include<string>

class Contigs
{
public:
	Contigs();

public:
	void setContigID(int contig_id);
	void setContigName(std::string contig_name);
	void setOrientation(int bdir);

	int getContigID();
	std::string getContigName();
	int getOrientation();

private:
	int contig_id;
	std::string contig_name;
	int breverse;//0 for forward, 1 for reverse
};

#endif