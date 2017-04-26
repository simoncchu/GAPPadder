#ifndef _H_RESRAP_IAF_
#define _H_RESRAP_IAF_

#include<string>
#include<vector>
#include<utility>

class ChromInfo
{
public:
	int id;
	std::string cname;
	int length;
	unsigned long startpos;
	int size_chars;
	int size_ascii;
};

class FaiParser
{
public:
	FaiParser();
	FaiParser(std::string path);

public:
	void parseFai();
	int getChromLen(std::string chrom);
	void setPath(std::string path);

public:
	std::vector<ChromInfo> vchroms;//save all the chrom info in vectors. 

private:
	std::string path;
};

#endif

