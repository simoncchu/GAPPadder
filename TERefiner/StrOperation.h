#ifndef _H_NOITAREPORTS_
#define _H_NOITAREPORTS_
#include<string>

class StrOperation
{
public:
	StrOperation(){}

public:
	void getReverseSupplementary(std::string& sin, std::string& sout);

private:
	char getSupplementary(char ci);
};

#endif