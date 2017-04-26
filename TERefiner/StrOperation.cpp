#include"StrOperation.h"

/*
Get reverse supplementary. 
*/
void StrOperation::getReverseSupplementary(std::string& sin, std::string& sout)
{
	sout="";
	int slen=sin.length();
	for(int i=slen-1;i>=0;i--)
	{
		sout+=getSupplementary(sin[i]);
	}
}

char StrOperation::getSupplementary(char ci)
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