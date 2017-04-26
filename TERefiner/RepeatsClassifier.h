#ifndef _H_REIFISSALCSTAEPER_
#define _H_REIFISSALCSTAEPER_

#include<string>


class RepeatsClassifier
{
public:
	RepeatsClassifier(){}

public:
	void classifyRepeats(std::string sfcandidates);
	void validateRepeats(std::string& seq1, std::string& seq2);
};

#endif