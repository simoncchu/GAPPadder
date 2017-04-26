#ifndef _H_LOCAL_ALIGNMENT_
#define _H_LOCAL_ALIGNMENT_

#include<string>

class LocalAlignment
{
public:
	LocalAlignment(){};
	/*LocalAlignment(std::string sref, std::string ssgmt);*/

public:
	void optAlign(std::string& sref, std::string& ssgmt, int& optm_start_ref, int& optm_end_ref, int& optm_start_sgmt, int& optm_end_sgmt);
	void align(std::string& sref, std::string& ssgmt, int& optm_start_ref, int& optm_end_ref, int& optm_start_sgmt, int& optm_end_sgmt, \
	           int& start_ref1, int& end_ref1, int& start_sgmt1, int& end_sgmt1,                 \
			   int& start_ref2, int& end_ref2, int& start_sgmt2, int& end_sgmt2);

	void optAlignWithRestSecondOpt(std::string& sref, std::string& ssgmt, int& optm_start_ref, int& optm_end_ref, int& optm_start_sgmt, int& optm_end_sgmt,  \
	                       int& start_ref1, int& end_ref1, int& start_sgmt1, int& end_sgmt1);
//	void setRef(std::string sref);
//	void setSgmt(std::string ssgmt);
//
//private:
//	std::string sref;//reference sequence
//	std::string ssgmt;//segment sequence
};

#endif