#ifndef _H_NAMRETAW_HTIMS_
#define _H_NAMRETAW_HTIMS_

#include<string>
#include<vector>
#include<utility>

class SmithWaterman
{
public:
	SmithWaterman();
	SmithWaterman(std::string sref, std::string ssgmt);

public:
	void align();//align segment with reference sequence
	void outputCigar();

private:
	void loadMatrix(int** &matrix, int size_ref, int size_sgmt);//allocate memory for score_matrix and path_matrix
	void initMatrix(int size_ref, int size_sgmt);//initialize score matrix
	void freeMatrix(int** &matrix, int size_sgmt);//release memory of score_matrix and path_matrix
	void traceback(int max_row, int max_col);//trace back 

private:
	std::string sref;//reference sequence
	std::string ssgmt;//segment sequence
	int** score_matrix;//score mathrix
	int** path;//matrix save moving paths
	std::vector<std::pair<int,std::string> > cigar;//save cigar information
	int lm_map_pos;//left most mapping position
	int max_score;//maximum mapping socre
};

#endif

