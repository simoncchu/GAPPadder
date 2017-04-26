#include"smith-waterman.h"
#include<iostream>

const int START=0;
const int DIAG=1;
const int LEFT=2;
const int TOP=3;

const std::string SINSERTION="I";
const std::string SDELETION="D";
const std::string SMATCH="M";
const std::string SMISMATCH="U";
const std::string SSTART="S";

SmithWaterman::SmithWaterman()
{
	this->sref="";
	this->ssgmt="";
	cigar.clear();
	max_score=0;
	lm_map_pos=0;
}

SmithWaterman::SmithWaterman(std::string sref, std::string ssgmt)
{
	this->sref=sref;
	this->ssgmt=ssgmt;
	cigar.clear();
	max_score=0;
	lm_map_pos=0;//left most mapping position 
}


/*
Description:
	Align segment with reference.
Input:
	string: sref, reference.
	string: ssgmt, segment.
Output:
	int: pos, leftmost mapping position.
	string: cigar, mapping cigar.
*/
void SmithWaterman::align()
{
	int size_ref=sref.length();
	int size_sgmt=ssgmt.length();

	loadMatrix(this->score_matrix,size_ref,size_sgmt);
	initMatrix(size_ref,size_sgmt);

	int ipref=1;
	int ipsgmt=1;
	
	int max_col=0;
	int max_row=0;
	for(int i=1;i<=size_sgmt;i++)
	{
		for(int j=1;j<=size_ref;j++)
		{
			int diag,left,top;
			if(ssgmt[i-1]==sref[j-1])
			{
				diag=score_matrix[i-1][j-1]+2;
			}
			else
			{
				diag=score_matrix[i-1][j-1]-1;
			}
			left=score_matrix[i][j-1]-1;
			top=score_matrix[i-1][j]-1;

			if(diag<0 && left<0 && top<0)
			{
				score_matrix[i][j]=0;
				path[i][j]=START;
			}
			else if(top>left && top>diag)
			{
				score_matrix[i][j]=top;
				path[i][j]=TOP;
			}
			else if(left>diag && left>top)
			{
				score_matrix[i][j]=left;
				path[i][j]=LEFT;
			}
			else
			{
				score_matrix[i][j]=diag;
				path[i][j]=DIAG;
			}
			if(score_matrix[i][j]>max_score)
			{
				max_score=score_matrix[i][j];
				max_row=i;
				max_col=j;
			}
		}
	}
	
	//output socore_matrix/////////////////////////////////////////////////////////////////////
	//for(int i=1;i<=size_sgmt;i++)
	//{
	//	for(int j=1;j<=size_ref;j++)
	//	{
	//		std::cout<<path[i][j]<<" ";
	//	}
	//	std::cout<<std::endl;
	//}
	//std::cout<<max_row<<" "<<max_col<<std::endl; 



	traceback(max_row, max_col);

	freeMatrix(this->score_matrix,size_sgmt);
}

void SmithWaterman::outputCigar()
{
	int size=cigar.size();
	for(int i=0;i<size;i++)
	{
		std::cout<<cigar[i].first<<cigar[i].second;
	}
	std::cout<<std::endl;
}

//------------private functions---------------------------------------------------------------------------
/*
Description:
	trace back
Input:
	int: max_row, row of the highest score in score_matrix;
	int: max_col, column of the highest score in score_matrix.
Output:
	string: cigar, mapping information;
	int: pos, left-most mapping position.
*/
void SmithWaterman::traceback(int max_row, int max_col)
{
	int iprow=max_row;
	int ipcol=max_col;
	int size_sgmt=ssgmt.length();
	
	int flag=-1;
	int len=1;
	std::string spre_op="";
	std::string sop;
	while(iprow>0 && ipcol>0)
	{
		if(path[iprow][ipcol]==DIAG)
		{
			//check is decrease or increase 
			if(score_matrix[iprow-1][ipcol-1]<score_matrix[iprow][ipcol])//increase
				sop=SMATCH;
			else//decrease
				sop=SMISMATCH;
			iprow--;
			ipcol--;
		}
		else if(path[iprow][ipcol]==TOP)
		{
			sop=SINSERTION;
			iprow--;
		}
		else if(path[iprow][ipcol]==LEFT)
		{
			sop=SDELETION;
			ipcol--;
		}
		else
		{//START
			sop=SSTART;
			break;
		}

		

		if(sop!=spre_op)
		{
			//save into vector
			if(spre_op!="")
				cigar.push_back(std::make_pair(len,spre_op));
			len=1;
		}
		else
			len++;
	

		spre_op=sop;
	}//end of while 
	cigar.push_back(std::make_pair(len,spre_op));
}


//create a size_ptn* size_ref matrix
void SmithWaterman::loadMatrix(int** &matrix, int size_ref, int size_sgmt)
{
	score_matrix=new int*[size_sgmt+1];
	path=new int*[size_sgmt+1];
	for(int i=0;i<=size_sgmt;i++)
	{
		score_matrix[i]=new int[size_ref+1];
		path[i]=new int[size_ref+1];
	}
}

void SmithWaterman::initMatrix(int size_ref, int size_sgmt)
{
	for(int i=0;i<size_sgmt;i++)
		score_matrix[i][0]=0;

	for(int j=0;j<size_ref;j++)
		score_matrix[0][j]=0;
}

//free memory
void SmithWaterman::freeMatrix(int** &matrix, int size_sgmt)
{
	for(int i=0;i<=size_sgmt;i++)
	{
		delete[] score_matrix[i];
		delete[] path[i];
	}
	delete[] score_matrix;
	delete[] path;
}