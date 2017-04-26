#ifndef _H_TNEMNGILA_
#define _H_TNEMNGILA_

#include"bam_alignment_record.h"
#include<string>
#include<vector>
#include<utility>

class Alignment
{
public:
	Alignment(){}
	Alignment(BamAlignmentRecord*& bar);

	void setBar(BamAlignmentRecord*& bar);

public:
	int getReadType();
	int getPEReadType();
	int getClipType(int& pos1, int& len1, int& pos2, int& len2);
	std::string getCigar();
	void str2Cigar(std::string cigar, std::vector<std::pair<std::string, int> >& vcigar);
	
	bool isFullyMapped();
	bool isClipped();
	bool isHardClipped();
	bool isDuplicate();//whether optimal duplicate 
	bool isPrimaryAlign();//whether is primary alignment
	bool passQualityCK();//whether fails quality check 
	bool isMateMapped();//whether mate read is mapped or not
	bool isFirstInPair();//for pair-end reads, it is the first in pair.
	bool isSecondInPair();//for pair-end reads, it is the second in pair.
	bool isSecondaryAlignment();//whether it is secondary alignment
	bool isReverseStrand();
	bool isMateReverseStrand();

	int getReadType(int rlength);
	bool isFullyMapped(int rlength);//for those with varied read length
	bool isPerfectMapped(int rlength);//perfect mapped
	bool isClipped(int rlength);

private:
	BamAlignmentRecord* bar;
	std::string cigar;
};


#endif
