//
//  GenSeqsUtils.h
//  
//
//  Created by Yufeng Wu on 11/2/14.
//  Generic sequencing related utils
//

#ifndef ____GenSeqsUtils__
#define ____GenSeqsUtils__

#include <string>
#include <vector>
#include <set>
using namespace std;

bool IsMissing(char nt);
bool IsGap(char nt);
char GetComplement(char bp);
char GetBaseUpper(char bp);
int ConvBaseToInt(char bp);
char ConvIntToBase(int bpInt);
void GetBasesFromSeq(const string &seqIn, int posStart, int numBases, bool fRight, bool fSkipGap, string &strSub);
void GetHammingDist1NgbrForSeq(const string &seqIn, set<string> &setNgbrs);
string ConsConsensusSeq( const vector<pair<string,double> > &listSeqsToMerge );
void GetAllSeqsShiftByOne( const string &seqIn, set<string> &setSeqsShift1 );
string GetReverseCompSeq(const string &seqIn);

#endif /* defined(____GenSeqsUtils__) */
