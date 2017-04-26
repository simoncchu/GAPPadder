//
//  fastaMultiSeqs.h
//  
//
//  Created by Yufeng Wu on 11/2/14.
//  a list of fasta sequences
//	
//  Edit: add setFastqSeqs function by Chong Chu on 08/16/15
//

#ifndef ____fastaMultiSeqs__
#define ____fastaMultiSeqs__

#include "fastareader.h"
#include <vector>
using namespace std;

// ******************************************************************
// represent a list of fasta sequences

class MultiFastqSeqs
{
public:
    MultiFastqSeqs();
    MultiFastqSeqs(const MultiFastqSeqs &rhs);
    ~MultiFastqSeqs();
	void setFastqSeqs(const MultiFastqSeqs &rhs);
    void ReadFromFile(const char *fileName);
    void Dump(bool fSeqOnly = false, int lenDump=60);
    int GetNumOfSeqs() const { return listFastaSeqPtrs.size(); }
    FastaSequence *GetSeq(int i) { return listFastaSeqPtrs[i]; }
    void EraseSeq(FastaSequence *pseq);
    void Append( FastaSequence *pseq ) { listFastaSeqPtrs.push_back( pseq ); }
    
private:
    void Reset();
    
    vector<FastaSequence *> listFastaSeqPtrs;
};

// ******************************************************************
// read from a file but instead of holding all sequences, only read
// one each time

class MultiFastqSeqsStream
{
public:
    MultiFastqSeqsStream(const char *fileName);
    ~MultiFastqSeqsStream();
    bool ReadOneSeq();  // return false when done
    void Dump();
    FastaSequence *GetCurSeq() { return ptrCurSeq; }
    
private:
    void Reset();
    
    const char *fileNameSeqs;
    ifstream *ptrfstream;
    FastaSequence * ptrCurSeq;
};


#endif /* defined(____fastaMultiSeqs__) */
