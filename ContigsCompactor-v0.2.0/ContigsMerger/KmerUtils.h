//
//  KmerUtils.h
//  
//
//  Created by Yufeng Wu on 8/31/14.
//
//  Process kmer related stuff
//

#ifndef ____KmerUtils__
#define ____KmerUtils__

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include "stdint.h"
using namespace std;


//////////////////////////////////////////////////////////////////////
// define kmer that is short than 32 nts

typedef char ReadNtType;
typedef uint64_t KmerTypeShort;
typedef string KmerTypeLong;

// define the map for k-mer frequency
typedef map<KmerTypeShort, double> MapShortKmerFreq;
typedef map<KmerTypeLong, double> MapLongKmerFreq;

//////////////////////////////////////////////////////////////////////
// init kmer by a pointer to an array, and a position
void FormKmerTypeShortSeg( const ReadNtType *pArrayIn, int posKmer, int kmerLen,KmerTypeShort &kmer );
// or create by another kmer and then left-shift by one nt
void FormKmerTypeShortShift( const KmerTypeShort &kmerTypePrev, int kmerLen, ReadNtType ntRightMost, KmerTypeShort &kmerNew);
// get all kmers from a short sequence
void GetAllKmersFromSeq( const ReadNtType *pSeq, int seqLen, int kmerLen, vector<KmerTypeShort> &listKmers);
// dump kmer out
void DumpKmer( const KmerTypeShort &kmer, int kmerLen );
// convert integer-based kmer to string representation
void ConvKmerToString( const KmerTypeShort &kmer, int szKmer, char *strKmer );
// read a set of kmer from file and store in a hash table; return kmer length
int ReadInKmerToHashMap( const char *fileName, MapShortKmerFreq &mapKmerFreq );
// test whether a new seq (i.e. read) contain at least certain number
// of frequent k-mer as in the map
bool IsReadContainingFreqKmers(  const ReadNtType *pSeq, int seqLen, int kmerLen, int minOccurKmerThres, const MapShortKmerFreq &mapKmerFreqIn);
// write out a read
void OutputReadFastq(ostream &outputStream, const ReadNtType *pSeq, int seqLen, const char *idReadDesc, const char *qualityRead);
bool IsKmerInKmerHashMap(const KmerTypeShort &kmer, MapShortKmerFreq &mapKmerFreq, double &freqKmer);
int GetKmerLength( MapShortKmerFreq &mapKmerFreq );
void AddShortKmerToHashMap( const KmerTypeShort &kmer, MapShortKmerFreq &mapKmerFreq, double freq );


//////////////////////////////////////////////////////////////////////
// Long kmer utilities
int ReadInLongKmerToHashMap( const char *fileName, MapLongKmerFreq &mapKmerFreq );
void AddLongKmerToHashMap( const KmerTypeLong &kmer, MapLongKmerFreq &mapKmerFreq, double freq );
bool IsKmerInKmerHashMap(const KmerTypeLong &kmer, MapLongKmerFreq &mapKmerFreq, double &freqKmer);
int GetLongKmerLength( MapLongKmerFreq &mapKmerFreq );
void GetLongKmersFromSeq( const ReadNtType *pSeq, int seqLen, int kmerLen, set<KmerTypeLong> &setKmers );

//////////////////////////////////////////////////////////////////////
// common class

//class AbstractKmerType
//{
//public:
//    virtual ;
//};



#endif /* defined(____KmerUtils__) */
