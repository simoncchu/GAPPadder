//
//  KmerUtils.cpp
//  
//
//  Created by Yufeng Wu on 8/31/14.
//
//

#include "KmerUtils.h"
#include <bitset>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib> 
#include <cstring>


//////////////////////////////////////////////////////////////////////
// define kmer that is short than 32 nts
// kmer: correspond nts from left to right; that is, nt[0] is the MSB

// set the kmer to be the specific nt at position
static void SetKmerTypeForNtAt( KmerTypeShort &kmer, int pos, ReadNtType nt )
{
    // each position: 2 bits in the representation
    // YW: if nt is bogus, we simply treat it as a good nt of random type (say A)
    const int szKmerType = 8*sizeof(KmerTypeShort);
    int posSet = szKmerType-2*pos-1;
    int bit1 = 0, bit2 = 0;
    if( nt == 'c' || nt == 'C'  )
    {
        bit2 = 1;
    }
    else if( nt == 'g' || nt == 'G'  )
    {
        bit1 = 1;
    }
    else if( nt == 't' || nt == 'T'  )
    {
        bit1 = 1;
        bit2 = 1;
    }
    
    uint64_t isfhit1 = 1;
    isfhit1 = ~(isfhit1 << posSet);
    uint64_t isfhit2 = 1;
    isfhit2 = ~(isfhit2 << (posSet-1));
    uint64_t isfhit3 = bit1;
    isfhit3 = (isfhit3 << posSet);
    uint64_t isfhit4 = bit2;
    isfhit4 = isfhit4 << (posSet-1);

    kmer &= isfhit1;
    kmer &= isfhit2;
    kmer |= isfhit3;
    kmer |= isfhit4;
//cout << "szKmerType: " << szKmerType << ", pos = "<< pos << ", poset = " << posSet << ":  nt = " << nt << ", bit1 = " << bit1 << ", bit2 = " << bit2 << ", kmer = ";
    //DumpKmer( kmer, 8 );
}

// init kmer by a pointer to an array, and a position
void FormKmerTypeShortSeg( const ReadNtType *pArrayIn, int posKmer, int kmerLen, KmerTypeShort &kmer )
{
    // caution: don't check for out-of-bound
    kmer = 0;
    for(int i=0; i<kmerLen; ++i)
    {
        SetKmerTypeForNtAt( kmer, i, pArrayIn[i+posKmer] );
    }
}

// or create by another kmer and then left-shift by one nt
void FormKmerTypeShortShift( const KmerTypeShort &kmerTypePrev, int kmerLen, ReadNtType ntRightMost, KmerTypeShort &kmerNew)
{
    //
    kmerNew = kmerTypePrev;
    kmerNew = kmerNew << 2;
//cout << "*** kmerprev: ";
//DumpKmer( kmerTypePrev, 8);
//cout << ", ntNew = " << ntRightMost << endl;
    // clear the top2 bits
    //const int szKmerType = 8*sizeof(KmerTypeShort);
    //int leftmostposSet = szKmerType-1-2*kmerLen;
    //kmerNew &= ~(1 << leftmostposSet);
    //kmerNew &= ~(1 << (leftmostposSet-1));
    // set the rightmost positions
    SetKmerTypeForNtAt( kmerNew, kmerLen-1, ntRightMost );
}

// get all kmers from a short sequence
void GetAllKmersFromSeq( const ReadNtType *pSeq, int seqLen, int kmerLen, vector<KmerTypeShort> &listKmers)
{
    // YW: caution: assume kmer len is < seqLen!!!
    if( seqLen < kmerLen)
    {
        cout << "Fatal error: read is too short for kmer\n";
        exit(1);
    }
    
    
    listKmers.clear();
    // get the first by directly getting it
    KmerTypeShort ktFirst;
    FormKmerTypeShortSeg( pSeq, 0, kmerLen, ktFirst );
    listKmers.push_back( ktFirst );
    KmerTypeShort ktLast = ktFirst;
    for(int i=1; i<seqLen-kmerLen+1; ++i)
    {
        //
        KmerTypeShort ktcur;
        int posKRightMost = i + kmerLen-1;
        FormKmerTypeShortShift( ktLast, kmerLen, pSeq[ posKRightMost ], ktcur );
        listKmers.push_back(ktcur);
        ktLast = ktcur;
    }
}

// dump kmer out
void DumpKmer( const KmerTypeShort &kmer, int kmerLen )
{
    // for now assume it is 8 bit
    //std::bitset<8> y( kmer );
    //cout << y;
    cout << std::hex << setw(8) << kmer;
}

// convert integer-based kmer to string representation
void ConvKmerToString( const KmerTypeShort &kmer, int szKmer, char *strKmer )
{
    // YW: caution: the code assume buffer size is correct; that is, one
    // needs the one last additional \0
    //std::ostringstream o;
    //o << kmer;
    //string str;
    //str = o.str();
    //cout << "-- converted string: " << o.str() << "  ";
    //strncpy(strKmer, o.str().c_str(), sizeof(KmerTypeShort) );
    //strcpy(strKmer, o.str().c_str() );
    const int szKmerType = 8*sizeof(KmerTypeShort);
    for(int i=0; i<szKmer; ++i)
    {
        int posSet = szKmerType-2*i-1;
        // get the two bits of the kmer
        KmerTypeShort kt1 = 1;
        kt1 = kt1 << posSet;
        KmerTypeShort kt2 = 1;
        kt2 = kt2 << (posSet-1);
        
        bool flagk1Zero = ( kt1 & kmer ) == 0;
        bool flagk2Zero = ( kt2 & kmer ) == 0;
        
        if( flagk1Zero == true && flagk2Zero == true )
        {
            strKmer[i] = 'A';
        }
        else if( flagk1Zero == true && flagk2Zero == false )
        {
            strKmer[i] = 'C';
        }
        else if( flagk1Zero == false && flagk2Zero == true )
        {
            strKmer[i] = 'G';
        }
        else
        {
            strKmer[i] = 'T';
        }
    }
    strKmer[szKmer] = 0;
}

// read a set of kmer from file and store in a hash table 
int ReadInKmerToHashMap( const char *fileName, MapShortKmerFreq &mapKmerFreq )
{
    // input file (with filename) contains the list of kmer and its freq
    // in the format: <kmer in ACGT>, <freq in int or double>
    ifstream inKmerFile;
    inKmerFile.open( fileName );
    
    int resKmerLen = -1;
    
    // read line by line
    while( inKmerFile.eof() == false )
    {
        char buf[10240];
        inKmerFile.getline( buf, sizeof(buf) );
        if( strlen(buf) == 0 )
        {
            break;
        }
        std::stringstream bufStream(buf);
        string kmer;
        int kmerFreq;
        bufStream >> kmer >> kmerFreq;
        //cout << "Found one kmer: " << kmer << ", with freq: " << kmerFreq << endl;
        if( resKmerLen < 0 )
        {
            resKmerLen = kmer.length();
        }
        
        // now form the kmer from the string
        KmerTypeShort kmerConv;
        FormKmerTypeShortSeg( kmer.c_str(), 0, kmer.length(), kmerConv );
        
        if( mapKmerFreq.find( kmerConv ) != mapKmerFreq.end() )
        {
            cout << "Warning: kmer has been found before: " << kmer << endl;
        }
        mapKmerFreq.insert( MapShortKmerFreq::value_type( kmerConv, kmerFreq ) );
    }
    return resKmerLen;
}

// test whether a new seq (i.e. read) contain at least certain number
// of frequent k-mer as in the map 
bool IsReadContainingFreqKmers(  const ReadNtType *pSeq, int seqLen, int kmerLen, int minOccurKmerThres, const MapShortKmerFreq &mapKmerFreqIn)
{
    // first form all the k-mer in the desired length
    vector<KmerTypeShort> listKermInSeq;
    GetAllKmersFromSeq( pSeq, seqLen, kmerLen, listKermInSeq);
    int numOccurs = 0;
    for(int i=0; i<(int)listKermInSeq.size(); ++i)
    {
        if( mapKmerFreqIn.find(listKermInSeq[i] ) != mapKmerFreqIn.end() )
        {
//cout << "** This kmer is frequent: ";
//char kmerChar[1024];
//ConvKmerToString( listKermInSeq[i], kmerLen, kmerChar  );
//cout << kmerChar << endl;
            ++numOccurs;
        }
    }
    if( numOccurs >= minOccurKmerThres )
    {
        //cout << "^^^^^^ This read is frequent: " << pSeq << endl;
        return true;
    }
    else
    {
        return false;
    }
}

// write out a read
void OutputReadFastq(ostream &outputStream, const ReadNtType *pSeq, int seqLen, const char *idReadDesc, const char *qualityRead)
{
    //
    outputStream << "@" << idReadDesc << endl;
    for(int i=0; i<seqLen; ++i)
    {
        outputStream << pSeq[i];
    }
    outputStream << endl;
    outputStream << "+\n";
    outputStream << qualityRead << endl;
}

bool IsKmerInKmerHashMap(const KmerTypeShort &kmer, MapShortKmerFreq &mapKmerFreq, double &freqKmer)
{
    //
    bool res = mapKmerFreq.find(kmer) != mapKmerFreq.end();
    if(res == true)
    {
        freqKmer = mapKmerFreq[kmer];
    }
    return res;
}

void AddShortKmerToHashMap( const KmerTypeShort &kmer, MapShortKmerFreq &mapKmerFreq, double freq )
{
    double freqExist = 0.0;
    bool fExist = IsKmerInKmerHashMap( kmer, mapKmerFreq, freqExist );
    if( fExist == false )
    {
        // add an empty entry
        mapKmerFreq.insert( MapShortKmerFreq :: value_type(kmer, 0.0) );
    }
    mapKmerFreq[kmer] += freq;
}

//////////////////////////////////////////////////////////////////////
// Long kmer utilities
int ReadInLongKmerToHashMap( const char *fileName, MapLongKmerFreq &mapKmerFreq )
{
    //
    // input file (with filename) contains the list of kmer and its freq
    // in the format: <kmer in ACGT>, <freq in int or double>
    ifstream inKmerFile;
    inKmerFile.open( fileName );
    
    int resKmerLen = -1;
    
    // read line by line
    while( inKmerFile.eof() == false )
    {
        char buf[10240];
        inKmerFile.getline( buf, sizeof(buf) );
        if( strlen(buf) == 0 )
        {
            break;
        }
        std::stringstream bufStream(buf);
        string kmer;
        int kmerFreq;
        bufStream >> kmer >> kmerFreq;
//cout << "Found one kmer: " << kmer << ", with freq: " << kmerFreq << endl;
        if( resKmerLen < 0 )
        {
            resKmerLen = kmer.length();
        }
        
        // now form the kmer from the string        
        if( mapKmerFreq.find( kmer ) != mapKmerFreq.end() )
        {
            cout << "Warning: kmer has been found before: " << kmer << endl;
        }
        mapKmerFreq.insert( MapLongKmerFreq::value_type( kmer, kmerFreq ) );
    }
    return resKmerLen;

}

bool IsKmerInKmerHashMap(const KmerTypeLong &kmer, MapLongKmerFreq &mapKmerFreq, double &freqKmer)
{
    //
    bool res = mapKmerFreq.find(kmer) != mapKmerFreq.end();
    if(res == true)
    {
        freqKmer = mapKmerFreq[kmer];
    }
    return res;
}

void AddLongKmerToHashMap( const KmerTypeLong &kmer, MapLongKmerFreq &mapKmerFreq, double freq )
{
    //
    double freqExist = 0.0;
    bool fExist = IsKmerInKmerHashMap( kmer, mapKmerFreq, freqExist );
    if( fExist == false )
    {
        // add an empty entry
        mapKmerFreq.insert( MapLongKmerFreq :: value_type(kmer, 0.0) );
    }
    mapKmerFreq[kmer] += freq;
}

int GetLongKmerLength( MapLongKmerFreq &mapKmerFreq )
{
    //
    if( mapKmerFreq.size() == 0 )
    {
        //
        cout << "FATAL ERROR: map of kmers is empty.\n";
        exit(1);
    }
    return ( mapKmerFreq.begin()->first ).length();
}

void GetLongKmersFromSeq( const ReadNtType *pSeq, int seqLen, int kmerLen, set<KmerTypeLong> &setKmers )
{
    // YW: caution: assume kmer len is < seqLen!!!
    if( seqLen < kmerLen)
    {
        cout << "Fatal error: read is too short for kmer\n";
        exit(1);
    }
    
    //
    setKmers.clear();
    string strPrev;
    // get the first by directly getting it
    for(int i=0; i<seqLen-kmerLen+1; ++i)
    {
        // first string: just do it
        string strCur;
        int posKRightMost = i + kmerLen-1;
        if( i == 0 )
        {
            for( int kk=i; kk<=posKRightMost; ++kk  )
            {
                //
                strCur += pSeq[kk];
            }
        }
        else
        {
            //do by shift
            strCur = strPrev.substr( 1, kmerLen-1 );
            strCur.push_back( pSeq[posKRightMost] );
        }
        setKmers.insert(strCur);
        
        strPrev = strCur;
    }

}

