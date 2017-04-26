//
//  GenSeqsUtils.cpp
//  
//
//  Created by Yufeng Wu on 11/2/14.
//  Generic sequencing related utils
//

#include "GenSeqsUtils.h"
#include <cctype>
#include <iostream>

// is the nt missing or not
bool IsMissing(char nt)
{
    return nt == 'N' || nt == 'n';
}

bool IsGap(char nt)
{
    return nt == '-';
}

char GetComplement(char bp)
{
    //
    if( IsMissing(bp) )
    {
        return bp;
    }
    char bpUse = toupper(bp);
    if( bpUse == 'A')
    {
        return 'T';
    }
    //else if(bp == 'a')
    //{
    //    return 't';
    //}
    else if( bpUse == 'T')
    {
        return 'A';
    }
    //else if( bpUse == 't')
    //{
    //    return 'a';
    //}
    else if( bpUse == 'G')
    {
        return 'C';
    }
    //else if(bp == 'g')
    //{
    //    ;
    //}
    else if(bpUse == 'C')
    {
        return 'G';
    }
    return 'N';
}

void GetBasesFromSeq(const string &seqIn, int posStart, int numBases, bool fRight, bool fSkipGap, string &strSub)
{
    // get a portion of based from some position of certian length
    // seqIn: the original string
    // posStart: where to get bases
    // numBases: how many bases to get
    // fRight: move to right?
    // fSkipGap: only A/G/C/T/N or gap is allowed
    strSub.clear();
    int curpos = posStart;
    while( curpos >=0 && curpos <(int)seqIn.length() && (int)strSub.length() < numBases)
    {
        char nt = seqIn[curpos];
        if(fSkipGap == false || IsGap(nt) == false )
        {
            strSub.append(1, nt);
        }
        if( fRight == true )
        {
            ++curpos;
        }
        else
        {
            --curpos;
        }
    }
    
}

char GetBaseUpper(char bp)
{
    if( bp == 'a' )
    {
        return 'A';
    }
    if( bp == 'c' )
    {
        return 'C';
    }
    if( bp == 'g' )
    {
        return 'G';
    }
    if( bp == 't' )
    {
        return 'T';
    }
    // assume it is already upper case
    return bp;
}

void GetHammingDist1NgbrForSeq(const string &seqIn, set<string> &setNgbrs)
{
    // assume the seq is a DNA sequence (in all captial letters), get all distance (hamming) 1 neighbor
    setNgbrs.clear();
    char listBases[4] = {'A', 'T',  'C', 'G'};
    for(int i=0; i<(int)seqIn.length(); ++i)
    {
        // change this base to other
        if( IsMissing(seqIn[i]) || IsGap(seqIn[i]) )
        {
            //
            cout << "WARNING: GetHammingDist1NgbrForSeq does not support gap or missing values" << endl;
            //exit(1);
        }
        
        //
        char bp = GetBaseUpper(seqIn[i]);
        for(int j=0; j<4; ++j)
        {
            if( bp != listBases[j] )
            {
                // now create a new string
                string seqNgbr = seqIn;
                seqNgbr[i] = listBases[j];
                setNgbrs.insert( seqNgbr );
            }
        }
    }
}

int ConvBaseToInt(char bp)
{
    // assume captial
    if( bp == 'A')
    {
        return 0;
    }
    if( bp == 'T')
    {
        return 1;
    }
    if( bp == 'C')
    {
        return 2;
    }
    if( bp == 'G')
    {
        return 3;
    }
    return 4;
}

char ConvIntToBase(int bpInt)
{
    if( bpInt == 0 )
    {
        return 'A';
    }
    if( bpInt == 1 )
    {
        return 'T';
    }
    if( bpInt == 2 )
    {
        return 'C';
    }
    if( bpInt == 3 )
    {
        return 'G';
    }
    return 'N';
}


string ConsConsensusSeq( const vector<pair<string,double> > &listSeqsToMerge )
{
    if( listSeqsToMerge.size() == 0 )
    {
        cout << "ConsConsensusSeq: cannot merge empty sequences\n";
        return "";
    }
    
    // given a list of sequences (with weights), construct the consensus sequence
    // here, assume each sequence is of same length
    for(int i=0; i<(int)listSeqsToMerge.size(); ++i)
    {
        if( listSeqsToMerge[i].first.length() != listSeqsToMerge[0].first.length()  )
        {
            cout << "ConsConsensusSeq: the sequences must be of same length.\n";
            return "";
        }
    }
    
    //
    string res;
    for(int pos=0; pos<(int)listSeqsToMerge[0].first.length(); ++pos )
    {
        //
        vector<double> listFreqBases(5);
        for(int i=0; i<5; ++i)
        {
            listFreqBases[i] = 0.0;
        }
        for(int i=0; i<(int)listSeqsToMerge.size(); ++i)
        {
            int bpInt = ConvBaseToInt( listSeqsToMerge[i].first[pos] );
            listFreqBases[bpInt] += listSeqsToMerge[i].second;
        }
        int bpMax = 0;
        double valMax = listFreqBases[0];
        for(int j=1; j<5; ++j)
        {
            if( listFreqBases[j] > valMax)
            {
                valMax = listFreqBases[j];
                bpMax = j;
            }
        }
        char bpMaxConv = ConvIntToBase(bpMax);
        res += bpMaxConv;
    }
    return res;
}

void GetAllSeqsShiftByOne( const string &seqIn, set<string> &setSeqsShift1 )
{
    // shift by one to the left and fill in a base at the right end
    if( seqIn.length() == 0 )
    {
        // fatal error
        cout << "WARNING: empty string in GetAllSeqsShiftByOne" << endl;
        return;
    }
    string strShift = seqIn.substr( 1, seqIn.length()-1 );
    char bases[] = { 'A','T','C','G' };
    for(int i=0; i<=3; ++i)
    {
        string strNew = strShift;
        strNew += bases[i];
        setSeqsShift1.insert( strNew );
    }
}

string GetReverseCompSeq(const string &seqIn)
{
    //
    // do the reverse complement
    string strNew;
    for(int i=(int)seqIn.length()-1; i>=0; --i)
    {
        char bpComp = GetComplement( seqIn[i] );
        strNew.push_back(bpComp);
    }
    // assign the new string
    return strNew;
}

