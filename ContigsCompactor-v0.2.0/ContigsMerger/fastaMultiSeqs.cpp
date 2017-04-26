//
//  fastaMultiSeqs.cpp
//  
//
//  Created by Yufeng Wu on 11/2/14.
//
//  Edit: add setFastqSeqs function by Chong Chu on 08/16/15
//

#include "fastaMultiSeqs.h"
#include <iostream>
#include <fstream>

// ******************************************************************
// represent a list of fasta sequences


MultiFastqSeqs :: MultiFastqSeqs()
{
    // nothing for now
}


MultiFastqSeqs :: MultiFastqSeqs(const MultiFastqSeqs &rhs)
{
    for(int i=0; i<(int)rhs.listFastaSeqPtrs.size(); ++i )
    {
        Append( new FastaSequence( *rhs.listFastaSeqPtrs[i] ) );
    }
}


MultiFastqSeqs :: ~MultiFastqSeqs()
{
    Reset();
}

void MultiFastqSeqs :: setFastqSeqs(const MultiFastqSeqs &rhs)
{
	for(int i=0; i<(int)rhs.listFastaSeqPtrs.size(); ++i )
    {
        Append( new FastaSequence( *rhs.listFastaSeqPtrs[i] ) );
    }
}

void MultiFastqSeqs :: ReadFromFile(const char *fileName)
{
    Reset();
    
    //
    ifstream inFile;
    inFile.open( fileName );
    
    FastaReader fastaSeqsReader;
    while(true)
    {
        //
        FastaSequence *ps = fastaSeqsReader.toupper(true).next( inFile );
        if( ps != NULL )
        {
            // a new seq
            listFastaSeqPtrs.push_back(ps);
        }
        else
        {
            // that is it
            break;
        }
    }
    
    inFile.close();
}

void MultiFastqSeqs :: Dump(bool fSeqOnly, int lenDump)
{
    // dump out
    if( fSeqOnly == false )
    {
        cout << "** Number of fastq sequences: " << listFastaSeqPtrs.size() << endl;
    }
    for(int i=0; i<(int)listFastaSeqPtrs.size(); ++i)
    {
        //
        listFastaSeqPtrs[i]->printFasta(cout, lenDump);
    }
}

// free mem if needed
void MultiFastqSeqs :: Reset()
{
    for(int i=0; i<(int)listFastaSeqPtrs.size(); ++i)
    {
        //
        delete listFastaSeqPtrs[i];
    }
    listFastaSeqPtrs.clear();
}

// erase one element
void MultiFastqSeqs :: EraseSeq(FastaSequence *pseq)
{
    //
    vector<FastaSequence *>::iterator iter = listFastaSeqPtrs.begin();
    vector<FastaSequence *>::iterator endIter = listFastaSeqPtrs.end();
    for(; iter != endIter; ++iter)
    {
        if(*iter == pseq)
        {
            listFastaSeqPtrs.erase(iter);
            // also free the buffer as well
            delete pseq;
            
            return;
        }
    }

}

// ******************************************************************
// read from a file but instead of holding all sequences, only read
// one each time


MultiFastqSeqsStream :: MultiFastqSeqsStream(const char *fileName) : fileNameSeqs(fileName), ptrfstream(NULL), ptrCurSeq(NULL)
{
    //
}

MultiFastqSeqsStream :: ~MultiFastqSeqsStream()
{
    //
    Reset();
    // close file
    if( ptrfstream != NULL )
    {
        ptrfstream->close();
        delete ptrfstream;
    }
}

bool MultiFastqSeqsStream :: ReadOneSeq()
{
    // open file if not done yet
    if( ptrfstream == NULL )
    {
        ptrfstream = new ifstream;
        ptrfstream->open( fileNameSeqs );
    }
    
    //
    Reset();
    //
    FastaReader fastaSeqsReader;
    FastaSequence *ps = fastaSeqsReader.next( *ptrfstream );
    if( ps != NULL )
    {
        ptrCurSeq = ps;
        return true;
    }
    else
    {
        return false;
    }
}

void MultiFastqSeqsStream :: Dump()
{
    //
    if( ptrCurSeq != NULL )
    {
        ptrCurSeq->printFasta( cout );
    }
}
    

void MultiFastqSeqsStream :: Reset()
{
    //
    if( ptrCurSeq != NULL )
    {
        delete ptrCurSeq;
    }
    ptrCurSeq = NULL;
}
    


