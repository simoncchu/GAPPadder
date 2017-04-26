#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>

#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <cmath>
using namespace std;


#include <cstring>
#include <sstream>
//#include "fastareader.h"
#include "fastaMultiSeqs.h"
//#include "ScaffoldGapFiller.h"
#include "ContigsCompactor.h"

const char *versionCrecoal = "**************************************************************************\n            ContigsMerger ver. 0.1.4, released on March 11, 2015 \n**************************************************************************\n\n";
//static bool fTestMode = true;
static int repeatfileArgIndex = 1;
//static int asmfileArgIndex = 2;
static bool fVerbose = false;
static double minFracOverlap = 0.005;
static double minOverlapLen = 100000;
static double minOverlapLenWithScaffold=6;
//static double minOverlapLen01 = 12;
static double maxOverlapClipLen = 0;
static double maxFracScoreLoss = 0.01;
static int numOfThread=6;
static int minSupportKmer=5;
static int lengthContigOutputLine = 60;
static int defQuickKmerLen = 10;
//static bool fNewSeqsOnly = true;
static double scoreMismatch = -1.0;
static double scoreIndel = -1.0;
static string fileContigsInfo = "tmp.info";
static int maxContigPathLen = -1;
static int maxCountContigInPath = -1;
static string fileScaffoldInfo="";
static int supportPairsCutoff=2;

//SNPDistInfo SNPDefaultDistInfo;  // sites position info

static void Usage()
{
    cout << "Usage: ./ContigsMerger <repeat contigs>\n";
    exit(1);
}


static bool CheckArguments(int argc, char **argv) 
{
    if( argc <= 1 )
    {
        return true;
    }
    
    // Check argument one by one
    int argpos = 1;
    while( argpos < argc)
    {
        //
        if( argv[ argpos ][ 0 ] != '-' )
        {
            // must have this
            //return false;
            repeatfileArgIndex = argpos;
            ++argpos;
        }
        
        else if( argv[argpos][1] == 'V' )
        {
            argpos++;
            fVerbose = true;
cout << "Turn on Verbose\n";
            //cout << "Species tree name: " << fileNameSpecies << endl;
        }
        else if( argv[argpos][1] == 'l' )
        {
            //
            argpos ++;
            sscanf(argv[argpos], "%d", &lengthContigOutputLine );
            argpos++;
        }
        else if( argv[argpos][1] == 's' )
        {
            //
            argpos ++;
            float fval;
            sscanf(argv[argpos], "%f", &fval );
            maxFracScoreLoss = fval;
            argpos++;
        }
        else if( argv[argpos][1] == 'c' )
        {
            //
            argpos ++;
            float fval;
            sscanf(argv[argpos], "%f", &fval );
            minFracOverlap = fval;
            argpos++;
        }
        else if( argv[argpos][1] == 'x' )
        {
            //
            argpos ++;
            float fval;
            sscanf(argv[argpos], "%f", &fval );
            minOverlapLen = fval;
            argpos++;
        }
		else if (argv[argpos][1] == 'y')
		{
			//
			argpos++;
			float fval;
			sscanf(argv[argpos], "%f", &fval);
			maxOverlapClipLen = fval;
			argpos++;
		}
		else if( argv[argpos][1] == 'm' )
        {
            //
            argpos ++;
            int fval;
            sscanf(argv[argpos], "%d", &fval );
            minSupportKmer = fval;
            argpos++;
        }
		else if( argv[argpos][1] == 't' )
        {
            //
            argpos ++;
            int fval;
            sscanf(argv[argpos], "%d", &fval );
            numOfThread = fval;
            argpos++;
        }
		else if( argv[argpos][1] == 'z' )
        {
            //
            argpos ++;
            float fval;
            sscanf(argv[argpos], "%f", &fval );
            minOverlapLenWithScaffold = fval;
            argpos++;
        }
        else if( argv[argpos][1] == 'k' )
        {
            //
            argpos ++;
            int fval;
            sscanf(argv[argpos], "%d", &fval );
            defQuickKmerLen = fval;
            argpos++;
        }
        else if( argv[argpos][1] == 'i' && argv[argpos][2] == '1' )
        {
            //
            argpos ++;
            float fval;
            sscanf(argv[argpos], "%f", &fval );
            scoreMismatch = fval;
            argpos++;
        }
        else if( argv[argpos][1] == 'i' && argv[argpos][2] == '2' )
        {
            //
            argpos ++;
            float fval;
            sscanf(argv[argpos], "%f", &fval );
            scoreIndel = fval;
            argpos++;
        }
        else if( argv[argpos][1] == 'o' )
        {
            //
            argpos ++;
            fileContigsInfo = argv[argpos];
            argpos++;
        }
        else if( argv[argpos][1] == 'p'  && argv[argpos][2] == '1' )
        {
            //
            argpos ++;
            int fval;
            sscanf(argv[argpos], "%d", &fval );
            maxContigPathLen = fval;
            argpos++;
        }
        else if( argv[argpos][1] == 'p'  && argv[argpos][2] == '2' )
        {
            //
            argpos ++;
            int fval;
            sscanf(argv[argpos], "%d", &fval );
            maxCountContigInPath = fval;
            argpos++;
        }
		else if( argv[argpos][1] == 'e' )
		{
			argpos ++;
			fileScaffoldInfo = argv[argpos];
			argpos ++;
		}
		else if(argv[argpos][1] == 'u')
		{
			argpos ++;
			int fval;
            sscanf(argv[argpos], "%d", &fval );
			supportPairsCutoff=fval;
			argpos++;
		}
        else
        {
            cout << "Wrong input.\n";
            exit(1);
        }
        //else if( argv[argpos][1] == 'n' )
        //{
        //    argpos++;
        //    fNewSeqsOnly = true;
        //    //cout << "Species tree name: " << fileNameSpecies << endl;
        //}
    }
    

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TestKmerUtils( const char *repeatFileName)
{
    // Test read the two fasta file for now
    MultiFastqSeqs repeatSeqs;
    repeatSeqs.ReadFromFile( repeatFileName );
//repeatSeqs.Dump();
    MultiFastqSeqs listContigsOrig( repeatSeqs );
    
    
    // test scaffold filling
    //ScaffoldGapFiller filler;
    //filler.FillScaffold( asmSeqs, repeatSeqs );

//cout << "Score mismatch: " << scoreMismatch << ", scoreIndel: " << scoreIndel << endl;
    
    ContigsCompactor compactor;
    compactor.SetVerbose( fVerbose );
    compactor.SetFractionLossScore(maxFracScoreLoss);
    compactor.SetMinOverlap( minFracOverlap );
	compactor.SetMaxOverlapLenClip( maxOverlapClipLen );
//cout<<minOverlapLen<<" "<<minOverlapLenWithScaffold<<endl;
    compactor.SetMinOverlapLen( minOverlapLen );
	compactor.SetMinOverlapLenWithScaffold(minOverlapLenWithScaffold);
    compactor.SetQuickCheckKmerLen( defQuickKmerLen );
    compactor.SetMismatchScore(scoreMismatch);
    compactor.SetIndelScore(scoreIndel);
	compactor.SetNumOfThreads(numOfThread);
	compactor.SetMinSupportKmers(minSupportKmer);

    if( maxContigPathLen > 0 )
    {
        compactor.SetMaxContigPathLen( maxContigPathLen );
    }
    if( maxCountContigInPath > 0 )
    {
        compactor.SetMaxCountPerContigInPaths( maxCountContigInPath );
    }
    
    //compactor.CompactVer2( repeatSeqs );
	compactor.CompactVer3( repeatSeqs, fileScaffoldInfo.c_str(),supportPairsCutoff);
    
    //cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
    //cout << "After compacting, the set of fasta files: \n"; 
    //if( fNewSeqsOnly == true )
    //{
    //
    compactor.GetNewSeqs().Dump( true, lengthContigOutputLine );
    //}
    //else
    //{
    listContigsOrig.Dump(true, lengthContigOutputLine);
    //}
    // also output contigs info
    compactor.OutputContigsInfoVer2( fileContigsInfo.c_str() );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
    // info
    //cout << versionCrecoal << endl;
    if( CheckArguments( argc, argv) == false)
    {
        Usage();
    }
    if( repeatfileArgIndex < 0 )
    {
        cout << "Reads file not specified.\n";
        exit(1);
    }
    
//TestGraph();
//exit(1);
    
    // test code only for now 
    TestKmerUtils( argv[repeatfileArgIndex] );
    //TestSeqUtils( argv[ readfileArgIndex ] );
    //FindReadsOnRepeats( argv[repeatfileArgIndex], argv[asmfileArgIndex] );
    
    // for now, do nothing
    return 0;
}
