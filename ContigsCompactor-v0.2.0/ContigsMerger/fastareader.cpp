#include <cstdio>
#include <cctype>
#include <algorithm>
//#include "throw.h"
#include "fastareader.h"
#include "GenSeqsUtils.h"
using namespace std;

// **********************************************************************************

void THROW(const char *str)
{
    std::cout << "FATAL ERROR: " << str << std::endl;
    exit(1);
}

//bool IsMissing(char ch)
//{
//    //
//    return ch == 'N' || ch == 'n';
//}

// **********************************************************************************

const int32_t FastaSequence::DEFAULT_LINE_LENGTH=60;

FastaSequence::FastaSequence() : _name("uninitialized"), _seq("") {}
FastaSequence::FastaSequence(const FastaSequence& cp):_name(cp._name),_seq(cp._seq)
{
}
FastaSequence::~FastaSequence()
{
}
char FastaSequence::at(int32_t index) const
{
    return _seq.at(index);
}
int32_t FastaSequence::size() const
{
    return (int32_t)_seq.size();
}
const char* FastaSequence::name() const
{
    return _name.c_str();
}
const char* FastaSequence::c_str() const
{
    return _seq.c_str();
}

const char*  FastaSequence :: c_str(int posStart) const
{
    return _seq.c_str()+posStart;
}

FastaSequence& FastaSequence::operator=(const FastaSequence& cp)
{
    if(this!=&cp)
    {
        _name.assign(cp._name);
        _seq.assign(cp._seq);
    }
    return *this;
}
void FastaSequence::printFasta(std::ostream& out,int32_t lineLength) const
{
    out << ">" << _name;
    for(int32_t i=0;i< size();++i)
    {
        if(i%lineLength==0) out << std::endl;
        out << at(i);
    }
    out << std::endl;
}
void FastaSequence::printFasta(std::ostream& out)
{
    printFasta(out,DEFAULT_LINE_LENGTH);
}

void FastaSequence :: RevsereComplement()
{
    // do the reverse complement
    string strNew;
    for(int i=(int)_seq.length()-1; i>=0; --i)
    {
        char bpComp = GetComplement( _seq[i] );
        strNew.push_back(bpComp);
    }
    // assign the new string
    _seq = strNew;
}

void FastaSequence :: Reverse()
{
    // just reverse, no complememnt
    std::reverse(_seq.begin(),_seq.end());
}

void FastaSequence :: Clear()
{
    _seq.clear();
}
void FastaSequence :: AddChar( char base )
{
    _seq += base;
}

void FastaSequence :: AppendSeq(const std::string &str)
{
    _seq.append(str);
}

void FastaSequence :: GetSubSequence(int posStart, int posEnd, FastaSequence &strSub) const
{
    strSub.Clear();
    for(int pos=posStart; pos<=posEnd; ++pos)
    {
        strSub.AddChar( at(pos) );
    }
}

#if 0
char FastaSequence :: GetComplement(char bp)
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
}
#endif

// *******************************************************************************


FastaReader::FastaReader():_reserve(BUFSIZ),to_upper(false)
{
}
FastaReader::~FastaReader()
{
}

FastaReader& FastaReader::reserve(int32_t len)
{
    this->_reserve=(len<=0?BUFSIZ:len);
    return *this;
}

FastaReader& FastaReader::toupper(bool choice)
{
    this->to_upper=choice;
    return *this;
}

//std::auto_ptr<FastaSequence> FastaReader::next(std::istream& in)
FastaSequence* FastaReader::next(std::istream& in)
{
    FastaSequence *ret = NULL;
    //std::auto_ptr<FastaSequence> ret(0);
    if(!in.good() || in.eof()) return ret;
    int c;
    
    while((c=in.get())!=EOF)
    {
        if(c=='>')
        {
            //if(ret.get()!=0)
            if( ret != NULL )
            {
                in.unget();
                return ret;
            }
            ret = new FastaSequence;
            ret->ResetName();
            //ret.reset(new FastaSequence);
            ret->_seq.reserve(_reserve);
            while((c=in.get())!=EOF && c!='\n')
            {
                if(c=='\r') continue;
                ret->_name+=(char)c;
            }
            continue;
        }
        if(std::isspace(c)) continue;
        if(!std::isalpha(c))
        {
            char buf[10240];
            sprintf(buf, "Bad char in sequence %c", (char)c);
            THROW(buf ) ;
        }
        if( ret == NULL )
        {
            THROW("header missing");
        }
        //if(ret.get()==0) THROW("header missing");
        if(to_upper) c=std::toupper(c);
        ret->_seq+=(char)c;
    }
    return ret;
}

