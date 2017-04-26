#ifndef FASTAREADER_H
#define FASTAREADER_H
#include <iostream>
#include <string>
#include <memory>
#include <cstdio>
#include <cstdlib>
#include "abstractcharsequence.h"
using namespace std;

//bool IsMissing(char ch);
void THROW(const char *str);

class FastaReader;

class FastaSequence:public AbstractCharSequence
{
private:
    std::string _name;
    std::string _seq;

public:
    FastaSequence();
    static const int32_t DEFAULT_LINE_LENGTH;
    FastaSequence(const FastaSequence& cp);
    virtual ~ FastaSequence();
    virtual char at(int32_t index) const;
    virtual int32_t size() const;
    const char* name() const;
    void ResetName() { _name.clear(); }
    const char* c_str() const;
    const char* c_str(int posStart) const;
    void Clear();
    void SetName(const char *name) { _name = name; }
    void SetSeq(const string &str) { _seq = str; }
    void AddChar( char base );
    void AppendSeq(const std::string &str);
    void GetSubSequence(int posStart, int posEnd, FastaSequence &strSub) const;
    FastaSequence& operator=(const FastaSequence& cp);
    void printFasta(std::ostream& out,int32_t lineLength) const;
    void printFasta(std::ostream& out);
    //static char GetComplement(char bp);
    void RevsereComplement();
    void Reverse();
    
    friend class FastaReader;
};

class FastaReader
{
private:
    std::size_t _reserve;
    bool to_upper;
public:
    FastaReader();
    ~FastaReader();
    FastaReader& reserve(int32_t len);
    FastaReader& toupper(bool choice);
    FastaSequence* next(std::istream& in);
};

#endif
