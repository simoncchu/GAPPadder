#ifndef _H_SDAERPILCDRAH_
#define _H_SDAERPILCDRAH_
#include<string>

class HardClipReads
{
public:
	HardClipReads(){}
	
public:
	void traceRawReadById(std::string fleft_fastq, std::string fright_fastq);
};

#endif
