#include"contigs.h"
#include<string>

using namespace std;

Contigs::Contigs()
{
	contig_name="";
	breverse=0;
}

void Contigs::setContigName(string contig_name)
{
	this->contig_name=contig_name;
}

void Contigs::setOrientation(int breverse)
{
	this->breverse=breverse;
}

std::string Contigs::getContigName()
{
	return this->contig_name;
}

int Contigs::getOrientation()
{
	return this->breverse;
}

void Contigs::setContigID(int contig_id)
{
	this->contig_id=contig_id;
}

int Contigs::getContigID()
{
	return this->contig_id;
}