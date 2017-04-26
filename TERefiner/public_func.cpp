
#include"public_func.h"
#include<fstream>
#include<vector>

/*
Description: Convert int to string 
*/
std::string PubFuncs::cvtInt2Str(int i)
{
	std::stringstream ss;
	ss<<i;
	return ss.str();
}


/*
Description:  Convert string to int. 
*/
int PubFuncs::cvtStr2Int(std::string str)
{
	int num=0;
	int len = str.length();
	
	for(int i=0;i<len;i++)
	{
		num*=10;
		num+=(str[i]-'0');
	}

	return num; 
}

///*
//Description: Check whether a file exist or not.
//*/
//bool PubFuncs::fileExist(std::string filename)
//{
//	std::ifstream freads(filename);
//	if(freads)
//		return true;
//	else
//		return false;
//}

std::vector<std::string> &PubFuncs::split(const std::string &s, char delim, std::vector<std::string> &elems) 
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> PubFuncs::split(const std::string &s, char delim) 
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
