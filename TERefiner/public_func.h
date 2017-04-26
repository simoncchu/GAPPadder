
#ifndef _CORE_PUBLIC_FUNC_H_
#define _CORE_PUBLIC_FUNC_H_

#include<string>
#include<sstream>
#include<string>
#include<vector>

class PubFuncs
{
public:
	static std::string cvtInt2Str(int i); //Convert int to string 
		
	static int cvtStr2Int(std::string str); //Convert string to int 

	static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	static std::vector<std::string> split(const std::string &s, char delim);
	//static bool fileExist(std::string fpath);// Check whether a file exist or not 
};

#endif