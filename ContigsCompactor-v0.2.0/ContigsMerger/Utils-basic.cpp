//
//  Utils-basic.cpp
//  
//
//  Created by Yufeng Wu on 11/21/14.
//
//

#include "Utils-basic.h"
#include <iostream>
#include <cstdlib>
using namespace std;

void YW_ASSERT_INFO2(bool f, const char *info)
{
    if( f == false )
    {
        cout << info << endl;
        exit(1);
    }
}

