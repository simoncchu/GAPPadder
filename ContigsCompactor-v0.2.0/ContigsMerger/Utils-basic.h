//
//  Utils-basic.h
//  
//
//  Created by Yufeng Wu on 11/21/14.
//  Some basic utilties that are applicable to everything
//

#ifndef ____Utils_basic__
#define ____Utils_basic__

#include <set>
using namespace std;

// ******************************************************************

// compare two values to see if one is within some fraction of the other
template<class T> 
bool IsValWithinDiffRange(const T &val, const T &valCmp, double fracDiff)
{
    // fracDiff: fraction of differences, must be within positive and between 0 and 1
    if( val > valCmp)
    {
        return val < (1.0+fracDiff)*valCmp;
    }
    else
    {
        return val > (1.0-fracDiff)*valCmp;
    }
}

// interval length
template<typename T>
T CalcIntervalLen( const T& iv1, const T& iv2 )
{
    return iv2-iv1+1;
}

void YW_ASSERT_INFO2(bool f, const char *info);


template<class TYPE>
bool IsSetContainerGen1(const set<TYPE> &container, const set<TYPE> &contained)
{
    //
    for(typename set<TYPE> ::iterator it = contained.begin(); it != contained.end(); ++it)
    {
        if(  container.find( *it ) == container.end()  )
        {
            return false;
        }
    }
    return true;
}


#endif /* defined(____Utils_basic__) */
