#include <sstream>
#include "abstractcharsequence.h"


AbstractCharSequence::AbstractCharSequence()
{
}
AbstractCharSequence::~AbstractCharSequence()
{
}

char AbstractCharSequence::operator[](int32_t index) const
{
    return at(index);
}

std::ostream& AbstractCharSequence::print(std::ostream& out) const
{
    return print(out,0,size());
}

std::ostream& AbstractCharSequence::print(std::ostream& out,int32_t beg, int32_t end) const
{
    while(beg<end)
    {
        out << at(beg);
        beg++;
    }
    return out;
}

std::ostream& AbstractCharSequence::print(std::ostream& out,int32_t beg) const
{
    return print(out,beg,size());
}

std::auto_ptr<std::string> AbstractCharSequence::toString() const
{
    std::ostringstream os;
    print(os);
    return std::auto_ptr<std::string>(new std::string(os.str()));
}
