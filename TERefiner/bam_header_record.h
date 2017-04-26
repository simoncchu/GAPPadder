#ifndef _CORE_BAMHEADERRECORD_H_
#define _CORE_BAMHEADERRECORD_H_

#include<cstdio>
#include<string>
#include<vector>
#include<utility>


/**
.Enum.BamHeaderRecordType
..cat:BAM I/O
..summary:Enumeration for the header record type.
..signature:BamHeaderRecordType
..value.BAM_HEADER_FIRST:@Class.BamHeaderRecord@ is of type $@HD$
..value.BAM_HEADER_REFERENCE:@Class.BamHeaderRecord@ is of type $@SQ$
..value.BAM_HEADER_READ_GROUP:@Class.BamHeaderRecord@ is of type $@RG$
..value.BAM_HEADER_PROGRAM:@Class.BamHeaderRecord@ is of type $@PG$
..value.BAM_HEADER_COMMENT:@Class.BamHeaderRecord@ is of type $@CO$    
..include:seqan/bam_io.h    
*/

enum BamHeaderRecordType
{
    BAM_HEADER_FIRST       = 0,
    BAM_HEADER_REFERENCE   = 1,
    BAM_HEADER_READ_GROUP  = 2,
    BAM_HEADER_PROGRAM     = 3,
    BAM_HEADER_COMMENT     = 4
};

/**
.Enum.BamSortOrder
..cat:BAM I/O
..summary:Enumeration for the header record type.
..signature:BamSortOrder
..value.BAM_SORT_UNKNOWN:BAM file sort order is unknown.
..value.BAM_SORT_UNSORTED:BAM file is unsorted.
..value.BAM_SORT_QUERYNAME:BAM file is sorted by query name.
..value.BAM_SORT_COORDINATE:BAM file is sorted by coordinates.
..include:seqan/bam_io.h
*/

enum BamSortOrder
{
    BAM_SORT_UNKNOWN    = 0,
    BAM_SORT_UNSORTED   = 1,
    BAM_SORT_QUERYNAME  = 2,
    BAM_SORT_COORDINATE = 3
};

/**
.Class.BamHeaderRecord
..cat:BAM I/O
..summary:Represents a header entry in a SAM file or the header section of the BAM header.
..signature:BamHeaderRecord
..remarks:Comment records are stored with one tag where the key is empty and the value is the comment. 
..include:seqan/bam_io.h

.Memfunc.BamHeaderRecord#BamHeaderRecord
..class:Class.BamHeaderRecord
..signature:BamHeaderRecord()
..summary:Constructor.
..remarks:Only the default constructor is provided.

.Typedef.BamHeaderRecord#TTagName
..class:Class.BamHeaderRecord
..summary:Type of the tag keys.

.Typedef.BamHeaderRecord#TTagValue
..class:Class.BamHeaderRecord
..summary:Type of the tag values.

.Typedef.BamHeaderRecord#TTag
..class:Class.BamHeaderRecord
..summary:@Class.Pair@ to use for storing tags.

.Typedef.BamHeaderRecord#TTags
..class:Class.BamHeaderRecord
..summary:Type of the string of tag @Class.Pair|Pairs@.

.Memvar.BamHeaderRecord#type
..summary:Type of the record.
..class:Class.BamHeaderRecord
..type:Enum.BamHeaderRecordType

.Memvar.BamHeaderRecord#tags
..summary:The header record's tags.
..class:Class.BamHeaderRecord
..type:Typedef.BamHeaderRecord#TTags 
*/

class BamHeaderRecord
{
public:
    typedef std::string TTagName;
    typedef std::string TTagValue;
    typedef std::pair<TTagName, TTagValue> TTag;
    typedef std::vector<TTag> TTags;

    BamHeaderRecordType type;
    std::vector<std::pair<TTagName, TTagValue> > tags;

    BamHeaderRecord() {}
};

/**
.Class.BamHeader
..cat:BAM I/O
..summary:Stores the information of the BAM header.
..signature:BamHeader
..see:Class.BamHeaderRecord
..include:seqan/bam_io.h

.Memvar.BamHeader#sequenceInfos
..class:Class.BamHeader
..summary:String of $(seqid, length)$ with reference name / length information.
..type:nolink:$String<Pair<CharString, __int32> >$

.Memvar.BamHeader#records
..class:Class.BamHeader    
..summary:String of @Class.BamHeaderRecord|BamHeaderRecords@.
..type:nolink:$String<BamHeaderRecord>$    
*/

class BamHeader
{
public:
    typedef std::pair<std::string, int> TSequenceInfo;
    
    std::vector<std::pair<std::string, int> > sequenceInfos;
    std::vector<BamHeaderRecord> records;
};

//} 

#endif