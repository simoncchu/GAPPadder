#ifndef _CORE_BAMALIGNMENTRECORD_H_
#define _CORE_BAMALIGNMENTRECORD_H_

#include<cstdio>
#include<string>
#include<vector>
#include<utility>


/**
.Class.BamAlignmentRecord
..cat:BAM I/O
..summary:Represent a record from a BAM/SAM file.
..remarks:While also used to represent SAM records, called $BamAlignmentRecord$ since the data directly reflects a BAM record (0-based positions, identify references by ids, not names, tags stored in BAM format.)
..include:seqan/bam_io.h
..see:Enum.BamFlags

.Memfunc.BamAlignmentRecord#BamAlignmentRecord
..class:Class.BamAlignmentRecord
..summary:Constructor.
..signature:BamAlignmentRecord()
..remarks:Only the default constructor is provided.

.Memvar.BamAlignmentRecord#INVALID_POS
..class:Class.BamAlignmentRecord
..summary:Static member with invalid/sentinel position value.
..type:nolink:$__uint32$

.Memvar.BamAlignmentRecord#INVALID_REFID
..class:Class.BamAlignmentRecord
..summary:Static member with invalid/sentinel reference id (-1 as in BAM/SAM).
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#INVALID_LEN
..class:Class.BamAlignmentRecord
..summary:Static member with invalid/sentinel position value.
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#qName
..class:Class.BamAlignmentRecord
..summary:The read/query name.
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#flag
..class:Class.BamAlignmentRecord
..summary:The flag of this mapping, see @Enum.BamFlags@ for flag constants and the $hasFlag*$ functions.
..type:nolink:$__uint16$

.Memvar.BamAlignmentRecord#rId
..class:Class.BamAlignmentRecord
..summary:ID of reference for this fragment mapping (0-based, $INVALID_REFID$ for '*').
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#pos
..class:Class.BamAlignmentRecord
..summary:The position of this fragment mapping (0-based, $INVALID_POS$ for '*').
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#mapQ
..class:Class.BamAlignmentRecord
..summary:The mapping quality (255 for '*').
..type:nolink:$__uint8$

.Memvar.BamAlignmentRecord#bin
..class:Class.BamAlignmentRecord
..summary:The bin of the alignment, automatically computed when writing BAM.
..type:nolink:$__uint16$

.Memvar.BamAlignmentRecord#cigar
..class:Class.BamAlignmentRecord
..summary:The CIGAR string as string of @Class.CigarElement@ objects (empty for '*').
..type:nolink:$String<CigarElement<> >$

.Memvar.BamAlignmentRecord#rNextId
..class:Class.BamAlignmentRecord
..summary:ID of reference for next fragment mapping (0-based, $INVALID_REFID$ for '*')
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#pNext
..class:Class.BamAlignmentRecord
..summary:Position of next fragment mapping (0-based, $INVALID_POS$ for '*')
..type:nolink:$__uint32$

.Memvar.BamAlignmentRecord#tLen
..class:Class.BamAlignmentRecord
..summary:The inferred template size ($INVALID_LEN$ for '*')
..type:nolink:$__int32$

.Memvar.BamAlignmentRecord#seq
..class:Class.BamAlignmentRecord
..summary:The sequence string (empty for '*').
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#qual
..class:Class.BamAlignmentRecord
..summary:String with Phred scores (as in SAM file, empty for '*').
..type:Shortcut.CharString

.Memvar.BamAlignmentRecord#tags
..class:Class.BamAlignmentRecord
..summary:Raw BAM tag string, use @Class.BamTagsDict@ for comfortable access.
..type:Shortcut.CharString
*/

class BamAlignmentRecord
{
public:
    std::string qName;          // QNAME
    int flag;                   // FLAG
    std::string rName;          // REF SEQ NAME
	int rID;                    // REF SEQ ID
    int pos;                    // POS
    int mapQ;                   // MAPQ mapping quality, 255 for */invalid
    int bin;                    // bin for indexing 
    std::vector< std::pair<std::string,int> > cigar;  // CIGAR string 
	int rNextID;                // Ref ID of the next segment 
	std::string rNext;          // RNEXT (0-based) 
    int pNext;                  // PNEXT (0-based) 
    int tLen;                   // TLEN
    std::string seq;            // SEQ, as in SAM/BAM file.
    std::string qual;           // Quality string as in SAM (Phred).
    std::string tags;           // Tags, raw as in BAM.

    // Constants for marking pos, reference id and length members invalid (== *).
 /*   static int const INVALID_POS = MaxValue<int>::VALUE;
    static int const INVALID_REFID = -1;
    static int const INVALID_LEN = MaxValue<int>::VALUE;*/
	BamAlignmentRecord(){}

};

//class BamAlignment
//{
//public:
//	std::vector<BamAlignmentRecord> records;
//	
//};

#endif