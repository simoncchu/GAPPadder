
extern int READ_LENGTH;
extern double MEAN_INSERT_SIZE;//mean insert size 
extern double SD_INSERT_SIZE;//standard derivation of insert size
extern double READ_FULL_MAPPED_CUTOFF;//larger than this cutoff, then a read will be considered as fully mapped

const int READ_TYPE_FULLMAP=1;
const int READ_TYPE_CLIP=2;
const int READ_TYPE_LEFT_SOFTCLIP=21;
const int READ_TYPE_RIGHT_SOFTCLIP=22;
const int READ_TYPE_BOTH_SOFTCLIP=23;
const int READ_TYPE_LEFT_HARDCLIP=24;
const int READ_TYPE_RIGHT_HARDCLIP=25;
const int READ_TYPE_BOTH_HARDCLIP=26;
const int READ_TYPE_UNMAP=3;
const int READ_TYPE_OTHER=4;

const int READ_PAIR_MAP_TYPE_11=11;//both mapped 
const int READ_PAIR_MAP_TYPE_10=10;//read map, mate unmap
const int READ_PAIR_MAP_TYPE_01=12;//read unmap, mate map
const int READ_PAIR_MAP_TYPE_00=0;//both unmapped 