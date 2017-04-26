CC = g++
CFLAGS =  -O3 -Wall -mcmodel=medium


all: ContigsMerger 

abstractcharsequence.o: abstractcharsequence.cpp abstractcharsequence.h 
	$(CC)  $(CFLAGS) -c abstractcharsequence.cpp

fastareader.o: fastareader.cpp fastareader.h
	$(CC)  $(CFLAGS) -c fastareader.cpp

fastaMultiSeqs.o: fastaMultiSeqs.cpp fastaMultiSeqs.h
	$(CC)  $(CFLAGS) -c fastaMultiSeqs.cpp

MurmurHash3.o: MurmurHash3.cpp MurmurHash3.h
	$(CC)  $(CFLAGS) -c MurmurHash3.cpp

ContigsCompactor.o: ContigsCompactor.cpp ContigsCompactor.h
	$(CC)  $(CFLAGS) -c ContigsCompactor.cpp
    
GenSeqsUtils.o: GenSeqsUtils.cpp GenSeqsUtils.h
	$(CC)  $(CFLAGS) -c GenSeqsUtils.cpp
    
KmerUtils.o: KmerUtils.cpp KmerUtils.h
	$(CC)  $(CFLAGS) -c KmerUtils.cpp

Utils-basic.o: Utils-basic.cpp Utils-basic.h 
	$(CC)  $(CFLAGS) -c Utils-basic.cpp
    
GraphUtils.o: GraphUtils.cpp GraphUtils.h 
	$(CC)  $(CFLAGS) -c GraphUtils.cpp  

ScaffoldUtils.o: ScaffoldUtils.cpp ScaffoldUtils.h
	$(CC)  $(CFLAGS) -c ScaffoldUtils.cpp


main.o: main.cpp
	$(CC)  $(CFLAGS) -c main.cpp
    

# Now final link

ContigsMerger:  abstractcharsequence.o main.o fastareader.o fastaMultiSeqs.o  GenSeqsUtils.o MurmurHash3.o ContigsCompactor.o KmerUtils.o \
  Utils-basic.o GraphUtils.o ScaffoldUtils.o
	$(CC) $(CFLAGS) -o ContigsMerger abstractcharsequence.o main.o fastareader.o fastaMultiSeqs.o  GenSeqsUtils.o MurmurHash3.o ContigsCompactor.o \
         KmerUtils.o Utils-basic.o GraphUtils.o ScaffoldUtils.o \
         -lz -lm -lpthread
