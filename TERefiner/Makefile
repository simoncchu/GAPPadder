BAMTOOLS=../bamtools-master/include
BAMTOOLS_LD=../bamtools-master/lib

CC = g++
CFLAGS =  -O3 -Wall -static -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)

all: TERefiner

public_func.o :	public_func.cpp public_func.h
	$(CC) $(CFLAGS) -c public_func.cpp

StrOperation.o: StrOperation.cpp StrOperation.h
	$(CC) $(CFLAGS) -c StrOperation.cpp

Alignment.o: Alignment.cpp Alignment.h
	$(CC) $(CFLAGS) -c Alignment.cpp
	
bam_parse.o : bam_parse.cpp bam_parse.h
	$(CC) $(CFLAGS) -c bam_parse.cpp

local_alignment.o: ./algorithms/local_alignment.cpp ./algorithms/local_alignment.h
	$(CC) $(CFLAGS) -c ./algorithms/local_alignment.cpp

fai_parser.o : fai_parser.cpp fai_parser.h
	$(CC) $(CFLAGS) -c fai_parser.cpp

contigs.o : contigs.cpp contigs.h
	$(CC) $(CFLAGS) -c contigs.cpp

Coverage.o : Coverage.cpp Coverage.h
	$(CC) $(CFLAGS) -c Coverage.cpp

scaffolding.o: scaffolding.cpp scaffolding.h
	$(CC) $(CFLAGS) -c scaffolding.cpp

RepeatsClassifier.o: RepeatsClassifier.cpp RepeatsClassifier.h
	$(CC) $(CFLAGS) -c RepeatsClassifier.cpp

refiner.o: refiner.cpp refiner.h
	$(CC) $(CFLAGS) -c refiner.cpp

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp

TERefiner: public_func.o StrOperation.o local_alignment.o Alignment.o bam_parse.o fai_parser.o contigs.o Coverage.o scaffolding.o RepeatsClassifier.o refiner.o main.o 
	$(CC) $(CFLAGS) -o TERefiner_1 public_func.o StrOperation.o local_alignment.o Alignment.o bam_parse.o fai_parser.o contigs.o Coverage.o scaffolding.o RepeatsClassifier.o refiner.o main.o \
	-lbamtools -lz -lm