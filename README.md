# GAPPadder
@2017 by Chong Chu, Xin Li and Yufeng Wu. This software is provided ``as is” without warranty of any
kind. In no event shall the author be held responsible for any damage resulting from the
use of this software. The program package, including source codes, executables, and this
documentation, is distributed free of charge.
If you use this program in a publication, please cite the following reference:  
Chong Chu, Xin Li, and Yufeng Wu. "GAPPadder: A Sensitive Approach for Closing Gaps on Draft Genomes with Short Sequence Reads." bioRxiv (2017): 125534.

## **Functionalities and Usage of GAPPadder**
GAPPadder is designed for closing gaps on the draft genomes with paired-end reads or mate-paired reads. 
The main advantages of GAPPadder is that (Refer to the paper for more detailed information): 
- It collects more reads, especially the repeats associated reads, to close the gaps. 
- It performs a two stage local assembly to construct the gap sequences. 
- It fully utilize the different insert size PE or MP reads. 

## **Dependencies**
The current released version of GAPPadder runs on Linux OS. And GAPPadder needs the following tools to be installed in the machine you are working on.

- Python 2.7 or higher version is required to run GAPPadder. Also need the Biopython package (http://biopython.org/wiki/Download) to be installed. 
- A k-mer counting tool. GAPPadder uses KMC program for performing k-mer counting. KMC can be downloaded from https://github.com/refresh-bio/KMC.
- A reads assembler. GAPPadder uses Velvet at this point. In the future, we may support different assembler. Velvet can be downloaded from: https://www.ebi.ac.uk/~zerbino/velvet/. 
Caution: if you want to assemble k-mers that are longer than 30 bp, you need to recompile Velvet to let it work with longer sequence length. 
For example, use the follow command to recompile, while will make velvet work for k-mer length up to 60. 
	```sh
	$ make clean
	$ make 'MAXKMERLENGTH=60'
	 ```
- Reads mapping. GAPPadder uses bwa mem. BWA (version 0.7 or later) can be downloaded from https://github.com/lh3/bwa.
- Sequence processing utilities. These include the commonly used samtools (v1.3.1 or later, note old version of samtools have different parameter setting with the latest one, now GAPPadder fully support the latest version (1.3.1) of samtools, but may fail when use older version). Our code also uses bamtools (https://github.com/pezmaster31/bamtools), but bamtools is not required to be installed.

## **Download and Install**
First, download the whole folder from https://github.com/Reedwarbler/GAPPadder, including the subfolder TERefiner and ContigsMerger-v0.2.0.

By default, users can directly run the tool and there is no need to install if you have all the dependencies installed. Before run, need to run the following command:
```sh
$ chmod +x ./TERefiner_1  &&  chmod +x ./ContigsMerger 
```

However, on some machines users may fail to run the pre-compiled tools TERefiner_1 and ContigsMerger, then users need to compile by themselves (Note, TERefiner needs bamtools to compile, and users need to set the bamtools path in the makefile) and run the follow commands:
```sh
$ cd TERefiner  &&  make  &&  cd .. 
$ cd ./ContigsCompactor-v0.2.0/ContigsMerger/  &&  make  &&  cd .. 
$ cp ./TERefiner/TERefiner_1 ./  &&  cp ./ContigsCompactor-v0.2.0/ContigsMerger/ContigsMerger ./
$ chmod +x ./TERefiner_1  &&  chmod +x ./ContigsMerger 
```
## **Preparing inputs**

GAPPadder needs a configuration file in JSON format. The configuration file tells GAPPadder the basic settings. Users can find one sample from the same folder in this github cite.
Once finish the configuration file, users can use this website (http://jsonlint.com/) to check whether there are errors. 
Here, we give an explanation on the parameters.

**draft_genome**

The path of the draft genome

**raw_reads**

The groups of paired end reads, with each pair one group.

**alignments**

The path of the alignment files (must be sorted bam/cram files). Should keep the same number as the group of PE reads.

**software_path**

- kmc: Set to the folder of KMC executable. If already added the PATH, then use "" in this field.
- velvet: Set to the folder of Velvet executable. If already added the PATH, then use "" in this field.
- bwa: Set to the path of BWA executable. If already added the PATH, then use "bwa" in this field.
- samtools: Set to the path of samtools executable. If already added the PATH, then use "samtools" in this field.
- TERefiner: set to the path of TERefiner_1 executable. If put at the same folder as main.py, then keep the same as the sample configuration file.
- ContigsMerger: set to the path of ContigsMerger executable. If put at the same folder as main.py, then keep the same as the sample configuration file.

**parameters**

- working_folder : Where to output the results.
- min_gap_size: Minimum gap size the tool focues on.
- flank_length: Length of the left and right flank regions.
- nthreads: Number of threads.
- verbose. If set to be 1, output more information about the current running states of GAPPadder.

**kmer_length**

The kmer lengths for assembly. For each k, there are several sub-k, which are the length of kmers going to be used by velvet.

## **Basic usage**

Preprocess the draft genome to get the gap positions and flank regions:
```sh
$ python ./main.py -c Preprocess -g configuration-file-name 
```

Collect reads for each gap:
```sh
$ python ./main.py -c Collect -g configuration-file-name 
```

Construct the gap sequence and pick the best one:
```sh
$ python ./main.py -c Assembly -g configuration-file-name
```

Clean the old data:
```sh
$ python ./main.py -c Clean -g configuration-file-name
```

## **Ouptput**

picked_seqs.fa contains the selected fully closed and extended gap sequences.
 