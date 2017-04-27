DRAFT_GENOME="draft_genome.fa"
ALIGNMENT_LIST=[]
RAW_READS_LIST=[]
BWA_PATH="bwa"
SAMTOOLS_PATH="samtools"
REFINER_PATH="./TERefiner_1"
MERGER_PATH="./ContigsMerger"
KMC_PATH=""
VELVET_PATH=""
THREADS=15
OUTPUT_FOLDER="./GAPPadder_Output/"
VERBOSE=0
MIN_GAP_LENGTH=100
FLANK_LENGTH=300
KMER_LIST=[]

def set_software_paths(pbwa, psamtools, prefiner, pmerger, pkmc, pvelvet):
    if pkmc[-1]!="/":
        pkmc=pkmc+"/"
    if pvelvet[-1]!="/":
        pvelvet=pvelvet+"/"
    global  BWA_PATH
    BWA_PATH=pbwa
    global  SAMTOOLS_PATH
    SAMTOOLS_PATH=psamtools
    global REFINER_PATH
    REFINER_PATH=prefiner
    global MERGER_PATH
    MERGER_PATH=pmerger
    global  KMC_PATH
    KMC_PATH=pkmc
    global VELVET_PATH
    VELVET_PATH=pvelvet

def set_parameters(pdraft, t, pout, vbs, min_gap_len, flank_len):
    global DRAFT_GENOME
    DRAFT_GENOME=pdraft
    global THREADS
    THREADS=t
    if pout[-1]!="/":
        pout=pout+"/"
    global OUTPUT_FOLDER
    OUTPUT_FOLDER=pout
    global VERBOSE
    VERBOSE=vbs
    global MIN_GAP_LENGTH
    MIN_GAP_LENGTH=min_gap_len
    global FLANK_LENGTH
    FLANK_LENGTH=flank_len

def set_alignment_list(algnmt_list):
    global ALIGNMENT_LIST
    for key in algnmt_list:
        ALIGNMENT_LIST.append(key)

def set_raw_reads_list(lraw_reads):
    global RAW_READS_LIST
    for key in lraw_reads:
        RAW_READS_LIST.append(key)

def set_kmer_list(kmer_list):
    global KMER_LIST
    for key in kmer_list:
        KMER_LIST.append(key)

def get_draft_genome():
    global DRAFT_GENOME
    return DRAFT_GENOME

def get_alignment_list():
    global ALIGNMENT_LIST
    return ALIGNMENT_LIST

def get_raw_reads_list():
    global RAW_READS_LIST
    return RAW_READS_LIST

def get_min_gap_length():
    global MIN_GAP_LENGTH
    return MIN_GAP_LENGTH

def get_flank_length():
    global FLANK_LENGTH
    return FLANK_LENGTH

def get_kmer_list():
    global KMER_LIST
    return KMER_LIST

def get_bwa_path():
    global BWA_PATH
    return BWA_PATH

def get_samtools_path():
    global SAMTOOLS_PATH
    return SAMTOOLS_PATH

def get_refiner_path():
    global REFINER_PATH
    return REFINER_PATH

def get_merger_path():
    global MERGER_PATH
    return MERGER_PATH

def get_kmc_path():
    global KMC_PATH
    return KMC_PATH

def get_velvet_path():
    global VELVET_PATH
    return VELVET_PATH

def get_threads_num():
    global THREADS
    return THREADS

def get_output_folder():
    global OUTPUT_FOLDER
    return OUTPUT_FOLDER

def get_verbose():
    global VERBOSE
    return VERBOSE

def print_command(cmd):
    if VERBOSE!=0:
        print "Running command: "+cmd+" ..."