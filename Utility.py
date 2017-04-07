BWA_PATH="bwa"
SAMTOOLS_PATH="samtools"
REFINER_PATH="./TERefiner_1"
JELLYFISH_PATH=""
VELVET_PATH=""
THREADS=15
OUTPUT_FOLDER="./REPdenovo_Output/"
VERBOSE=True

def set_parameters(pbwa,psamtools,prefiner,pjellyfish,pvelvet,t,pout,vbs):
    global  BWA_PATH
    BWA_PATH=pbwa
    global  SAMTOOLS_PATH
    SAMTOOLS_PATH=psamtools
    global REFINER_PATH
    REFINER_PATH=prefiner
    global  JELLYFISH_PATH
    JELLYFISH_PATH=pjellyfish
    global VELVET_PATH
    VELVET_PATH=pvelvet
    global THREADS
    THREADS=t
    global OUTPUT_FOLDER
    OUTPUT_FOLDER=pout
    global VERBOSE
    VERBOSE=vbs

def get_bwa_path():
    global BWA_PATH
    return BWA_PATH

def get_samtools_path():
    global SAMTOOLS_PATH
    return SAMTOOLS_PATH

def get_refiner_path():
    global REFINER_PATH
    return REFINER_PATH

def get_jellyfish_path():
    global JELLYFISH_PATH
    return JELLYFISH_PATH

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