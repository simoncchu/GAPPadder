__author__ = 'Chong Chu'

import sys
import os
from subprocess import *
from Utility import print_command
from Utility import BWA_PATH
from Utility import get_bwa_path
from Utility import SAMTOOLS_PATH
from Utility import get_samtools_path
from Utility import REFINER_PATH
from Utility import get_refiner_path

def remove_duplicate_contained(fcontig, foutput, cutoff, rm_contained):
    BWA_PATH=get_bwa_path()
    SAMTOOLS_PATH=get_samtools_path()
    REFINER_PATH=get_refiner_path()

    #remove duplicate or contained contigs
    cmd="{0} faidx {1}".format(SAMTOOLS_PATH,fcontig)
    #print_command("Running command: "+cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} -U -r {1} -o {2}".format(REFINER_PATH,fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="{0} faidx {1}".format(SAMTOOLS_PATH,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} index {1}".format(BWA_PATH,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} mem -a {1} {2} > {3}.itself.sam".format(BWA_PATH,fcontig,fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="{0} view -h -S -b {1}.itself.sam > {2}.itself.bam".format(SAMTOOLS_PATH,fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} sort {1}.itself.bam -o {2}.itself.sort.bam".format(SAMTOOLS_PATH,fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} index {1}.itself.sort.bam".format(SAMTOOLS_PATH,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    if rm_contained==0:
        cmd="{0} -P -b {1}.itself.sort.bam -r {2} -o {3} -c {4} -g".format(REFINER_PATH,fcontig,fcontig,foutput,cutoff)
    elif rm_contained==1:
        cmd="{0} -P -b {1}.itself.sort.bam -r {2} -o {3} -c {4}".format(REFINER_PATH,fcontig,fcontig,foutput,cutoff)
    else:
        cmd="{0} -K -b {1}.itself.sort.bam -r {2} -o {3} -c {4}".format(REFINER_PATH,fcontig,fcontig,foutput,cutoff)
    print_command(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    ##clean all the temporary files
    cmd="rm {0}.sa".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.pac".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.bwt".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.ann".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.amb".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.sam".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.bam".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.sort.bam".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.sort.bam.bai".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.fai".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def merge_contigs(contigs_merger_path, fout_folder, nthreads, cutoff_dup_bf_merge, cutoff_dup_af_merge):
    fcontig=fout_folder+"contigs.fa"
    foutput=fcontig+"_no_dup.fa"
    remove_duplicate_contained(fcontig, foutput, cutoff_dup_bf_merge, False)

    if os.path.exists(foutput)==True:
        if int(os.path.getsize(foutput))>1000000:
            os.rename(fout_folder+"contigs.fa", fout_folder+"original_contigs_before_merging.fa")
            os.rename(foutput,fout_folder+"contigs.fa")
            return

    cmd="{0} -s 0.4 -i1 -2.0 -i2 -2.0 -x 12 -y 50 -k 10 -t {1} -m 1 -o {2}.merge.info {3} > {4}.merged.fa".format(contigs_merger_path,nthreads,foutput,foutput,foutput)
    # ###cmd="./ContigsMerger_v0.1.7  -s 0.2 -i1 -6.0 -i2 -6.0  -x 15 -k 5 -o {0}.merge.info {1} > {2}.merged.fa".format(foutput,foutput,foutput)
    print_command(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
        #
    fcontig="{0}.merged.fa".format(foutput)
    foutput="{0}.no_dup.fa".format(fcontig)
    remove_duplicate_contained(fcontig,foutput, cutoff_dup_af_merge, False)
    fcontig=foutput
    foutput=fcontig+".no_contained.fa"
    remove_duplicate_contained(fcontig, foutput, cutoff_dup_af_merge, True)
    #
    # #rename contigs.fa
    os.rename(fout_folder+"contigs.fa", fout_folder+"original_contigs_before_merging.fa")
    os.rename(foutput,fout_folder+"contigs.fa")


def rm_dup_contain(fout_folder, cutoff_dup_af_merge):
    fcontig=fout_folder+"contigs.fa"
    foutput="{0}.no_dup.fa".format(fcontig)
    remove_duplicate_contained(fcontig,foutput, cutoff_dup_af_merge, 0)
    fcontig=foutput
    foutput=fcontig+".no_contained.fa"
    remove_duplicate_contained(fcontig, foutput, cutoff_dup_af_merge, 1)

def rm_contain(fout_folder, cutoff_ctn):
    fcontig=fout_folder+"contigs.fa"
    foutput="{0}.no_ctn.fa".format(fcontig)
    remove_duplicate_contained(fcontig,foutput, cutoff_ctn, 2)