import sys
import os
from subprocess import *
from multiprocessing import Pool

def run_clear(id):
    sf_ori_contig="./velvet_temp/{0}/original_contigs_before_merging.fa".format(id)
    if os.path.exists(sf_ori_contig):
        print id

    sf_contig0="./velvet_temp/{0}/contigs.fa".format(id)
    sf_contig1="./velvet_temp/{0}/contigs_31_29.fa".format(id)
    sf_contig2="./velvet_temp/{0}/contigs_41_37.fa".format(id)
    sf_contig3="./velvet_temp/{0}/contigs_41_39.fa".format(id)
    sf_contig5="./velvet_temp/{0}/contigs_51_47.fa".format(id)
    sf_contig6="./velvet_temp/{0}/contigs_61_57.fa".format(id)
    sf_seq="./velvet_temp/{0}/Sequences".format(id)
    if os.path.exists(sf_contig0) and os.path.exists(sf_contig1) and os.path.exists(sf_contig2) \
            and os.path.exists(sf_contig3) and os.path.exists(sf_contig5) and os.path.exists(sf_contig6)\
            and os.path.exists(sf_seq):
        cmd="rm ./temp/{0}.*".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="rm ./velvet_temp/{0}/Sequences".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm ./velvet_temp/{0}/Roadmaps".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm ./velvet_temp/{0}/PreGraph".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm ./velvet_temp/{0}/Graph".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm ./velvet_temp/{0}/LastGraph".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm ./velvet_temp/{0}/Log".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm ./velvet_temp/{0}/stats.txt".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="rm ./kmers/{0}_*".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()


def run_statistic(sf_fai, sf_gap_pos):
    m_scaffold_id={}
    cnt=0
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            scaffold_id=fields[0]
            m_scaffold_id[scaffold_id]=cnt
            cnt=cnt+1
    cnt=1
    pre_id=0

    with open(sf_gap_pos) as fin_gap_pos:#for each gap
        for line in fin_gap_pos:
            fields=line.split()
            scaffold=fields[3]
            scaffold_id=m_scaffold_id[scaffold]

            if pre_id!=scaffold_id:
                cnt=1

            sf_fq="{0}_{1}".format(scaffold_id,cnt)
            sf_contig0="./velvet_temp/{0}/contigs.fa_no_dup.fa".format(sf_fq)

            if os.path.exists(sf_contig0):
                sz=os.path.getsize(sf_contig0)
                print sf_fq, sz

            cnt=cnt+1
            pre_id=scaffold_id


def run_cat(id):
    cmd="cat ./velvet_temp/{0}/contigs_31_29.fa ./velvet_temp/{0}/contigs_41_37.fa ./velvet_temp/{0}/contigs_41_39.fa" \
        " ./velvet_temp/{0}/contigs_51_47.fa ./velvet_temp/{0}/contigs_61_57.fa > ./velvet_temp/{0}/contigs.fa"\
        .format(id,id,id,id,id,id)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def clear_temp(sf_fai, sf_gap_pos, n):
    m_scaffold_id={}
    cnt=0
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            scaffold_id=fields[0]
            m_scaffold_id[scaffold_id]=cnt
            cnt=cnt+1

    fq_list=[]
    cnt=1
    pre_id=0

    with open(sf_gap_pos) as fin_gap_pos:#for each gap
        for line in fin_gap_pos:
            fields=line.split()
            scaffold=fields[3]
            scaffold_id=m_scaffold_id[scaffold]

            if pre_id!=scaffold_id:
                cnt=1

            sf_fq="{0}_{1}".format(scaffold_id,cnt)
            cnt=cnt+1
            pre_id=scaffold_id

            fq_list.append(sf_fq)

    pool = Pool(n)
    pool.map(run_cat, fq_list)
    pool.close()
    pool.join()

if __name__ == "__main__":
    sf_fai=sys.argv[1]
    sf_pos=sys.argv[2]
    n_jobs=int(sys.argv[3])
    clear_temp(sf_fai, sf_pos, n_jobs)
    #run_statistic(sf_fai, sf_pos)