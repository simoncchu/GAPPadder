import sys
from subprocess import *
from multiprocessing.dummy import Pool as ThreadPool

sf_bam=""
sf_fai=""
sf_gap_pos=""
anchor_mapq=""

def run_cmd(cmd1):
    print cmd1
    check_output(cmd1, shell=True)

def dispath_collect_jobs(n):
    global sf_fai
    global sf_bam
    global sf_gap_pos
    global anchor_mapq

    m_scaffold_has_gap_list={}
    with open(sf_gap_pos) as fin_gap_pos:
        for line in fin_gap_pos:
            fields=line.split()
            m_scaffold_has_gap_list[fields[3]]=1

    cmd_list=[]
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            if m_scaffold_has_gap_list.has_key(fields[0])==False:
                continue
            cmd="samtools view {0} \"{1}\" | python collect_reads_for_gaps_short_insert.py {2} {3} -".format(sf_bam, fields[0], sf_gap_pos, anchor_mapq)
            cmd_list.append(cmd)

    pool = ThreadPool(n)
    pool.map(run_cmd, cmd_list)
    pool.close()
    pool.join()


if __name__ == "__main__":
    sf_fai=sys.argv[1]
    sf_bam=sys.argv[2]
    sf_gap_pos=sys.argv[3]
    anchor_mapq=int(sys.argv[4])

    dispath_collect_jobs(int(sys.argv[5]))
