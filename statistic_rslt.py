import sys
import os
from subprocess import *
from multiprocessing import Pool
from Bio import SeqIO


def run_statistic(id):
    #sf_picked="./picked/{0}.fa".format(id)
    #sf_true_seq="./only_gap_seqs/{0}.fa".format(id)

    sf_true_seq="./picked/{0}.fa".format(id)
    sf_picked="./../only_gap_seqs/{0}.fa".format(id)

    if os.path.exists(sf_picked)==False or os.path.exists(sf_true_seq)==False:
        return

    cmd="bwa index {0}".format(sf_picked)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="bwa mem {0} {1} | samtools view -h -b -S - | samtools sort - -o ./alignment/{2}.sort.bam".format(sf_picked,sf_true_seq,id)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="samtools view ./alignment/{0}.sort.bam > ./alignment/{1}.sam".format(id,id)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def do_statistic(sf_fai, sf_gap_pos, n, sf_picked):
    m_scaffold_id={}
    m_gap_length={}

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
            m_gap_length[sf_fq]=fields[2] ####################################save gap length
            cnt=cnt+1
            pre_id=scaffold_id

            fq_list.append(sf_fq)

    #return m_gap_length #########################################################################################
    for record in SeqIO.parse(sf_picked, "fasta"):
        s_id=str(record.id)
        fields=s_id.split("_")
        s_new_id=fields[0]+"_"+fields[1]
        s_p="./picked/{0}.fa".format(s_new_id)
        with open(s_p,"w") as fout_picked:
            fout_picked.write(">"+s_id+"\n")
            fout_picked.write(str(record.seq)+"\n")

    pool = Pool(n)
    pool.map(run_statistic, fq_list)
    pool.close()
    pool.join()

    cmd="cat ./alignment/*.sam > hit.sam"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    return m_gap_length

def is_qualified_clipped(cigar, cutoff_len):
    l=len(cigar)
    signal=[]
    lenth=[]
    temp=""
    for i in range(l):
        if cigar[i]>="0" and cigar[i]<="9":
            temp=temp+cigar[i]
        else:
            signal.append(cigar[i])
            lenth.append(int(temp))
            temp=""
    if (signal[0]=="S" or signal[0]=="H") and lenth[0]>=cutoff_len:
        return True
    if (signal[len(signal)-1]=="S" or signal[len(signal)-1]=="H") and lenth[len(signal)-1]>=cutoff_len:
        return True
    return False


def cnt(sf_algn, m_gap_length):
    cnt=0
    m_id={}
    m_hit={}
    with open(sf_algn) as fin_algn:
        for line in fin_algn:
            fields=line.split()
            rname=fields[0]
            m_id[rname]=1
            cigar=fields[5]
            if cigar=="*":
                continue

            if is_qualified_clipped(cigar,20) ==False:
                m_hit[fields[2]]=0
                #cnt=cnt+1
                #print rname, fields[2]
    print len(m_hit),len(m_id)

    with open("hit_list.txt", "w") as fout_hit_list:
        for id in m_hit:
            fout_hit_list.write(id+"\n")

    with open("closed_gap_length.txt", "w") as fout_gap_lenth:
        fout_gap_lenth.write("[") #################################################
        for id in m_hit:
            fields=id.split("_")
            fout_gap_lenth.write(m_gap_length[fields[0]+"_"+fields[1]]+",")
        fout_gap_lenth.write("0]") #################################################

def output_hit_list(sf_algn):
    cnt=0
    m_id={}
    m_hit={}
    with open(sf_algn) as fin_algn:
        for line in fin_algn:
            fields=line.split()
            rname=fields[0]
            m_id[rname]=1
            cigar=fields[5]
            if cigar=="*":
                continue

            if is_qualified_clipped(cigar,20) ==False:
                m_hit[fields[2]]=0
                #cnt=cnt+1
                #print rname, fields[2]
    #print len(m_hit),len(m_id)

    with open("hit_list.txt", "w") as fout_hit_list:
        for id in m_hit:
            fout_hit_list.write(id+"\n")



if __name__ == "__main__":
    # sf_fai=sys.argv[1]
    # sf_pos=sys.argv[2]
    # n_jobs=int(sys.argv[3])
    # sf_picked=sys.argv[4]
    #
    # sf_temp="./alignment"
    # if os.path.exists(sf_temp)==False:
    #     cmd="mkdir {0}".format(sf_temp)
    #     Popen(cmd, shell = True, stdout = PIPE).communicate()
    # sf_temp="./picked"
    # if os.path.exists(sf_temp)==False:
    #     cmd="mkdir {0}".format(sf_temp)
    #     Popen(cmd, shell = True, stdout = PIPE).communicate()

    #m_gap_length=do_statistic(sf_fai, sf_pos, n_jobs, sf_picked)
    #cnt(sys.argv[1])
    #cnt("hit.sam", m_gap_length)

    output_hit_list(sys.argv[1])
