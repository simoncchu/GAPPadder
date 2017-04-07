import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def put_gap_seq_back_to_scaffold(sf_scf, sf_gap_pos, sf_gap_seq, sf_new_scf):
    m_gap_seq={}
    #sf_picked="picked_seqs.fa"
    for record in SeqIO.parse(sf_gap_seq, "fasta"):
        gap_name=str(record.id)
        gap_seq=str(record.seq)
        gap_name_fields=gap_name.split("_")
        #print int(gap_id_fields[0]),int(gap_id_fields[1])##################################################
        scf_id=int(gap_name_fields[0])
        gap_id=int(gap_name_fields[1])
        if m_gap_seq.has_key(scf_id)==False:
            m_gap_seq[scf_id]={}
        m_gap_seq[scf_id][gap_id]=gap_seq

    m_gap_pos={}
    m_scaffold_id={}
    cnt=0
    sf_fai=sf_scf+".fai"
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
            start=int(fields[0])
            end=int(fields[1])
            scaffold=fields[3]

            scaffold_id=int(m_scaffold_id[scaffold])

            if pre_id!=scaffold_id:
                cnt=1

            #gap_id="{0}_{1}".format(scaffold_id,cnt)
            #m_gap_pos[gap_id]=(start,end)
            if m_gap_pos.has_key(scaffold_id)==False:
                m_gap_pos[scaffold_id]={}
            m_gap_pos[scaffold_id][cnt]=(start,end)

            cnt=cnt+1
            pre_id=scaffold_id

    #print m_gap_pos #################################################################################

    handle=open(sf_new_scf,"w")
    for record in SeqIO.parse(sf_scf, "fasta"):
        scf_name=str(record.id)
        scf_id=int(m_scaffold_id[scf_name]) ##start from 0

        #print m_gap_pos[scf_id]##################################################################3
        if m_gap_pos.has_key(scf_id)==False:
            SeqIO.write(record,handle,"fasta")
            #print "Scaffold ", scf_id, "does not exist!"
            continue

        scf_gap_num=len(m_gap_pos[scf_id])
        if scf_gap_num<=0:
            SeqIO.write(record,handle,"fasta")
            #print "Scaffold ", scf_id, "has no gaps!"
            continue
        if m_gap_seq.has_key(scf_id)==False:
            SeqIO.write(record,handle,"fasta")
            #print "Scaffold ", scf_id, "'s gaps are not filled!"
            continue
        pre_pos=0
        s_final_scf=""
        for i in range(scf_gap_num): #for each gap of this scaffold
            gap_id = i+1 #start from 1
            if m_gap_seq[scf_id].has_key(gap_id)==False:
                print "Gap ", gap_id, "of Scaffold ", scf_id, "does not exist!!!"
                continue

            gap_start = m_gap_pos[scf_id][gap_id][0]
            gap_end = m_gap_pos[scf_id][gap_id][1]

            s_final_scf=s_final_scf+str(record.seq[pre_pos:gap_start])
            s_final_scf=s_final_scf+m_gap_seq[scf_id][gap_id]
            pre_pos=gap_end+1
        s_final_scf=s_final_scf+str(record.seq[pre_pos:])

        new_scf=SeqRecord(Seq(s_final_scf), id = record.id, description="")
        SeqIO.write(new_scf,handle,"fasta")
    handle.close()

if __name__ == "__main__":
    sf_scf=sys.argv[1]
    sf_gap_pos=sys.argv[2]
    sf_gap_seq=sys.argv[3]
    sf_new_scf=sys.argv[4]

    put_gap_seq_back_to_scaffold(sf_scf, sf_gap_pos, sf_gap_seq, sf_new_scf)

