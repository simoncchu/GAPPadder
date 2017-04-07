import sys
import os
from subprocess import *
from multiprocessing import Pool
from Bio import SeqIO

with open(sys.argv[1]) as fin_sam:
    pre_rname=""
    cigar=""
    ref=""
    pre_pos=0
    m_pos={}
    for line in fin_sam:
        fields=line.split()
        rname=fields[0]
        ref=fields[2]
        cigar=fields[5]
        if cigar=="*":
            continue
        pos=int(fields[3])
        rname_fields=rname.split("_")
        rname_id=rname_fields[0]+"_"+rname_fields[1]

        ref_fieds=ref.split("_")
        ref_id=ref_fieds[1]

        if rname_fields[0]==ref_id and rname_id==pre_rname:
            if m_pos.has_key(rname_fields[0])==False:
                m_pos[rname_fields[0]]={}
            m_pos[rname_fields[0]][rname_fields[1]]=(pre_pos, pos)

        pre_rname=rname_id
        pre_pos=pos

    with open("only_seq.fa","w") as fout_seq:
        for record in SeqIO.parse(sys.argv[2], "fasta"):

            scf_fields=str(record.id).split("_")
            scf_id=scf_fields[1]
            seq=str(record.seq)
            if m_pos.has_key(scf_id):
                for key in m_pos[scf_id]:
                    start=m_pos[scf_id][key][0]+295
                    end=m_pos[scf_id][key][1]

                    if start>end:
                        print "error!", scf_id, key, start, end
                        continue

                    fout_seq.write(">"+scf_id+"_"+key+"\n")
                    fout_seq.write(seq[start:end+1]+"\n")

