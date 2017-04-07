import sys
import os
import user
from Bio import SeqIO
from subprocess import *

def gnrt_gap_positions(ref_path, sf_pos):
    with open(sf_pos,"w") as fout_pos:
        for record in SeqIO.parse(ref_path, "fasta"):
            seq1=record.seq
            #search for N and NX
            pos=0
            while True:
                start=seq1.find("N",pos)
                if start==-1:
                    break
                pos=start+1
                end1=seq1.find("A", pos)
                end2=seq1.find("C", pos)
                end3=seq1.find("G", pos)
                end4=seq1.find("T", pos)

                min_pos=-1
                if end1!=-1:
                    min_pos=end1
                    if end2!=-1 and end2<min_pos:
                        min_pos=end2
                    if end3!=-1 and end3<min_pos:
                        min_pos=end3
                    if end4!=-1 and end4<min_pos:
                        min_pos=end4
                elif end2!=-1:
                    min_pos=end2
                    if end3!=-1 and end3<min_pos:
                        min_pos=end3
                    if end4!=-1 and end4<min_pos:
                        min_pos=end4
                elif end3!=-1:
                    min_pos=end3
                    if end4!=-1 and end4<min_pos:
                        min_pos=end4
                elif end4!=-1:
                    min_pos=end4

                if min_pos==-1:
                    break

                #print start, min_pos, min_pos-start, str(record.id)
                str1="{0} {1} {2} {3}\n".format(start, min_pos, min_pos-start, str(record.id))
                fout_pos.write(str1)
                pos=min_pos+2

def get_gap_flank_seqs(ref_path, sf_gap_pos, frame_length, sf_fai):
    m_scaffold_id={}
    l_scaffold_id=[]
    cnt=0##scaffold id, started from 0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1
                l_scaffold_id.append(scaffold_id)

    num=1
    for record in SeqIO.parse(ref_path, "fasta"):
        seq1=record.seq
        with open(sf_gap_pos) as fin_gap_pos:
            for line in fin_gap_pos:
                fields=line.split()
                scf_id=fields[3]
                if scf_id!=str(record.id):
                    num=1
                    continue

                cnt=m_scaffold_id[scf_id]
                gap_id="{0}_{1}".format(cnt,num)
                sf_flank="flank_regions/{0}.fa".format(gap_id)

                with open(sf_flank,"w") as fout_flank:
                    start=int(fields[0])
                    end=int(fields[1])
                    fout_flank.write(">"+gap_id+"_left\n")
                    if start<frame_length:
                        fout_flank.write(str(seq1[0:start-5])+"\n")
                    else:
                        fout_flank.write(str(seq1[start-frame_length:start-5])+"\n")
                    fout_flank.write(">"+gap_id+"_right\n")
                    fout_flank.write(str(seq1[end+5:end+frame_length])+"\n")
                    num=num+1

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

def gnrt_gap_seqs(sf_sam, sf_ref, flank_length):
    ##first cat all the seqs, and align to true genome, sort and covert to sam
    pre=""
    pre_pos=-1
    pre_ref=""
    pre_ori=""

    gap_seq_pos={}
    with open(sf_sam) as fin_sam:
        for line in fin_sam:
            fields=line.split()
            pos=int(fields[3])
            ref=fields[2]
            cigar=fields[5]

            sname_fields=fields[0].split("_")
            scf_id=sname_fields[0]
            scf_gap_id=sname_fields[1]
            scf_ori=sname_fields[2]
            sid=scf_id+"_"+scf_gap_id

            if cigar=="*" or is_qualified_clipped(cigar,10)==True:#unmapped
                pre=sid
                continue

            if pre != sid or ref != pre_ref or scf_ori==pre_ori:
                pre=sid
                pre_pos=pos
                pre_ref=ref
                pre_ori=scf_ori
                continue

            # start=pre_pos+flank_length-1
            # end=pos

            start=pre_pos+flank_length-1
            end=pos

            #print start,end            ################################################################################

            if start<end:
                if gap_seq_pos.has_key(ref)==False:
                    gap_seq_pos[ref]={}
                gap_seq_pos[ref][str(start)+"_"+str(end)]=sid
                #print sid, start, end #################################################################################
            pre=sid
            pre_pos=pos
            pre_ref=ref
            pre_ori=scf_ori

    for record in SeqIO.parse(sf_ref, "fasta"):
        scaffold_id=str(record.id)
        if gap_seq_pos.has_key(scaffold_id)==False:
            continue

        for s_pos in gap_seq_pos[scaffold_id]:
            gap_id=gap_seq_pos[scaffold_id][s_pos]
            fields=s_pos.split("_")

            start=int(fields[0])
            end=int(fields[1])

            seq=str(record.seq[start:end])
            sf_gap_seq="only_gap_seqs/{0}.fa".format(gap_id)
            #sf_gap_seq="gap_seqs/{0}.fa".format(gap_id)
            with open(sf_gap_seq,"w") as fout_gap_seq:
                fout_gap_seq.write(">"+gap_id+"\n")
                fout_gap_seq.write(seq+"\n")


def rm_bias(sf_gap_pos):
    with open(sf_gap_pos) as fin_gap_pos:#for each gap
        for line in fin_gap_pos:
            fields=line.split()
            id=fields[0]

            sf_old="./gap_seqs/{0}.fa".format(id)
            if os.path.exists(sf_old)==False:
                continue

            sf_new="./only_gap_seqs/{0}.fa".format(id)
            with open(sf_new,"w") as fout_new:
                for record in SeqIO.parse(sf_old, "fasta"):
                    seq_len=len(str(record.seq))
                    if seq_len<(390+105):
                        continue
                    fout_new.write(">"+str(record.id)+"\n")
                    fout_new.write(str(record.seq)[389:seq_len-105])



def test(sf_seq):
    for record in SeqIO.parse(sf_seq, "fasta"):
        #85689374 85690678
        print str(record.seq)[85689174:85690878]


ref_path=sys.argv[1]
#sf_gap_pos=sys.argv[2]
sf_sam=sys.argv[2]
flank_length=int(sys.argv[3])
#sf_fai=sys.argv[4]

#gnrt_gap_positions(ref_path,sf_gap_pos)
#get_gap_flank_seqs(ref_path, sf_gap_pos, flank_length, sf_fai)
gnrt_gap_seqs(sf_sam, ref_path, flank_length)
#rm_bias(sys.argv[1])
