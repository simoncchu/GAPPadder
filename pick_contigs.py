import sys
import os
from subprocess import *
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from Bio import SeqIO
from Utility import *

BOTH_CLIP=1
LEFT_CLIP=2
RIGTH_CLIP=3
UNCLIP=4
m_flanks={}
bwa_min_score=30
working_folder=""
bwa_path=""
samtools_path=""

def gnrt_reverse_complementary(s):
    lth=len(s)
    s_rc=""
    for i in range(lth-1,-1,-1):
        if s[i]=="A" or s[i]=="a":
            s_rc=s_rc+"T"
        elif s[i]=="T" or s[i]=="t":
            s_rc=s_rc+"A"
        elif s[i]=="C" or s[i]=="c":
            s_rc=s_rc+"G"
        elif s[i]=="G" or s[i]=="g":
            s_rc=s_rc+"C"
        else:
            s_rc=s_rc+s[i]
    return s_rc

def get_clip_type_length(cigar):
    global LEFT_CLIP
    global RIGTH_CLIP
    global UNCLIP
    global BOTH_CLIP
    l=len(cigar)
    signal=[]
    lenth=[]
    temp=""
    total_M=0
    for i in range(l):
        if cigar[i]>="0" and cigar[i]<="9":
            temp=temp+cigar[i]
        else:
            if cigar[i]=="M":
                total_M=total_M+int(temp)
            signal.append(cigar[i])
            lenth.append(int(temp))
            temp=""

    clip_type=UNCLIP
    if (signal[0]=="S" or signal[0]=="H") and (signal[len(signal)-1]=="S" or signal[len(signal)-1]=="H"):
        clip_type=BOTH_CLIP
    elif signal[0]=="S" or signal[0]=="H":
        clip_type=LEFT_CLIP
    elif signal[len(signal)-1]=="S" or signal[len(signal)-1]=="H":
        clip_type=RIGTH_CLIP
    return clip_type, total_M

def run_pick_full_constructed_contig(id):
    global working_folder
    global LEFT_CLIP
    global RIGTH_CLIP
    global UNCLIP
    global BOTH_CLIP
    global bwa_path
    global samtools_path
    global bwa_min_score

    sf_flank=working_folder+"flank_regions/{0}.fa".format(id)
    if os.path.exists(sf_flank)==False:
        print "Wrong flank regions: ", id
        return

    sf_contig=working_folder+"velvet_temp/{0}/contigs.fa".format(id)
    cmd="{0} index {1}".format(bwa_path, sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    sf_flank_alnmt=working_folder+"velvet_temp/{0}/flanks.sam".format(id)
    cmd="{0} mem -T {1} -a {2} {3} | {4} view -S - > {5}"\
        .format(bwa_path, bwa_min_score, sf_contig, sf_flank, samtools_path, sf_flank_alnmt)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    sf_pk=working_folder+"velvet_temp/{0}/picked_seqs.fa".format(id)
    if os.path.exists(sf_pk):
        cmd="rm {0}".format(sf_pk)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    sf_pk=working_folder+"velvet_temp/{0}/picked_contigs.fa".format(id)
    if os.path.exists(sf_pk):
        cmd="rm {0}".format(sf_pk)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    m_hit_id={}
    m_picked_id={}
    with open(sf_flank_alnmt) as fin_algnmt:
        for line in fin_algnmt:
            fields=line.split()
            cigar=fields[5]
            flag=int(fields[1])
            if cigar=="*":
                continue
            clip_type, map_length=get_clip_type_length(cigar)

            if clip_type==BOTH_CLIP:
                continue

            map_pos=int(fields[3])

            rname=fields[2]
            qname=fields[0]
            qname_field=qname.split("_")
            if qname_field[-1]=="left":
                if m_hit_id.has_key(rname)==False:
                    m_hit_id[rname]={}
                    m_hit_id[rname]["left"]={}
                    m_hit_id[rname]["left"][clip_type]=(flag, map_length, map_pos)
                elif m_hit_id[rname].has_key("left")==False:
                    m_hit_id[rname]["left"]={}
                    m_hit_id[rname]["left"][clip_type]=(flag, map_length, map_pos)
                else:
                    if m_hit_id[rname]["left"].has_key(clip_type):
                        ori=m_hit_id[rname]["left"][clip_type][1]
                        if map_length>ori:
                            m_hit_id[rname]["left"][clip_type]=(flag,map_length, map_pos)
                    else:
                        m_hit_id[rname]["left"][clip_type]=(flag,map_length, map_pos)
            elif qname_field[-1]=="right":
                if m_hit_id.has_key(rname)==False:
                    m_hit_id[rname]={}
                    m_hit_id[rname]["right"]={}
                    m_hit_id[rname]["right"][clip_type]=(flag,map_length, map_pos)
                elif m_hit_id[rname].has_key("right")==False:
                    m_hit_id[rname]["right"]={}
                    m_hit_id[rname]["right"][clip_type]=(flag,map_length, map_pos)
                else:
                    if m_hit_id[rname]["right"].has_key(clip_type):
                        ori=m_hit_id[rname]["right"][clip_type][1]
                        if map_length>ori:
                            m_hit_id[rname]["right"][clip_type]=(flag,map_length, map_pos)
                    else:
                        m_hit_id[rname]["right"][clip_type]=(flag,map_length, map_pos)
            else:
                print "Error with flank id ", id

    #print m_hit_id
    for rname in m_hit_id:
        if len(m_hit_id[rname])==2:
            b_l_u=m_hit_id[rname]["left"].has_key(UNCLIP)
            b_l_l=m_hit_id[rname]["left"].has_key(LEFT_CLIP)
            b_l_r=m_hit_id[rname]["left"].has_key(RIGTH_CLIP)
            b_r_u=m_hit_id[rname]["right"].has_key(UNCLIP)
            b_r_l=m_hit_id[rname]["right"].has_key(LEFT_CLIP)
            b_r_r=m_hit_id[rname]["right"].has_key(RIGTH_CLIP)

            left_clip_type=-2
            right_clip_type=-2
            match_length=-1
            max_match_length=-1
            left_match_length=-1
            right_match_length=-1
            left_map_pos=-1
            right_map_pos=-1
            left_pick_type=-2
            right_pick_type=-2
            b_rc=False

            if b_l_u and b_r_u:
                left_clip_type=UNCLIP
                right_clip_type=UNCLIP
                left_flag=m_hit_id[rname]["left"][left_clip_type][0]
                right_flag=m_hit_id[rname]["right"][right_clip_type][0]
                match_length=m_hit_id[rname]["left"][left_clip_type][1]+m_hit_id[rname]["right"][right_clip_type][1]
                if (max_match_length < match_length) and ((left_flag & 16) == (right_flag & 16)):
                    max_match_length=match_length
                    left_pick_type=left_clip_type
                    right_pick_type=right_clip_type
                    left_match_length=m_hit_id[rname]["left"][left_clip_type][1]
                    right_match_length=m_hit_id[rname]["right"][right_clip_type][1]
                    left_map_pos=m_hit_id[rname]["left"][left_clip_type][2]
                    right_map_pos=m_hit_id[rname]["right"][right_clip_type][2]
                    if left_flag & 16 != 0:
                        b_rc=True
            if b_l_u and b_r_l:
                left_clip_type=UNCLIP
                right_clip_type=LEFT_CLIP
                left_flag=m_hit_id[rname]["left"][left_clip_type][0]
                right_flag=m_hit_id[rname]["right"][right_clip_type][0]
                match_length=m_hit_id[rname]["left"][left_clip_type][1]+m_hit_id[rname]["right"][right_clip_type][1]
                if (max_match_length < match_length) and ((left_flag & 16) == (right_flag & 16)):
                    max_match_length=match_length
                    left_pick_type=left_clip_type
                    right_pick_type=right_clip_type
                    left_match_length=m_hit_id[rname]["left"][left_clip_type][1]
                    right_match_length=m_hit_id[rname]["right"][right_clip_type][1]
                    left_map_pos=m_hit_id[rname]["left"][left_clip_type][2]
                    right_map_pos=m_hit_id[rname]["right"][right_clip_type][2]
                    if left_flag & 16 != 0:
                        b_rc=True
            if b_l_u and b_r_r:
                left_clip_type=UNCLIP
                right_clip_type=RIGTH_CLIP
                left_flag=m_hit_id[rname]["left"][left_clip_type][0]
                right_flag=m_hit_id[rname]["right"][right_clip_type][0]
                match_length=m_hit_id[rname]["left"][left_clip_type][1]+m_hit_id[rname]["right"][right_clip_type][1]
                if (max_match_length < match_length) and ((left_flag & 16) == (right_flag & 16)):
                    max_match_length=match_length
                    left_pick_type=left_clip_type
                    right_pick_type=right_clip_type
                    left_match_length=m_hit_id[rname]["left"][left_clip_type][1]
                    right_match_length=m_hit_id[rname]["right"][right_clip_type][1]
                    left_map_pos=m_hit_id[rname]["left"][left_clip_type][2]
                    right_map_pos=m_hit_id[rname]["right"][right_clip_type][2]
                    if left_flag & 16 != 0:
                        b_rc=True
            if b_l_l and b_r_u:
                left_clip_type=LEFT_CLIP
                right_clip_type=UNCLIP
                left_flag=m_hit_id[rname]["left"][left_clip_type][0]
                right_flag=m_hit_id[rname]["right"][right_clip_type][0]
                match_length=m_hit_id[rname]["left"][left_clip_type][1]+m_hit_id[rname]["right"][right_clip_type][1]
                if (max_match_length < match_length) and ((left_flag & 16) == (right_flag & 16)):
                    max_match_length=match_length
                    left_pick_type=left_clip_type
                    right_pick_type=right_clip_type
                    left_match_length=m_hit_id[rname]["left"][left_clip_type][1]
                    right_match_length=m_hit_id[rname]["right"][right_clip_type][1]
                    left_map_pos=m_hit_id[rname]["left"][left_clip_type][2]
                    right_map_pos=m_hit_id[rname]["right"][right_clip_type][2]
                    if left_flag & 16 != 0:
                        b_rc=True
            if b_l_l and b_r_r:
                left_clip_type=LEFT_CLIP
                right_clip_type=RIGTH_CLIP
                left_flag=m_hit_id[rname]["left"][left_clip_type][0]
                right_flag=m_hit_id[rname]["right"][right_clip_type][0]
                match_length=m_hit_id[rname]["left"][left_clip_type][1]+m_hit_id[rname]["right"][right_clip_type][1]
                if (max_match_length < match_length) and ((left_flag & 16) == (right_flag & 16)):
                    max_match_length=match_length
                    left_pick_type=left_clip_type
                    right_pick_type=right_clip_type
                    left_match_length=m_hit_id[rname]["left"][left_clip_type][1]
                    right_match_length=m_hit_id[rname]["right"][right_clip_type][1]
                    left_map_pos=m_hit_id[rname]["left"][left_clip_type][2]
                    right_map_pos=m_hit_id[rname]["right"][right_clip_type][2]
                    if left_flag & 16 != 0:
                        b_rc=True
            if b_l_r and b_r_u:
                left_clip_type=RIGTH_CLIP
                right_clip_type=UNCLIP
                left_flag=m_hit_id[rname]["left"][left_clip_type][0]
                right_flag=m_hit_id[rname]["right"][right_clip_type][0]
                match_length=m_hit_id[rname]["left"][left_clip_type][1]+m_hit_id[rname]["right"][right_clip_type][1]
                if (max_match_length < match_length) and ((left_flag & 16) == (right_flag & 16)):
                    max_match_length=match_length
                    left_pick_type=left_clip_type
                    right_pick_type=right_clip_type
                    left_match_length=m_hit_id[rname]["left"][left_clip_type][1]
                    right_match_length=m_hit_id[rname]["right"][right_clip_type][1]
                    left_map_pos=m_hit_id[rname]["left"][left_clip_type][2]
                    right_map_pos=m_hit_id[rname]["right"][right_clip_type][2]
                    if left_flag & 16 != 0:
                        b_rc=True
            if b_l_r and b_r_l:
                left_clip_type=RIGTH_CLIP
                right_clip_type=LEFT_CLIP
                left_flag=m_hit_id[rname]["left"][left_clip_type][0]
                right_flag=m_hit_id[rname]["right"][right_clip_type][0]
                match_length=m_hit_id[rname]["left"][left_clip_type][1]+m_hit_id[rname]["right"][right_clip_type][1]
                if (max_match_length < match_length) and ((left_flag & 16) == (right_flag & 16)):
                    max_match_length=match_length
                    left_pick_type=left_clip_type
                    right_pick_type=right_clip_type
                    left_match_length=m_hit_id[rname]["left"][left_clip_type][1]
                    right_match_length=m_hit_id[rname]["right"][right_clip_type][1]
                    left_map_pos=m_hit_id[rname]["left"][left_clip_type][2]
                    right_map_pos=m_hit_id[rname]["right"][right_clip_type][2]
                    if left_flag & 16 != 0:
                        b_rc=True

            if left_pick_type >= 0 and right_pick_type>=0:
                m_picked_id[rname]=[]
                m_picked_id[rname].append(left_map_pos)
                m_picked_id[rname].append(right_map_pos)
                m_picked_id[rname].append(left_match_length)
                m_picked_id[rname].append(right_match_length)
                m_picked_id[rname].append(b_rc)

    s_picked=""
    b_rc=True
    left_map_pos=-1
    right_map_pos=-1
    left_map_length=-1
    right_map_length=-1
    min_lenth=-1

    for key in m_picked_id: ##choose one out of the satisfied ones
        #if (m_picked_id[key][2]+m_picked_id[key][3])>min_lenth:
        #    s_picked=key
        if m_picked_id.has_key(key)==False:
            return
        left_map_pos=m_picked_id[key][0]
        right_map_pos=m_picked_id[key][1]
        left_map_length=m_picked_id[key][2]
        right_map_length=m_picked_id[key][3]
        b_rc=m_picked_id[key][4]

        start_pos=-1
        end_pos=-1
        if b_rc==True:
            end_pos=left_map_pos
            start_pos=right_map_pos+right_map_length
        else:
            start_pos=left_map_pos+left_map_length
            end_pos=right_map_pos
        if (end_pos-start_pos) > min_lenth:
            s_picked=key
            min_lenth=end_pos-start_pos

    if s_picked!="": #if find one
        if m_picked_id.has_key(s_picked)==False:
            return
        left_map_pos=m_picked_id[s_picked][0]
        right_map_pos=m_picked_id[s_picked][1]
        left_map_length=m_picked_id[s_picked][2]
        right_map_length=m_picked_id[s_picked][3]
        b_rc=m_picked_id[s_picked][4]

        sf_picked=working_folder+"velvet_temp/{0}/picked_seqs.fa".format(id)
        sf_picked_contigs=working_folder+"velvet_temp/{0}/picked_contigs.fa".format(id)
        f_picked_contigs=open(sf_picked_contigs,"w")
        with open(sf_picked,"w") as fout_picked:
            for record in SeqIO.parse(sf_contig, "fasta"):
                if str(record.id)== s_picked:
                    s_gap=""
                    s_ori=str(record.seq)
                    s_contig=s_ori
                    if b_rc==True:
                        end_pos=left_map_pos
                        start_pos=right_map_pos+right_map_length
                        s_gap=gnrt_reverse_complementary(s_ori[start_pos-1:end_pos]) ##generate reverse complementary of the gap seq
                        s_contig=gnrt_reverse_complementary(s_ori)
                    else:
                        start_pos=left_map_pos+left_map_length
                        end_pos=right_map_pos
                        s_gap=s_ori[start_pos-1:end_pos]

                    if s_gap!="":
                        fout_picked.write(">"+id+"_"+str(record.id)+"\n")
                        fout_picked.write(s_gap+"\n")

                    if s_contig!="":
                        f_picked_contigs.write(">"+id+"_"+str(record.id)+"\n")
                        f_picked_contigs.write(s_contig+"\n")
        f_picked_contigs.close()


def run_pick_extended_contig(id):
    #print id ##########################################################################################################################3
    global working_folder
    global LEFT_CLIP
    global RIGTH_CLIP
    global UNCLIP
    global BOTH_CLIP

    sf_contig=working_folder+"velvet_temp/{0}/contigs.fa".format(id)
    if os.path.exists(sf_contig)==False:
        return
    m_contigs={}
    for record in SeqIO.parse(sf_contig, "fasta"):
        m_contigs[str(record.id)]=str(record.seq)

    sf_flank_alnmt=working_folder+"velvet_temp/{0}/flanks.sam".format(id)
    m_hit_id={}
    with open(sf_flank_alnmt) as fin_algnmt:
        for line in fin_algnmt:
            fields=line.split()
            cigar=fields[5]
            flag=int(fields[1])
            if cigar=="*":
                continue
            b_rc=False
            if flag*16!=0:
                b_rc=True

            clip_type, map_length=get_clip_type_length(cigar)
            if clip_type==UNCLIP or clip_type==BOTH_CLIP:
                continue

            map_pos=int(fields[3])
            rname=fields[2]
            qname=fields[0]

            qname_field=qname.split("_")
            if qname_field[-1]=="left":
                if (b_rc==True and clip_type==LEFT_CLIP) or (b_rc==False and clip_type==RIGTH_CLIP):
                    continue
                if m_hit_id.has_key(qname)==False:
                    m_hit_id[qname]={}
                    m_hit_id[qname]["left"]={}
                    m_hit_id[qname]["left"][rname]=(map_pos, map_length, b_rc)
                else:
                    if m_hit_id[qname]["left"].has_key(rname)==True:
                        ori_len=m_hit_id[qname]["left"][rname][1]
                        if map_length>ori_len:
                            m_hit_id[qname]["left"][rname]=(map_pos, map_length, b_rc)
                    else:
                        m_hit_id[qname]["left"][rname]=(map_pos, map_length, b_rc)

            elif qname_field[-1]=="right":
                if (b_rc==True and clip_type==RIGTH_CLIP) or (b_rc==False and clip_type==LEFT_CLIP):
                    continue

                if m_hit_id.has_key(qname)==False:
                    m_hit_id[qname]={}
                    m_hit_id[qname]["right"]={}
                    m_hit_id[qname]["right"][rname]=(map_pos, map_length, b_rc)
                else:
                    if m_hit_id[qname]["right"].has_key(rname)==True:
                        ori_len=m_hit_id[qname]["right"][rname][1]
                        if map_length>ori_len:
                            m_hit_id[qname]["right"][rname]=(map_pos, map_length, b_rc)
                    else:
                        m_hit_id[qname]["right"][rname]=(map_pos, map_length, b_rc)
            else:
                print "Error with flank id ", id


    s_left_picked=""
    s_right_picked=""
    min_lenth=0
    for qname in m_hit_id:
        min_lenth=0
        if m_hit_id[qname].has_key("left"):
            for rname in m_hit_id[qname]["left"]:
                if m_hit_id[qname]["left"].has_key(rname)==False:
                    continue
                if m_hit_id[qname]["left"][rname][1]>min_lenth:
                    min_lenth=m_hit_id[qname]["left"][rname][1]
                    s_left_picked=rname
                elif m_hit_id[qname]["left"][rname][1]==min_lenth and s_left_picked!="" \
                        and m_contigs.has_key(s_left_picked) and m_contigs.has_key(rname):
                    if len(m_contigs[rname]) > m_contigs[s_left_picked]:
                        s_left_picked=rname

        min_lenth=0
        if m_hit_id[qname].has_key("right"):#
            for rname in m_hit_id[qname]["right"]:
                if m_hit_id[qname]["right"].has_key(rname)==False:
                    continue
                if m_hit_id[qname]["right"][rname][1]>min_lenth:
                    min_lenth=m_hit_id[qname]["right"][rname][1]
                    s_right_picked=rname
                elif m_hit_id[qname]["right"][rname][1]==min_lenth and s_right_picked!="" \
                        and m_contigs.has_key(s_right_picked) and m_contigs.has_key(rname):
                    if len(m_contigs[rname])> m_contigs[s_right_picked]:
                        s_right_picked=rname

    #print "L_"+s_left_picked, "R_"+s_right_picked ############################################################################
    #print "222222222222222222222222222222"##################################################################
    s_left_seq=""
    s_right_seq=""
    b_rc_left=True
    b_rc_right=True

    s_contig=""
    if s_left_picked!="" and s_left_picked==s_right_picked:
        l_map_length=m_hit_id[id+"_left"]["left"][s_left_picked][1]
        b_rc_left=m_hit_id[id+"_left"]["left"][s_left_picked][2]
        b_rc_right=m_hit_id[id+"_right"]["right"][s_right_picked][2]
        r_map_length=m_hit_id[id+"_right"]["right"][s_right_picked][1]

        if l_map_length> r_map_length:
            if m_contigs.has_key(s_left_picked)==False:
                return
            map_pos=m_hit_id[id+"_left"]["left"][s_left_picked][0]
            if b_rc_left==True:
                s_left_seq=m_contigs[s_left_picked][0:map_pos]
            else:
                s_left_seq=m_contigs[s_left_picked][map_pos+l_map_length-1:]
            s_contig=m_contigs[s_left_picked]
        else:
            if m_contigs.has_key(s_right_picked)==False:
                return
            map_pos=m_hit_id[id+"_right"]["right"][s_right_picked][0]
            if b_rc_right==False:
                s_right_seq=m_contigs[s_right_picked][0:map_pos]
            else:
                s_right_seq=m_contigs[s_right_picked][map_pos+r_map_length-1:]
            s_contig=m_contigs[s_right_picked]
    else:
        if s_left_picked!="":
            if m_contigs.has_key(s_left_picked)==False:
                return
            map_pos=m_hit_id[id+"_left"]["left"][s_left_picked][0]
            map_length=m_hit_id[id+"_left"]["left"][s_left_picked][1]
            b_rc_left=m_hit_id[id+"_left"]["left"][s_left_picked][2]
            if b_rc_left==True:
                s_left_seq=m_contigs[s_left_picked][0:map_pos]
            else:
                s_left_seq=m_contigs[s_left_picked][map_pos+map_length-1:]
            s_contig=m_contigs[s_left_picked]

        if s_right_picked!="":
            if m_contigs.has_key(s_right_picked)==False:
                return
            map_pos=m_hit_id[id+"_right"]["right"][s_right_picked][0]
            map_length=m_hit_id[id+"_right"]["right"][s_right_picked][1]
            b_rc_right=m_hit_id[id+"_right"]["right"][s_right_picked][2]
            if b_rc_right==False:
                s_right_seq=m_contigs[s_right_picked][0:map_pos-1]
            else:
                s_right_seq=m_contigs[s_right_picked][map_pos+map_length-1:]
            s_contig=s_contig+"NN"+m_contigs[s_right_picked]

    s_seq=""
    if b_rc_left==False and b_rc_right==False:
        s_seq=s_left_seq+"NN"+s_right_seq
    elif b_rc_left==True and b_rc_right==True:
        s_seq=gnrt_reverse_complementary(s_left_seq)+"NN"+gnrt_reverse_complementary(s_right_seq)
    elif b_rc_left==False and b_rc_right==True:
        s_seq=s_left_seq+"NN"+gnrt_reverse_complementary(s_right_seq)
    elif b_rc_left==True and b_rc_right==False:
        s_seq=gnrt_reverse_complementary(s_left_seq)+"NN"+s_right_seq

    if s_seq!="" and s_seq!="NN":
        sf_picked=working_folder+"velvet_temp/{0}/picked_seqs.fa".format(id)
        with open(sf_picked,"w") as fout_picked:
            fout_picked.write(">"+id+"_"+s_left_picked+"_"+s_right_picked+"_extended"+"\n")
            fout_picked.write(s_seq+"\n")

    if s_contig!="" and s_contig!="NN":
        sf_contig=working_folder+"velvet_temp/{0}/picked_contigs.fa".format(id)
        with open(sf_contig,"w") as fout_picked:
            fout_picked.write(">"+id+"_"+s_left_picked+"_"+s_right_picked+"_extended"+"\n")
            fout_picked.write(s_contig+"\n")


class ContigsSelection():
    def __init__(self, working_space):
        global working_folder
        global bwa_path
        global samtools_path
        working_folder=working_space
        bwa_path=get_bwa_path()
        samtools_path=get_samtools_path()
        self.nthreads=get_threads_num()


    def pick_full_constructed_contigs(self, bwa_score, fa_list, sf_picked):
        global working_folder
        global bwa_min_score
        bwa_min_score=bwa_score
        pool = Pool(self.nthreads)
        pool.map(run_pick_full_constructed_contig, fa_list, 1)
        pool.close()
        pool.join()

        #with open("picked_seqs_round0.fa","a") as fout_seqs:
        #print fa_list ##########################################################################################################
        for key in fa_list:
            sf_tmp=working_folder+"velvet_temp/{0}/picked_seqs.fa".format(key)
            if os.path.exists(sf_tmp)==True:
                cmd="cat {0} >> {1}".format(sf_tmp, sf_picked)
                Popen(cmd, shell = True, stdout = PIPE).communicate()
            sf_tmp=working_folder+"velvet_temp/{0}/picked_contigs.fa".format(key)
            if os.path.exists(sf_tmp)==True:
                cmd="cat {0} >> {1}".format(sf_tmp, sf_picked+"_ori.txt")
                Popen(cmd, shell = True, stdout = PIPE).communicate()

    def get_already_picked(self, sf_picked):
        l_picked={}
        with open(sf_picked) as fin_picked:
            for line in fin_picked:
                if line[0]==">":
                    fields=line[1:].split("_")
                    l_picked[fields[0]+"_"+fields[1]]=1
        return l_picked


    def pick_extended_contigs(self, bwa_score, fa_list, sf_picked):
        global working_folder
        global bwa_min_score
        bwa_min_score=bwa_score
        pool = Pool(self.nthreads)
        pool.map(run_pick_extended_contig, fa_list, 1)
        pool.close()
        pool.join()

        for key in fa_list:
            sf_tmp=working_folder+"velvet_temp/{0}/picked_seqs.fa".format(key)
            if os.path.exists(sf_tmp)==True:
                cmd="cat {0} >> {1}".format(sf_tmp, sf_picked)
                print cmd
                Popen(cmd, shell = True, stdout = PIPE).communicate()

            sf_tmp=working_folder+"velvet_temp/{0}/picked_contigs.fa".format(key)
            if os.path.exists(sf_tmp)==True:
                cmd="cat {0} >> {1}".format(sf_tmp, sf_picked+"_ori.txt")
                print cmd
                Popen(cmd, shell = True, stdout = PIPE).communicate()
