import sys
import os
from Bio import SeqIO
from subprocess import *
from multiprocessing.dummy import Pool as ThreadPool


def gnrt_discordant_read_list(sf_fai,sf_out):
    sf_left=sf_out+"_left.txt"
    sf_right=sf_out+"_right.txt"

    with open(sf_left,"w") as fout_discordant_left:
        with open(sf_right,"w") as fout_discordant_right:
            with open(sf_fai) as fin_fai:
                for line in fin_fai:
                    fields=line.split()
                    scaffold_id=fields[0]

                    sf_scaffold_right_list="./scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(scaffold_id)
                    sf_scaffold_left_list="./scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(scaffold_id)
                    with open(sf_scaffold_left_list) as fin_left_list:
                        for record in fin_left_list:
                            info_fields=record.split()
                            if info_fields[3]=="discordant":
                                read_id=info_fields[0]
                                fout_discordant_left.write(read_id+"\n")

                    with open(sf_scaffold_right_list) as fin_right_list:
                        for record in fin_right_list:
                            info_fields=record.split()
                            if info_fields[3]=="discordant":
                                read_id=info_fields[0]
                                fout_discordant_right.write(read_id+"\n")


def get_discordant_of_scaffold():
    scaffold_id="gi|317021970|gb|GL583003.1|"
    left_reads_list={}
    right_reads_list={}

    reads_temp={}

    sf_id_left_list="./scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(scaffold_id)
    sf_id_right_list="./scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(scaffold_id)

    with open(sf_id_left_list) as fin_left:
        for line in fin_left:
            fields=line.split()
            read_id=fields[0]
            gap_id=fields[1]

            left_reads_list[read_id]=gap_id

    with open(sf_id_right_list) as fin_right:
        for line in fin_right:
            fields=line.split()
            read_id=fields[0]
            gap_id=fields[1]

            right_reads_list[read_id]=gap_id

    print "Parse list done!!"

    sf_reads_left="discordant_all_left.fastq"#"discordant_left_anchor_mapq30.fastq"
    sf_reads_right="discordant_all_right.fastq"#"discordant_right_anchor_mapq30.fastq"

    #parse left reads
    for record in SeqIO.parse(sf_reads_left, "fastq"):
        read_id=str(record.id)
        if left_reads_list.has_key(read_id)==False:
            continue
        gap_id=left_reads_list[read_id]

        if reads_temp.has_key(gap_id):
            reads_temp[gap_id].append(">"+read_id+"_1\n"+str(record.seq)+"\n")
            temp_size=len(reads_temp[gap_id])
            if temp_size > 5000:
                sf_temp="discordant_alignment_all/{0}.fa".format(gap_id)
                f_temp=open(sf_temp,"a+")
                for i in range(temp_size):
                    f_temp.write(reads_temp[gap_id][i])
                f_temp.close()
                del reads_temp[gap_id][:]
        else:
            reads_temp[gap_id]=[]
            reads_temp[gap_id].append(">"+read_id+"_1\n"+str(record.seq)+"\n")

    print "Parse left reads done!!"

    #parse right reads
    for record in SeqIO.parse(sf_reads_right, "fastq"):
        read_id=str(record.id)
        if right_reads_list.has_key(read_id)==False:
            continue
        gap_id=right_reads_list[read_id]

        if reads_temp.has_key(gap_id):
            reads_temp[gap_id].append(">"+read_id+"_1\n"+str(record.seq)+"\n")
            temp_size=len(reads_temp[gap_id])
            if temp_size > 5000:
                sf_temp="discordant_alignment_all/{0}.fa".format(gap_id)
                f_temp=open(sf_temp,"a+")
                for i in range(temp_size):
                    f_temp.write(reads_temp[gap_id][i])
                f_temp.close()
                del reads_temp[gap_id][:]
        else:
            reads_temp[gap_id]=[]
            reads_temp[gap_id].append(">"+read_id+"_1\n"+str(record.seq)+"\n")

    print "Parse right reads done!!"

    #clear the ones in memory
    total=0
    total_aligned=0
    for key in reads_temp:
        temp_size=len(reads_temp[key])
        sf_temp="discordant_alignment_all/{0}.fa".format(key)
        f_temp=open(sf_temp,"a+")
        for i in range(temp_size):
            f_temp.write(reads_temp[key][i])
        f_temp.close()
        del reads_temp[key][:]

        ref="/data2/chongchu/GapFilling/Assemble_gap_gi_317021970/gap_ref/{0}.fasta".format(key)
        cmd="bwa mem {0} discordant_alignment_all/{1}.fa | " \
            "samtools view -S -F 4 - > ./discordant_alignment_all/{2}.sam".format(ref,key,key)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="wc -l ./discordant_alignment_all/{0}.sam".format(key)
        rtn=Popen(cmd, shell = True, stdout = PIPE).communicate()
        scnt=str(rtn[0]).split()
        cnt=int(scnt[0])
        total_aligned = total_aligned+cnt

        sf_fa="./discordant_alignment_all/{0}.fa".format(key)
        cmd="wc -l {0} ".format(sf_fa)
        rtn=Popen(cmd, shell = True, stdout = PIPE).communicate()
        scnt=str(rtn[0]).split()
        cnt=int(scnt[0])/2
        total = total+cnt

    print total_aligned,total


def cnt_mapped():
    total_aligned=0
    total=0
    for i in range(1500):
        sf_sam="./discordant_alignment/{0}.sam".format(i)
        if os.path.exists(sf_sam)==False:
            continue

        cmd="samtools view ./discordant_alignment/{0}.sam | wc -l - ".format(i)
        rtn=Popen(cmd, shell = True, stdout = PIPE).communicate()
        scnt=str(rtn[0]).split()
        cnt=int(scnt[0])
        total_aligned=total_aligned+cnt

        sf_fa="./discordant_alignment/{0}.fa".format(i)
        if os.path.exists(sf_fa)==False:
            continue
        cmd="wc -l {0} ".format(sf_fa)
        rtn=Popen(cmd, shell = True, stdout = PIPE).communicate()
        scnt=str(rtn[0]).split()
        cnt=int(scnt[0])/2
        total=total+cnt

    print total_aligned, total


#/data2/chongchu/ParseFastqById/seqtk
#/data2/chongchu/GapFilling/Assemble_gap_all/scaffold_reads_list/"gi|317021970|gb|GL583003.1|_cluster_by_gap_reads_right.list"
#/data2/chongchu/GapFilling/Assemble_gap_all/scaffold_reads_list/"gi|317021970|gb|GL583003.1|_cluster_by_gap_reads_left.list"
#gnrt_discordant_read_list(sys.argv[1],sys.argv[2])

get_discordant_of_scaffold()
#cnt_mapped()
