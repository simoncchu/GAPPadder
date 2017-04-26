import sys
import os
from subprocess import *
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from Bio import SeqIO
from Utility import *

samtools_path=""
working_folder=""
m_both_unmap_ref={}
m_both_unmap_reads={}

def run_collect_both_unmapped(sf_bam): ###########################################################################################################
    global samtools_path
    sf_both_unmap=sf_bam+".both_unmapped.sam"
    cmd="{0} view -f 12 {1} > {2}".format(samtools_path, sf_bam, sf_both_unmap)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    sf_fq=sf_bam+".both_unmapped.fq"
    f_unmap=open(sf_fq,"w")
    with open(sf_both_unmap) as fin_sam:
        for line in fin_sam:
            fields=line.split()
            flag=int(fields[1])

            if flag>128:
                f_unmap.write("@"+fields[0]+"_2\n")
            else:
                f_unmap.write("@"+fields[0]+"_1\n")
            f_unmap.write(fields[9]+"\n")
            f_unmap.write("+\n")
            f_unmap.write(fields[10]+"\n")
    f_unmap.close()


def run_alignment_unmapped_reads(id):
    #if self.m_both_unmap_ref.has_key(id)==False:
    #    return
    global samtools_path
    global working_folder
    global m_both_unmap_reads
    global m_both_unmap_ref
    m_temp={}
    sf_algn=working_folder+"both_unmapped_2_gap_seqs.sort.bam"
    cmd="{0} view {1} ".format(samtools_path, sf_algn)

    for sctg in m_both_unmap_ref[id]:
        cmd=cmd+" {0}".format(id+"-"+sctg)

    succd=True
    try:
        buffer = Popen(cmd,shell = True,stdout=PIPE,stderr=PIPE)
        output,error = buffer.communicate()
        if output:
            lines=output.split("\n")
            for line in lines:
                #print line ##########################################################
                fields=line.split()
                if len(fields)<11:
                    continue
                mapq=int(fields[4])
                flag=int(fields[1])
                name_field=fields[0]
                #if mapq>=30 and mapq!=255:
                if mapq!=255:
                    #name_fields=fields[0].split("_")
                    bfirst=True
                    s_read_id=name_field+"_1"
                    #if name_fields[-1]=="2":
                    if flag&128!=0:
                        bfirst=False
                        s_read_id=name_field+"_2"

                    if m_both_unmap_reads.has_key(s_read_id)==False:
                        print s_read_id, "Error happen!!!"
                        return
                    m_temp[s_read_id]=m_both_unmap_reads[s_read_id]

                    read_ref=fields[2]
                    mate_ref=fields[6]

                    g_read=""
                    g_mate=""
                    if mate_ref!="=":
                        group=read_ref.split("-")
                        g_read=group[0]
                        group_mate=mate_ref.split("-")
                        g_mate=group_mate[0]

                    s_mate_head=""
                    if flag&8!=0 or g_read!=g_mate: #mate unmapped or not in same group
                        if bfirst==True:
                            s_mate_head=name_field+"_2"
                        else:
                            s_mate_head=name_field+"_1"


                        if m_both_unmap_reads.has_key(s_mate_head)==False:
                            #error
                            print s_mate_head, "Error happen!!!!!!!"
                            return

                        m_temp[s_mate_head]=m_both_unmap_reads[s_mate_head]
        if error:
            succd=False
        del buffer
    except:
        succd=False

    if succd==False:
        return

    sf_fq=working_folder+"gap_reads/{0}.fastq".format(id)
    with open(sf_fq,"a") as fout_fq:
        for shead in m_temp:
            fout_fq.write("@"+shead+"\n")
            fout_fq.write(m_both_unmap_reads[shead])

    sf_fq=working_folder+"unmapped_reads/{0}.fastq".format(id)
    with open(sf_fq,"w") as fout_fq:
        for shead in m_temp:
            fout_fq.write("@"+shead+"\n")
            fout_fq.write(m_both_unmap_reads[shead])

class BothUnmappedReadsCollector():
    def __init__(self, working_space):
        global working_folder
        working_folder=working_space
        global samtools_path
        samtools_path=get_samtools_path()
        self.bwa_path=get_bwa_path()
        self.nthreads=get_threads_num()


    def align_unmapped_to_contigs(self, fa_list):
        global working_folder
        global samtools_path
        global m_both_unmap_ref
        global m_both_unmap_reads

        sf_ref=working_folder+"gap_contigs_all.fa"
        sf_left_reads=working_folder+"both_unmapped_1.fq"
        sf_right_reads=working_folder+"both_unmapped_2.fq"
        sf_algn=working_folder+"both_unmapped_2_gap_seqs.sort"

        with open(sf_ref,"w") as fout_ref:
            for key in fa_list:
                sf_ctg=working_folder+"velvet_temp/{0}/contigs.fa".format(key)
                if os.path.exists(sf_ctg)==False:
                    continue
                for record in SeqIO.parse(sf_ctg, "fasta"):
                    fout_ref.write(">"+key+"-"+str(record.id)+"\n")
                    fout_ref.write(str(record.seq)+"\n")

        cmd="{0} index {1}".format(self.bwa_path, sf_ref)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="{0} mem -a -t 10 {1} {2} {3} | {4} view -h -S -b - | {5} sort - -o {6}.bam"\
            .format(self.bwa_path, sf_ref, sf_left_reads, sf_right_reads, samtools_path, samtools_path, sf_algn)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="{0} index {1}".format(samtools_path, sf_algn+".bam")
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="{0} view -H {1}".format(samtools_path, sf_algn+".bam")
        buffer=Popen(cmd, shell = True, stdout = PIPE)
        fa_list_hit=[]
        for line in buffer.stdout:
            fields=line.split()
            if fields[0]!="@SQ":
                continue
            ref_fields=fields[1].split(":")
            s_contig_fields=ref_fields[1].split("-")
            if self.m_both_unmap_ref.has_key(s_contig_fields[0])==False:
                self.m_both_unmap_ref[s_contig_fields[0]]={}
                fa_list_hit.append(s_contig_fields[0])
            self.m_both_unmap_ref[s_contig_fields[0]][s_contig_fields[1]]=1

        pool1 = Pool(self.nthreads)
        pool1.map(run_alignment_unmapped_reads, fa_list_hit, 1)
        pool1.close()
        pool1.join()

        del m_both_unmap_reads
        del m_both_unmap_ref


    def collect_both_unmapped_reads(self, bam_list, id_list):
        global m_both_unmap_reads
        global m_both_unmap_ref
        global  working_folder

        pool = ThreadPool(3)
        pool.map(run_collect_both_unmapped, bam_list)
        pool.close()
        pool.join()

        sf_both_unmap=working_folder+"both_unmapped.fq"
        cmd="cat "
        for sf_bam in bam_list:
            cmd=cmd + sf_bam +".both_unmapped.fq "
        cmd=cmd+" > {0}".format(sf_both_unmap)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        #read in all the unmapped reads, and keep in map
        cnt=1
        sread=""
        shead=""
        with open(sf_both_unmap) as fin_reads:
            for line in fin_reads:
                if cnt==1:
                    stemp=line.rstrip()
                    shead=stemp[1:]
                else:
                    sread=sread+line.rstrip()+"\n"

                if cnt%4==0:
                    m_both_unmap_reads[shead]=sread
                    sread=""
                    cnt=0
                cnt=cnt+1

        sf_left=working_folder+"both_unmapped_1.fq"
        sf_right=working_folder+"both_unmapped_2.fq"
        f_left=open(sf_left,"w")
        f_right=open(sf_right,"w")
        for key in m_both_unmap_reads:
            read_id=key[:-2]
            if key[-1]=="1":
                f_left.write("@"+read_id+"\n")
                f_left.write(m_both_unmap_reads[key])

                f_right.write("@"+read_id+"\n")
                f_right.write(m_both_unmap_reads[read_id+"_2"])
        f_left.close()
        f_right.close()

        self.align_unmapped_to_contigs(id_list, self.nthreads)