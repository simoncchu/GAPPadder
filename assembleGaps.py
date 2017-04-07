import sys
import os
from subprocess import *
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from MergeContigs import *
from Bio import SeqIO

m_kmer_len=[]
m_both_unmap_reads={}
m_both_unmap_ref={}
m_flanks={}

BOTH_CLIP=1
LEFT_CLIP=2
RIGTH_CLIP=3
UNCLIP=4

BWA_MIN_SCORE=30
bunch_size=1

def run_cmd(cmd1):
    #print cmd1
    check_output(cmd1, shell=True)

def merge_reads(sf_fai, sf_gap_pos, n):
    m_scaffold_id={}
    cnt=0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1

    cmd_list=[]
    cnt=1
    pre_id=0
    with open(sf_gap_pos) as fin_gap_pos:#for each gap
        for line in fin_gap_pos:
            fields=line.split()
            scaffold=fields[3]
            scaffold_id=m_scaffold_id[scaffold]

            if pre_id!=scaffold_id:
                cnt=1

            sf_fq="{0}_{1}.fastq".format(scaffold_id,cnt)

            sf_fq1="/data2/chongchu/GapFilling/Assemble_gap_pe/gap_reads/{0}".format(sf_fq)
            sf_fq2="/data2/chongchu/GapFilling/Assemble_gap_mp_unkown/gap_reads/{0}".format(sf_fq)

            cnt=cnt+1
            pre_id=scaffold_id

            cmd=""
            if os.path.exists(sf_fq1)==True and os.path.exists(sf_fq2)==True:
                cmd="cat {0} {1} > gap_reads/{2}".format(sf_fq1, sf_fq2, sf_fq)
            elif os.path.exists(sf_fq1)==True and os.path.exists(sf_fq2)==False:
                cmd="cp {0} gap_reads/{1}".format(sf_fq1,sf_fq)
                #print "1"
            elif os.path.exists(sf_fq1)==False and os.path.exists(sf_fq2)==True:
                cmd="cp {0} gap_reads/{1}".format(sf_fq2,sf_fq)
                #print "2"

            if cmd!="":
                cmd_list.append(cmd)

    pool = ThreadPool(n)
    pool.map(run_cmd, cmd_list, bunch_size)
    pool.close()
    pool.join()


def merge_reads_v2(sf_fai, sf_gap_pos, merge_folder_list, reads_folder, n):
    m_scaffold_id={}
    cnt=0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1

    cmd_list=[]
    cnt=1
    pre_id=0
    with open(sf_gap_pos) as fin_gap_pos:#for each gap
        for line in fin_gap_pos:
            fields=line.split()
            scaffold=fields[3]
            scaffold_id=m_scaffold_id[scaffold]

            if pre_id!=scaffold_id:
                cnt=1

            sf_fq="{0}_{1}.fastq".format(scaffold_id,cnt)
            cnt=cnt+1
            pre_id=scaffold_id

            cmd="cat "
            for folder in merge_folder_list:
                sf_fq_x=folder+"/{0}/{1}".format(reads_folder, sf_fq)
                if os.path.exists(sf_fq_x)==True:
                    cmd=cmd+" {0}".format(sf_fq_x)

            if cmd!="cat ":
                cmd=cmd+" > {0}/{1}".format(reads_folder, sf_fq)
                cmd_list.append(cmd)

    pool = ThreadPool(n)
    pool.map(run_cmd, cmd_list, bunch_size)
    pool.close()
    pool.join()

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

def run_assembly(id):
    global m_kmer_len

    # fields=sf_fq.split("/")
    # fld_len=len(fields)
    # id_fields=fields[fld_len-1].split(".")
    # id=id_fields[0]

    for klen in m_kmer_len:
        klen_fields=klen.split("_")
        kmer_len=int(klen_fields[0])
        asm_len=int(klen_fields[1])

        pth="./kmc_temp/{0}".format(id)
        if os.path.exists(pth)==False:
            cmd="mkdir ./kmc_temp/{0}".format(id)
            #print cmd
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="/data2/chongchu/tools/kmc2.3/kmc -k{0} -cs10000000 -m52 gap_reads/{1}.fastq  " \
            "temp/{2}.res ./kmc_temp/{3}".format(kmer_len,id,id,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="/data2/chongchu/tools/kmc2.3/kmc_dump -ci0 temp/{0}.res temp/{1}.dump".format(id,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="awk -f ./cvtKMC_2_Fq.awk temp/{0}.dump > kmers/{1}_kmers.fq".format(id,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        pth="./velvet_temp/{0}".format(id)
        if os.path.exists(pth)==False:
            cmd="mkdir ./velvet_temp/{0}".format(id)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="/data2/chongchu/tools/velvet-master/velveth ./velvet_temp/{0} {1} -fastq " \
            "-short kmers/{2}_kmers.fq".format(id, asm_len, id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="/data2/chongchu/tools/velvet-master/velvetg ./velvet_temp/{0} -min_contig_lgth 50".format(id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="mv ./velvet_temp/{0}/contigs.fa ./velvet_temp/{1}/contigs_{2}_{3}.fa".format(id,id,kmer_len, asm_len)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    sf_merged="./velvet_temp/{0}/contigs.fa".format(id)
    with open(sf_merged,"w") as fout_merged:
        for klen in m_kmer_len:
            klen_fields=klen.split("_")
            ctg="./velvet_temp/{0}/contigs_{1}_{2}.fa".format(id, klen_fields[0], klen_fields[1])
            with open(ctg) as fin_ctg:
                for line in fin_ctg:
                    if line[0]==">":
                        fout_merged.write(">"+klen+"_"+line[1:])
                    else:
                        fout_merged.write(line)

    ##need to clear the temporary files here:
    run_clear(id)


def assembly(id_list, n):
    global m_kmer_len

    if os.path.exists("./empty_dir")==False:
        cmd="mkdir ./empty_dir"
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    ##clear old data
    cmd="rsync -a --delete empty_dir/ kmc_temp/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="rsync -a --delete empty_dir/ temp/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="rsync -a --delete empty_dir/ kmers/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="rsync -a --delete empty_dir/ velvet_temp/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    pool = Pool(n)
    pool.map(run_assembly, id_list, bunch_size)
    pool.close()
    pool.join()


def assembly_given_list(id_list, n):
    pool = Pool(n)
    pool.map(run_assembly, id_list, bunch_size)
    pool.close()
    pool.join()

def run_alignment_contig(id):
    sf_ref="./../gap_seqs/{0}.fa".format(id)
    if os.path.exists(sf_ref)==True:
        sf_bwt="./../gap_seqs/{0}.fa.bwt".format(id)
        if os.path.exists(sf_bwt)==False:
            cmd="bwa index {0}".format(sf_ref)
            Popen(cmd, shell = True, stdout = PIPE).communicate()
    else:
        return

    sf="./velvet_temp/{0}/contigs.fa".format(id)
    if os.path.exists(sf)==True:
        cmd="bwa mem {0} {1} | " \
            "samtools view -S -h -b - | samtools sort - -o merged_contigs_alignment/{2}.sort.bam " \
            "&& samtools index merged_contigs_alignment/{3}.sort.bam".format(sf_ref,sf,id,id)
        #print cmd#################
        Popen(cmd, shell = True, stdout = PIPE).communicate()

def run_alignment_reads(id):
    sf_ref="./../gap_seqs/{0}.fa".format(id)
    if os.path.exists(sf_ref)==True:
        sf_bwt="./../gap_seqs/{0}.fa.bwt".format(id)
        if os.path.exists(sf_bwt)==False:
            cmd="bwa index {0}".format(sf_ref)
            Popen(cmd, shell = True, stdout = PIPE).communicate()
    else:
        return

    sf="./gap_reads/{0}.fastq".format(id)
    if os.path.exists(sf)==True:
        cmd="bwa mem {0} {1} | " \
            "samtools view -S -h -b - | samtools sort - -o gap_reads_alignment/{2}.sort.bam " \
            "&& samtools index gap_reads_alignment/{3}.sort.bam".format(sf_ref,sf,id,id)
        #print cmd#################
        Popen(cmd, shell = True, stdout = PIPE).communicate()

def run_merge(sf_contig):
    #merge contigs
    sf_folder="./velvet_temp/{0}/".format(sf_contig)
    sf_contig=sf_folder+"contigs.fa"
    if os.path.exists(sf_contig)==False:
        return
    #print sf_folder
    #sf_folder="./{0}/".format(sf_contig)
    local_CONTIGS_MERGER_PATH="./ContigsMerger"
    merge_contigs(local_CONTIGS_MERGER_PATH, sf_folder, 5, 0.99, 0.99)


def prepare_list(sf_fai, sf_gap_pos):
    m_scaffold_id={}
    cnt=0
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            scaffold_id=fields[0]
            m_scaffold_id[scaffold_id]=cnt
            cnt=cnt+1

    fa_list=[]
    cnt=1
    pre_id=0

    with open(sf_gap_pos) as fin_gap_pos:#for each gap
        for line in fin_gap_pos:
            fields=line.split()
            scaffold=fields[3]
            scaffold_id=m_scaffold_id[scaffold]

            if pre_id!=scaffold_id:
                cnt=1

            sf_fa="{0}_{1}".format(scaffold_id,cnt)
            cnt=cnt+1

            pre_id=scaffold_id
            sf_fq="gap_reads/{0}.fastq".format(sf_fa)
            if os.path.exists(sf_fq):
                fa_list.append(sf_fa)
    return fa_list

def run_contigs_merge(n, fa_list):
    #print "Preparing list ..."

    print "Start merging..."
    pool=Pool(n)
    pool.map(run_merge, fa_list, bunch_size)
    pool.close()
    pool.join()


def run_alignment_unmapped_reads(id):
    global m_both_unmap_reads
    global m_both_unmap_ref

    #if m_both_unmap_ref.has_key(id)==False:
    #    return
    m_temp={}

    sf_algn="both_unmapped_2_gap_seqs.sort.bam"
    cmd="samtools view {0} ".format(sf_algn)

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

    sf_fq="gap_reads/{0}.fastq".format(id)
    with open(sf_fq,"a") as fout_fq:
        for shead in m_temp:
            fout_fq.write("@"+shead+"\n")
            fout_fq.write(m_both_unmap_reads[shead])

    sf_fq="unmapped_reads/{0}.fastq".format(id)
    with open(sf_fq,"w") as fout_fq:
        for shead in m_temp:
            fout_fq.write("@"+shead+"\n")
            fout_fq.write(m_both_unmap_reads[shead])


def run_collect_both_unmapped(sf_bam): ###########################################################################################################
    sf_both_unmap=sf_bam+".both_unmapped.sam"
    cmd="samtools view -f 12 {0} > {1}".format(sf_bam, sf_both_unmap)
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


def align_unmapped_to_contigs(fa_list, n):
    global m_both_unmap_reads
    global m_both_unmap_ref

    sf_ref="gap_contigs_all.fa"
    sf_left_reads="both_unmapped_1.fq"
    sf_right_reads="both_unmapped_2.fq"
    sf_algn="both_unmapped_2_gap_seqs.sort"

    with open(sf_ref,"w") as fout_ref:
        for key in fa_list:
            sf_ctg="./velvet_temp/{0}/contigs.fa".format(key)
            if os.path.exists(sf_ctg)==False:
                continue
            for record in SeqIO.parse(sf_ctg, "fasta"):
                fout_ref.write(">"+key+"-"+str(record.id)+"\n")
                fout_ref.write(str(record.seq)+"\n")

    cmd="bwa index {0}".format(sf_ref)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="bwa mem -a -t 10 {0} {1} {2} | samtools view -h -S -b - | " \
        "samtools sort - -o {3}.bam".format(sf_ref, sf_left_reads, sf_right_reads, sf_algn)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools index {0}".format(sf_algn+".bam")
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="samtools view -H {0}".format(sf_algn+".bam")
    buffer=Popen(cmd, shell = True, stdout = PIPE)
    fa_list_hit=[]
    for line in buffer.stdout:
        fields=line.split()
        if fields[0]!="@SQ":
            continue
        ref_fields=fields[1].split(":")
        s_contig_fields=ref_fields[1].split("-")
        if m_both_unmap_ref.has_key(s_contig_fields[0])==False:
            m_both_unmap_ref[s_contig_fields[0]]={}
            fa_list_hit.append(s_contig_fields[0])
        m_both_unmap_ref[s_contig_fields[0]][s_contig_fields[1]]=1

    pool1 = Pool(n)
    pool1.map(run_alignment_unmapped_reads, fa_list_hit, bunch_size)
    pool1.close()
    pool1.join()

    del m_both_unmap_reads
    del m_both_unmap_ref


def collect_both_unmapped_reads(bam_list, id_list, n):
    global m_both_unmap_reads

    pool = ThreadPool(3)
    pool.map(run_collect_both_unmapped, bam_list)
    pool.close()
    pool.join()

    sf_both_unmap="both_unmapped.fq"
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

    sf_left="both_unmapped_1.fq"
    sf_right="both_unmapped_2.fq"
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

    align_unmapped_to_contigs(id_list, n)


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

def run_collect_high_quality_unmap_to_contig_reads(id):
    sf_contig="velvet_temp/"+id+"/contigs.fa"
    #print sf_contig##########################################################################################
    # if os.path.exists(sf_contig+".bwt")==False:
    cmd="bwa index {0}".format(sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    #print cmd##############################################################################################33
    sf_reads="gap_reads_high_quality/{0}.fastq".format(id)
    if os.path.exists(sf_reads)==False:
        return
    sf_algmt="gap_reads_high_quality/{0}.sam".format(id)
    cmd="bwa mem -t 5 {0} {1} | samtools view -S - > {2}".format(sf_contig,sf_reads,sf_algmt)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    #print cmd##############################################################################################33
    #sf_picked_reads="gap_reads_high_quality/{0}.high_quality_unmap_to_contig.fa".format(id)
    #with open(sf_picked_reads,"w") as fout_picked:
    m_read_seq={}
    for record in SeqIO.parse(sf_reads, "fastq"):
        read_id=str(record.id)
        if m_read_seq.has_key(read_id)==False:
            m_read_seq[read_id]=str(record.seq)

    m_read_ref={}
    with open(sf_algmt) as fin_algnmt:
        for line in fin_algnmt:
            fields=line.split()
            cigar=fields[5]
            read_id=fields[0]
            if cigar=="*":
                continue

            rname=fields[2]
            if is_qualified_clipped(cigar, 1)==True: ##################################need a parameter to replace 20
                if m_read_ref.has_key(read_id)==False:
                    m_read_ref[read_id]={}
                m_read_ref[read_id][rname]=1

    sf_ori_contig="velvet_temp/"+id+"/original_contigs_before_merging.fa"
    cmd="rm {0}".format(sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mv {0} {1}".format(sf_ori_contig, sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    with open(sf_contig,"a") as fout_picked:
        for read_id in m_read_ref:
            if len(m_read_ref[read_id])>=2: ##at least clipped at two contigs
                fout_picked.write(">"+read_id+"\n")
                fout_picked.write(m_read_seq[read_id]+"\n")


def collect_high_quality_unmap_to_contigs_reads(n, id_list):
    #print id_list###############################################################################################3
    #id_list=[]
    #id_list=prepare_list(sf_fai, sf_gap_pos)
    #id_list.append("0_1139")
    pool = Pool(n)
    pool.map(run_collect_high_quality_unmap_to_contig_reads, id_list, bunch_size)
    pool.close()
    pool.join()


def run_filter(id):
    fcontig="velvet_temp/"+id+"/contigs.fa"
    cmd='bwa index {0}'.format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    #align reads to contigs
    sf_reads="gap_reads_alignment/"+id+".fastq"
    sf_sorted_bam="velvet_temp/"+id+"/alignment_for_coverage"
    cmd="bwa mem -a -t 5 {0} {1} | samtools view -h -S -b - | " \
        "samtools sort - -o {2}.sort.bam".format(fcontig, sf_reads,sf_sorted_bam)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="samtools faidx {0}".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #calculate coverage cover
    cmd="TERefiner_1 -G -r {0} -b {1} -c 0"
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    #"_cov_info_with_cutoff.txt"


def filter_contigs(sf_fai, sf_gap_pos, n):
    #align the reads back to contigs
    #calculate the average coverage
    #output the filtered contigs
    fa_list=prepare_list(sf_fai, sf_gap_pos)

    pool = Pool(n)
    pool.map(run_merge, fa_list,bunch_size)
    pool.close()
    pool.join()

    pool1 = Pool(n)
    pool1.map(run_filter, fa_list,bunch_size)
    pool1.close()
    pool1.join()


def get_clip_type_length(cigar):
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

def run_pick_full_constructed_contig(id):
    global m_flanks
    global BWA_MIN_SCORE
    sf_flank="./../flank_regions/{0}.fa".format(id)
    #if m_flanks.has_key(id)==False or len(m_flanks[id])!=2:
    if os.path.exists(sf_flank)==False:
        print "Wrong flank regions: ", id
        return

    sf_contig="velvet_temp/{0}/contigs.fa".format(id)
    cmd="bwa index {0}".format(sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #sf_flank="velvet_temp/{0}/flanks.fa".format(id)

    sf_flank_alnmt="velvet_temp/{0}/flanks.sam".format(id)
    cmd="bwa mem -T {0} -a {1} {2} | samtools view -S - > {3}".format(BWA_MIN_SCORE, sf_contig, sf_flank, sf_flank_alnmt)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    sf_pk="velvet_temp/{0}/picked_seqs.fa".format(id)
    if os.path.exists(sf_pk):
        cmd="rm {0}".format(sf_pk)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    sf_pk="velvet_temp/{0}/picked_contigs.fa".format(id)
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

        sf_picked="velvet_temp/{0}/picked_seqs.fa".format(id)
        sf_picked_contigs="velvet_temp/{0}/picked_contigs.fa".format(id)
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


def pick_full_constructed_contigs(bwa_min_score, n, fa_list, sf_picked):
    global m_flanks
    global BWA_MIN_SCORE

    BWA_MIN_SCORE=bwa_min_score
    # with open(sf_flank) as fin_flank:
    #     shead=""
    #     for line in fin_flank:
    #         if line[0]!=">":
    #             fields=shead[1:].split("_")
    #             id=fields[0]+"_"+fields[1]
    #             if m_flanks.has_key(id)==False:
    #                 m_flanks[id]=[]
    #             m_flanks[id].append(line.rstrip())
    #             if len(m_flanks[id])==2:
    #                 if os.path.exists("velvet_temp/{0}".format(id)):
    #                     sf_one_group="velvet_temp/{0}/flanks.fa".format(id)
    #                     with open(sf_one_group,"w") as fout_flank:
    #                         fout_flank.write(">"+id+"_left\n")
    #                         fout_flank.write(m_flanks[id][0]+"\n")
    #                         fout_flank.write(">"+id+"_right\n")
    #                         fout_flank.write(m_flanks[id][1]+"\n")
    #         else:
    #             shead=line.rstrip()

    #
    # #run the overlap program
    # #first pick out those both end overlapped with the two flank region ones
    #fa_list=prepare_list(sf_fai, sf_gap_pos)
    #fa_list=[]
    #fa_list.append("0_1221")
    #
    pool = Pool(n)
    pool.map(run_pick_full_constructed_contig, fa_list, bunch_size)
    pool.close()
    pool.join()

    #with open("picked_seqs_round0.fa","a") as fout_seqs:
    #print fa_list ##########################################################################################################
    for key in fa_list:
        sf_tmp="velvet_temp/{0}/picked_seqs.fa".format(key)
        if os.path.exists(sf_tmp)==True:
            cmd="cat {0} >> {1}".format(sf_tmp, sf_picked)
            Popen(cmd, shell = True, stdout = PIPE).communicate()
        sf_tmp="velvet_temp/{0}/picked_contigs.fa".format(key)
        if os.path.exists(sf_tmp)==True:
            cmd="cat {0} >> {1}".format(sf_tmp, sf_picked+"_ori.txt")
            Popen(cmd, shell = True, stdout = PIPE).communicate()

def get_already_picked(sf_picked):
    l_picked={}
    with open(sf_picked) as fin_picked:
        for line in fin_picked:
            if line[0]==">":
                fields=line[1:].split("_")
                l_picked[fields[0]+"_"+fields[1]]=1
    return l_picked


def run_pick_extended_contig(id):
    #print id ##########################################################################################################################3
    sf_contig="velvet_temp/{0}/contigs.fa".format(id)
    if os.path.exists(sf_contig)==False:
        return
    m_contigs={}
    for record in SeqIO.parse(sf_contig, "fasta"):
        m_contigs[str(record.id)]=str(record.seq)

    sf_flank_alnmt="velvet_temp/{0}/flanks.sam".format(id)
    m_hit_id={}
    m_picked_id={}
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

    #print m_hit_id ####################################################################################################
    #print "111111111111111111"
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

        # if m_hit_id.has_key(id+"_left")==False or m_hit_id.has_key(id+"_right")==False:
        #     return
        # if m_hit_id[id+"_left"]["left"].has_key(s_left_picked)==False \
        #         or m_hit_id[id+"_right"]["right"].has_key(s_right_picked)==False:
        #     return

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
        # if m_hit_id.has_key(id+"_left")==False or m_hit_id.has_key(id+"_right")==False:
        #     return
        # if m_hit_id[id+"_left"]["left"].has_key(s_left_picked)==False \
        #         or m_hit_id[id+"_right"]["right"].has_key(s_right_picked)==False:
        #     return

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

    #print "33333333333333333333333" #################################################################################3

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
        sf_picked="velvet_temp/{0}/picked_seqs.fa".format(id)
        with open(sf_picked,"w") as fout_picked:
            fout_picked.write(">"+id+"_"+s_left_picked+"_"+s_right_picked+"_extended"+"\n")
            fout_picked.write(s_seq+"\n")

    if s_contig!="" and s_contig!="NN":
        sf_contig="velvet_temp/{0}/picked_contigs.fa".format(id)
        with open(sf_contig,"w") as fout_picked:
            fout_picked.write(">"+id+"_"+s_left_picked+"_"+s_right_picked+"_extended"+"\n")
            fout_picked.write(s_contig+"\n")


def pick_extended_contigs(n, fa_list, sf_picked):
    pool = Pool(n)
    pool.map(run_pick_extended_contig, fa_list, bunch_size)
    pool.close()
    pool.join()

    for key in fa_list:
        sf_tmp="velvet_temp/{0}/picked_seqs.fa".format(key)
        if os.path.exists(sf_tmp)==True:
            cmd="cat {0} >> {1}".format(sf_tmp, sf_picked)
            print cmd
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        sf_tmp="velvet_temp/{0}/picked_contigs.fa".format(key)
        if os.path.exists(sf_tmp)==True:
            cmd="cat {0} >> {1}".format(sf_tmp, sf_picked+"_ori.txt")
            print cmd
            Popen(cmd, shell = True, stdout = PIPE).communicate()


def validate_results(sf_fai, sf_gap_pos, n):
    fa_list=prepare_list(sf_fai, sf_gap_pos)
    pool = Pool(n)
    pool.map(run_alignment_contig, fa_list, bunch_size)
    pool.map(run_alignment_reads, fa_list, bunch_size)
    pool.close()
    pool.join()

    sf_merged_contig="merged_contigs_alignment.sam"
    if os.path.exists(sf_merged_contig)==True:
        cmd="rm {0}".format(sf_merged_contig)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    sf_merged_reads="merged_gap_reads_alignment.sam"
    if os.path.exists(sf_merged_reads)==True:
        cmd="rm {0}".format(sf_merged_reads)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    for id in fa_list:
        sf="merged_contigs_alignment/{0}.sort.bam".format(id)
        if os.path.exists(sf):
            cmd="samtools view -F 4 {0} >> {1}".format(sf, sf_merged_contig)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        # sf="gap_reads_alignment/{0}.sort.bam".format(id)
        # if os.path.exists(sf):
        #     cmd="samtools view -F 4 {0} >> {1}".format(sf, sf_merged_reads)
        #     Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="cat chrom14_gap_seqs_header.txt merged_contigs_alignment.sam | samtools view -S -h -b - | " \
        "samtools sort - -o merged_contigs_alignment.sort.bam"
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools index merged_contigs_alignment.sort.bam"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    # cmd="cat chrom14_gap_seqs_header.txt merged_gap_reads_alignment.sam | samtools view -S -h -b - | " \
    #     "samtools sort - -o merged_gap_reads_alignment.sort.bam"
    # Popen(cmd, shell = True, stdout = PIPE).communicate()
    # cmd="samtools index merged_gap_reads_alignment.sort.bam"
    # Popen(cmd, shell = True, stdout = PIPE).communicate()


def run_contigs_assembly_test():
    global m_kmer_len

    fq_list=[]
    fq_list.append("0_729")

    for klen in m_kmer_len:
        klen_fields=klen.split("_")

        kmer_len=int(klen_fields[0])
        asm_len=int(klen_fields[1])

        pool = ThreadPool(1)
        pool.map(run_assembly, fq_list)
        pool.close()
        pool.join()


def run_contigs_merge_test():
    fa_list=[]
    fa_list.append("0_729")

    pool = ThreadPool(1)
    pool.map(run_merge, fa_list)
    pool.close()
    pool.join()

def initial():
    global m_kmer_len
    m_kmer_len.append("31_29")
    m_kmer_len.append("41_37")
    m_kmer_len.append("41_39")
    m_kmer_len.append("51_47")
    m_kmer_len.append("61_57")


if __name__ == "__main__":
    initial()
    sf_fai=sys.argv[1]
    sf_pos=sys.argv[2]
    n_jobs=int(sys.argv[3])

    l_temp={}
    fa_list=prepare_list(sf_fai, sf_pos)
    bwa_min_score=30
    # #
    # print "First round assembly and merger..."
    assembly(fa_list, n_jobs)
    run_contigs_merge(n_jobs, fa_list)
    # #
    sf_picked="picked_seqs.fa"
    # bwa_min_score=30######################################################################################

    #print fa_list
    pick_full_constructed_contigs(bwa_min_score, n_jobs, fa_list, sf_picked)
    m_picked=get_already_picked(sf_picked)

    id_remain=[]
    for key in fa_list:#
        if m_picked.has_key(key)==False:
            id_remain.append(key)

    bam_list=[]
    bam_list.append("./../insert_size_500/align_2_ref.sort.bam")
    bam_list.append("./../insert_size_750/align_2_ref.sort.bam")
    # bam_list.append("./../long_jump_nx/align_2_ref.sort.bam")
    # bam_list.append("./../short_insert/align_2_ref.sort.bam")
    # bam_list.append("./../short_jump/align_2_ref.sort.bam")

    print "Collect both unmapped reads..."
    collect_both_unmapped_reads(bam_list, id_remain, n_jobs)
    print "Second round assembly..."
    assembly_given_list(id_remain,n_jobs)

    pick_full_constructed_contigs(bwa_min_score, n_jobs, id_remain, sf_picked)
    m_picked=get_already_picked(sf_picked)
    #
    del id_remain[:]
    for key in fa_list:
        if m_picked.has_key(key)==False:
            id_remain.append(key)

    #print "Second round merging..."
    run_contigs_merge(n_jobs, id_remain)

    ##
    pick_full_constructed_contigs(bwa_min_score, n_jobs, id_remain, sf_picked)
    m_picked=get_already_picked(sf_picked)
    #
    del id_remain[:]
    for key in fa_list:
        if m_picked.has_key(key)==False:
            id_remain.append(key)
    #
    collect_high_quality_unmap_to_contigs_reads(n_jobs, id_remain)
    run_contigs_merge(n_jobs, id_remain)
    # #
    # # # # #
    bwa_min_score=15
    pick_full_constructed_contigs(bwa_min_score, n_jobs, id_remain, sf_picked)
    #
    m_picked=get_already_picked(sf_picked)
    id_remain=[]
    for key in fa_list:
        if m_picked.has_key(key)==False:
            id_remain.append(key)
    pick_extended_contigs(n_jobs, id_remain, sf_picked)
