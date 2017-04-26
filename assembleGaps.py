import sys
import os
from subprocess import *
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from MergeContigs import *
from Bio import SeqIO
from Utility import *
from pick_contigs import ContigsSelection
from collect_both_unmapped_reads import BothUnmappedReadsCollector

kmer_len_list=[]
working_folder=""
kmc_path=""
velvet_path=""
bwa_path=""
samtools_path=""

def run_clear(id):
    global working_folder
    sf_ori_contig=working_folder+"velvet_temp/{0}/original_contigs_before_merging.fa".format(id)
    if os.path.exists(sf_ori_contig):
        print id

    sf_contig0=working_folder+"velvet_temp/{0}/contigs.fa".format(id)
    sf_contig1=working_folder+"velvet_temp/{0}/contigs_31_29.fa".format(id)
    sf_contig2=working_folder+"velvet_temp/{0}/contigs_41_37.fa".format(id)
    sf_contig3=working_folder+"velvet_temp/{0}/contigs_41_39.fa".format(id)
    sf_contig5=working_folder+"velvet_temp/{0}/contigs_51_47.fa".format(id)
    sf_contig6=working_folder+"velvet_temp/{0}/contigs_61_57.fa".format(id)
    sf_seq=working_folder+"velvet_temp/{0}/Sequences".format(id)
    if os.path.exists(sf_contig0) and os.path.exists(sf_contig1) and os.path.exists(sf_contig2) \
            and os.path.exists(sf_contig3) and os.path.exists(sf_contig5) and os.path.exists(sf_contig6)\
            and os.path.exists(sf_seq):
        cmd="rm {0}temp/{1}.*".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="rm {0}velvet_temp/{1}/Sequences".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm {0}velvet_temp/{1}/Roadmaps".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm {0}velvet_temp/{1}/PreGraph".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm {0}velvet_temp/{1}/Graph".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm {0}velvet_temp/{1}/LastGraph".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm {0}velvet_temp/{1}/Log".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="rm {0}velvet_temp/{1}/stats.txt".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="rm {0}kmers/{0}_*".format(working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

def run_assembly(id):
    global kmer_len_list
    global working_folder

    for klen in kmer_len_list:
        klen_fields=klen.split("_")
        kmer_len=int(klen_fields[0])
        asm_len=int(klen_fields[1])

        pth=working_folder+"kmc_temp/{0}".format(id)
        if os.path.exists(pth)==False:
            cmd="mkdir {0}kmc_temp/{1}".format(working_folder,id)
            #print cmd
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="{0}kmc -k{1} -cs10000000 -m52 {2}gap_reads/{3}.fastq {4}temp/{5}.res {6}kmc_temp/{7}"\
            .format(kmc_path, kmer_len,working_folder,id,working_folder,id,working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="{0}kmc_dump -ci0 {1}temp/{2}.res {3}temp/{4}.dump"\
            .format(kmc_path, working_folder, id, working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="awk -f {0}../cvtKMC_2_Fq.awk {1}temp/{2}.dump > {3}kmers/{4}_kmers.fq"\
            .format(working_folder,working_folder,id,working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        pth=working_folder+"velvet_temp/{0}".format(id)
        if os.path.exists(pth)==False:
            cmd="mkdir {0}velvet_temp/{1}".format(working_folder,id)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="{0}velveth {1}velvet_temp/{2} {3} -fastq -short {4}kmers/{5}_kmers.fq"\
            .format(velvet_path, working_folder,id, asm_len, working_folder,id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="{0}velvetg {1}velvet_temp/{2} -min_contig_lgth 40".format(velvet_path, working_folder, id)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="mv {0}velvet_temp/{1}/contigs.fa {2}velvet_temp/{3}/contigs_{4}_{5}.fa"\
            .format(working_folder,id,working_folder,id,kmer_len, asm_len)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    sf_merged=working_folder+"velvet_temp/{0}/contigs.fa".format(id)
    with open(sf_merged,"w") as fout_merged:
        for klen in l_kmer_len:
            klen_fields=klen.split("_")
            ctg=working_folder+"velvet_temp/{0}/contigs_{1}_{2}.fa".format(id, klen_fields[0], klen_fields[1])
            with open(ctg) as fin_ctg:
                for line in fin_ctg:
                    if line[0]==">":
                        fout_merged.write(">"+klen+"_"+line[1:])
                    else:
                        fout_merged.write(line)
    run_clear(id)

def run_merge(sf_contig):
    global working_folder
    #merge contigs
    sf_folder=working_folder+"velvet_temp/{0}/".format(sf_contig)
    sf_contig=sf_folder+"contigs.fa"
    if os.path.exists(sf_contig)==False:
        return
    merge_contigs(sf_folder, 5, 0.99, 0.99)

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
    global working_folder
    global bwa_path
    global samtools_path

    sf_contig=working_folder+"velvet_temp/"+id+"/contigs.fa"
    #print sf_contig##########################################################################################
    # if os.path.exists(sf_contig+".bwt")==False:
    cmd="{0} index {1}".format(bwa_path, sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    #print cmd##############################################################################################33
    sf_reads=working_folder+"gap_reads_high_quality/{0}.fastq".format(id)
    if os.path.exists(sf_reads)==False:
        return
    sf_algmt=working_folder+"gap_reads_high_quality/{0}.sam".format(id)
    cmd="{0} mem -t 5 {1} {2} | {3} view -S - > {4}".format(bwa_path, sf_contig, sf_reads, samtools_path, sf_algmt)
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

    sf_ori_contig=working_folder+"velvet_temp/"+id+"/original_contigs_before_merging.fa"
    cmd="rm {0}".format(sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mv {0} {1}".format(sf_ori_contig, sf_contig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    with open(sf_contig,"a") as fout_picked:
        for read_id in m_read_ref:
            if len(m_read_ref[read_id])>=2: ##at least clipped at two contigs
                fout_picked.write(">"+read_id+"\n")
                fout_picked.write(m_read_seq[read_id]+"\n")


class GapAssembler():
    def __init__(self, sf_fai, sf_pos, n_jobs, working_space):
        global l_kmer_len
        l_kmer_len=get_kmer_list()
        self.sf_fai=sf_fai
        self.sf_pos=sf_pos
        self.n_jobs=int(n_jobs)
        global working_folder
        working_folder=working_space

        ##
        global velvet_path
        velvet_path=get_velvet_path()#"/gpfs/scratchfs1/chc12015/tools/velvet-master/" #/data2/chongchu/tools/velvet-master/
        global kmc_path
        kmc_path=get_kmc_path()#"/gpfs/scratchfs1/chc12015/tools/kmc2.3/KMC/bin/" #/data2/chongchu/tools/kmc2.3/
        global samtools_path
        samtools_path=get_samtools_path()
        global bwa_path
        bwa_path=get_bwa_path()
        self.bunch_size=1
        self.refiner_path=get_refiner_path()
        self.merger_path=get_merger_path()


    def prepare_list(self):
        global working_folder
        m_scaffold_id={}
        cnt=0
        with open(self.sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1

        fa_list=[]
        cnt=1
        pre_id=0

        with open(self.sf_gap_pos) as fin_gap_pos:#for each gap
            for line in fin_gap_pos:
                fields=line.split()
                scaffold=fields[3]
                scaffold_id=m_scaffold_id[scaffold]

                if pre_id!=scaffold_id:
                    cnt=1

                sf_fa="{0}_{1}".format(scaffold_id,cnt)
                cnt=cnt+1

                pre_id=scaffold_id
                sf_fq="{0}gap_reads/{1}.fastq".format(working_folder,sf_fa)
                if os.path.exists(sf_fq):
                    fa_list.append(sf_fa)
        return fa_list

    def assembly(self, id_list):
        global working_folder
        if os.path.exists(working_folder+"empty_dir")==False:
            cmd="mkdir {0}empty_dir".format(working_folder)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        ##clear old data
        cmd="rsync -a --delete {0}empty_dir/ {1}kmc_temp/".format(working_folder,working_folder)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="rsync -a --delete {0}empty_dir/ {1}temp/".format(working_folder,working_folder)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="rsync -a --delete {0}empty_dir/ {1}kmers/".format(working_folder,working_folder)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        cmd="rsync -a --delete {0}empty_dir/ {1}velvet_temp/".format(working_folder,working_folder)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        pool = Pool(self.n_jobs)
        pool.map(run_assembly, id_list, self.bunch_size)
        pool.close()
        pool.join()

    def run_contigs_merge(self, fa_list):
        print "Start merging..."
        pool=Pool(self.n_jobs)
        pool.map(run_merge, fa_list, self.bunch_size)
        pool.close()
        pool.join()

    def assembly_given_list(self, id_list):
        pool = Pool(self.n_jobs)
        pool.map(run_assembly, id_list, self.bunch_size)
        pool.close()
        pool.join()

    def collect_high_quality_unmap_to_contigs_reads(self, id_list):
        pool = Pool(self.n_jobs)
        pool.map(run_collect_high_quality_unmap_to_contig_reads, id_list, self.bunch_size)
        pool.close()
        pool.join()

    def pick_already_constructed(self, contigs_select, bwa_min_score, fa_list, sf_picked):
        contigs_select.pick_full_constructed_contigs(bwa_min_score, fa_list, sf_picked)
        m_picked=contigs_select.get_already_picked(sf_picked)

        id_remain=[]
        for key in fa_list:#
            if m_picked.has_key(key)==False:
                id_remain.append(key)
        return id_remain

    def assemble_pipeline(self):
        global working_folder
        fa_list=self.prepare_list()
        print "First round assembly and merger..."
        self.assembly(fa_list)
        self.run_contigs_merge(fa_list)

        sf_picked=working_folder+"picked_seqs.fa"
        bwa_min_score=30
        contigs_select=ContigsSelection(working_folder)

        id_remain=self.pick_already_constructed(contigs_select, bwa_min_score, fa_list, sf_picked)
        algnmt_list=get_alignment_list()
        bam_list=[]
        for algnmt in algnmt_list:
            bam_list.append(algnmt[0])

        print "Collect both unmapped reads..."
        burc=BothUnmappedReadsCollector(working_folder)
        burc.collect_both_unmapped_reads(bam_list, id_remain)

        print "Second round assembly..."
        self.assembly_given_list(id_remain)
        id_remain2=self.pick_already_constructed(contigs_select, bwa_min_score, id_remain, sf_picked)

        print "Second round merging..."
        self.run_contigs_merge(id_remain2)
        id_remain3=self.pick_already_constructed(contigs_select, bwa_min_score, id_remain2, sf_picked)

        print "Collecting high quality reads to improve the merging step..."
        self.collect_high_quality_unmap_to_contigs_reads(id_remain3)
        self.run_contigs_merge(id_remain3)

        # # # # # #
        print "Pick extended gap sequences..."
        bwa_min_score=15
        id_remain4=self.pick_already_constructed(contigs_select, bwa_min_score, id_remain3, sf_picked)
        contigs_select.pick_extended_contigs(bwa_min_score, id_remain4, sf_picked)

# def run_filter(id):
#     fcontig="velvet_temp/"+id+"/contigs.fa"
#     cmd='bwa index {0}'.format(fcontig)
#     Popen(cmd, shell = True, stdout = PIPE).communicate()
#     #align reads to contigs
#     sf_reads="gap_reads_alignment/"+id+".fastq"
#     sf_sorted_bam="velvet_temp/"+id+"/alignment_for_coverage"
#     cmd="bwa mem -a -t 5 {0} {1} | samtools view -h -S -b - | " \
#         "samtools sort - -o {2}.sort.bam".format(fcontig, sf_reads,sf_sorted_bam)
#     Popen(cmd, shell = True, stdout = PIPE).communicate()
#
#     cmd="samtools faidx {0}".format(fcontig)
#     Popen(cmd, shell = True, stdout = PIPE).communicate()
#
#     #calculate coverage cover
#     cmd="TERefiner_1 -G -r {0} -b {1} -c 0"
#     Popen(cmd, shell = True, stdout = PIPE).communicate()
#     #"_cov_info_with_cutoff.txt"
#
#
# def filter_contigs(sf_fai, sf_gap_pos, n):
#     #align the reads back to contigs
#     #calculate the average coverage
#     #output the filtered contigs
#     fa_list=prepare_list(sf_fai, sf_gap_pos)
#
#     pool = Pool(n)
#     pool.map(run_merge, fa_list,self.bunch_size)
#     pool.close()
#     pool.join()
#
#     pool1 = Pool(n)
#     pool1.map(run_filter, fa_list,self.bunch_size)
#     pool1.close()
#     pool1.join()
