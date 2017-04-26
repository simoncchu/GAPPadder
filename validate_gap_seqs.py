def run_alignment_reads(id):
    sf_ref=self.working_folder+"../gap_seqs/{0}.fa".format(id)
    if os.path.exists(sf_ref)==True:
        sf_bwt=self.working_folder+"../gap_seqs/{0}.fa.bwt".format(id)
        if os.path.exists(sf_bwt)==False:
            cmd="bwa index {0}".format(sf_ref)
            Popen(cmd, shell = True, stdout = PIPE).communicate()
    else:
        return

    sf=self.working_folder+"gap_reads/{0}.fastq".format(id)
    if os.path.exists(sf)==True:
        cmd="bwa mem {0} {1} | " \
            "samtools view -S -h -b - | samtools sort - -o gap_reads_alignment/{2}.sort.bam " \
            "&& samtools index gap_reads_alignment/{3}.sort.bam".format(sf_ref,sf,id,id)
        #print cmd#################
        Popen(cmd, shell = True, stdout = PIPE).communicate()

def run_alignment_contig(id):
    sf_ref=self.working_folder+"../gap_seqs/{0}.fa".format(id)
    if os.path.exists(sf_ref)==True:
        sf_bwt=self.working_folder+"../gap_seqs/{0}.fa.bwt".format(id)
        if os.path.exists(sf_bwt)==False:
            cmd="bwa index {0}".format(sf_ref)
            Popen(cmd, shell = True, stdout = PIPE).communicate()
    else:
        return

    sf=self.working_folder+"velvet_temp/{0}/contigs.fa".format(id)
    if os.path.exists(sf)==True:
        cmd="bwa mem {0} {1} | " \
            "samtools view -S -h -b - | samtools sort - -o merged_contigs_alignment/{2}.sort.bam " \
            "&& samtools index merged_contigs_alignment/{3}.sort.bam".format(sf_ref,sf,id,id)
        #print cmd#################
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
