import sys
import os
import argparse
from subprocess import *
import json
from Utility import *
from gnrt_pos_true_seqs import DGProcessor
from run_multi_threads_collect_reads import MultiThrdReadsCollector
from run_multi_threads_discordant import DiscordantReadsCollector
from merge_reads import ReadsMerger
from assemble_gaps import GapAssembler

MERGE_FOLDER="merged/"

def usage():
    print 'Usage: python {0} -c Options -g configure_file\n'.format(sys.argv[0]),
    print 'Options:\n',
    print '    Clean          Remove all the temporary files\n',
    print '    All            Run the whole pipeline\n',
    print '    Preprocess     Get gap positions from draft genomes\n',
    print '    Collect        Collect and merge the related reads from alignments\n',
    print '    Assembly       Assemble the collect reads (including collect both unmapped reads)\n',
    print 'Example of running the whole pipeline: python main.py -c All -g config.txt'
########################################################################################################################

def get_args():
    # Assign description to the help doc
    parser = argparse.ArgumentParser(
        description='Run the pipeline of GAPPadder')
    # Add arguments
    parser.add_argument(
        '-g', '--config', type=str, help='Configuration file name', required=True)
    parser.add_argument(
        '-c', '--command', type=str, help='Specific command', required=True)

    args = parser.parse_args()
    sfconfig=args.config
    scommand=args.command
    return sfconfig, scommand


def parse_configuration(sf_config):
    sf_draft, min_gap_len, flank_len, nthreads, working_folder, vbs="","","","","",""
    pbwa, psamtools, pvelvet, pkmc, prefiner, pmerger="","","","","",""
    lalgnmt=[]
    lkmers=[]
    lraw_reads=[]
    with open(sf_config) as data_file:
        data = json.load(data_file)

        ##first check whether required settings are set
        if "draft_genome"not in data or "alignments" not in data:
            return False
        ##required settings
        #input draft genome
        sf_draft=data["draft_genome"]["fa"]
        ##alignment list
        for record in data["alignments"]:
            bam=record["bam"]
            insert_size=record["is"]
            std_derivation=record["std"]
            lalgnmt.append((bam,insert_size,std_derivation))

        ##optinal settings
        #parameters
        if "parameters" in data:
            min_gap_len=data["parameters"]["min_gap_size"]
            flank_len=data["parameters"]["flank_length"]
            nthreads=data["parameters"]["nthreads"]
            working_folder=data["parameters"]["working_folder"]
            vbs=data["parameters"]["verbose"]

        #raw reads
        for record in data["raw_reads"]:
            left_fq=record["left"]
            rigth_fq=record["right"]
            lraw_reads.append((left_fq, rigth_fq))

        ##kmer list
        if "kmer_length" in data:
            for record in data["kmer_length"]:
                k=record["k"]
                for subrecord in record["k_velvet"]:
                    sub_k=subrecord["k"]
                    lkmers.append((k,sub_k))

        #software paths
        if "software_path" in data:
            pbwa=data["software_path"]["bwa"]
            psamtools=data["software_path"]["samtools"]
            pvelvet=data["software_path"]["velvet"]
            prefiner=data["software_path"]["TERefiner"]
            pmerger=data["software_path"]["ContigsMerger"]

    set_software_paths(pbwa, psamtools, prefiner, pmerger, pkmc, pvelvet)
    #print flank_len###############################################################################################
    set_parameters(sf_draft, nthreads, working_folder , vbs, min_gap_len, flank_len)
    set_alignment_list(lalgnmt)
    set_kmer_list(lkmers)
    return True


def prepare_sub_folders(spath):
    cmd="mkdir {0}/scaffold_reads_list_all".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/gap_reads".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/gap_reads_for_alignment".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/gap_reads_high_quality".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/discordant_reads_list".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/discordant_temp".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def prepare_sub_folders_merged(spath):
    if os.path.exists("{0}/kmc_temp".format(spath))==True:
        return
    cmd="mkdir {0}/gap_reads".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/gap_reads_for_alignment".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/gap_reads_high_quality".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="mkdir {0}/kmc_temp".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/temp".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/kmers".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/velvet_temp".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/both_unmapped".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="mkdir {0}/unmapped_reads".format(spath)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def prepare_folders(algmt_list, working_space):
    global MERGE_FOLDER
    #first prepare folders
    cnt=1
    folder_list=[]
    for algmt in algmt_list:
        folder="{0}_is{1}".format(cnt, algmt[1])
        folder_list.append(folder)
        spath=working_space+folder
        if os.path.exists(spath)==False:
            cmd="mkdir {0}".format(spath)
            Popen(cmd, shell = True, stdout = PIPE).communicate()
            if os.path.exists("{0}/scaffold_reads_list_all".format(spath)):
                continue
            prepare_sub_folders(spath)

    #folder="merged" #prepare for the merged folder
    folder=MERGE_FOLDER
    spath=working_space+folder
    prepare_sub_folders_merged(spath)
    if os.path.exists(spath)==False:
        cmd="mkdir {0}".format(spath)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        if os.path.exists("{0}/scaffold_reads_list_all".format(spath))==False:
            prepare_sub_folders(spath)
    return folder_list

def clean_all(working_folder):
    empty_folder="empty_dir"
    if os.path.exists(empty_folder):
        cmd="mkdir {0}".format(empty_folder)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="rsync -a --delete {0}/ {1}".format(empty_folder, working_folder)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def main_func(scommand, sf_config):
    global MERGE_FOLDER
    brtn=parse_configuration(sf_config)
    if brtn==False:
        print "The required parameters are not set in the configuration file!!!"
        return

    min_gap_lenth=get_min_gap_length()
    flank_length=get_flank_length()
    sf_draft=get_draft_genome()
    sf_fai=sf_draft+".fai"
    samtools_path=get_samtools_path()
    if os.path.exists(sf_fai)==False:
        cmd="{0} faidx {1}".format(samtools_path, sf_draft)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    working_folder=get_output_folder()
    if working_folder[-1]!="/":
        working_folder=working_folder+"/"
    sf_gap_pos=working_folder+"gap_positions.txt"
    anchor_mapq=30
    clip_dist=250
    nthreads=get_threads_num()

    algmt_list=get_alignment_list()
    if scommand=="Clean" or scommand=="All":
        clean_all(working_folder)
    if scommand=="Preprocess" or scommand=="All":
        dgp=DGProcessor(sf_draft, sf_gap_pos)
        dgp.gnrt_gap_positions(min_gap_lenth, sf_gap_pos)
        dgp.get_gap_flank_seqs(sf_draft, sf_gap_pos, flank_length, sf_fai)
    if scommand=="Collect" or scommand=="All":
        folder_list=prepare_folders(algmt_list, working_folder)#first collect reads for each alignment
        raw_reads_list=get_raw_reads_list()
        for algmt, folder, raw_reads in algmt_list, folder_list, raw_reads_list:
            sf_bam=algmt[0]
            insert_size=algmt[1]
            derivation=algmt[2]
            mtrc=MultiThrdReadsCollector(sf_fai, sf_bam, sf_gap_pos, anchor_mapq)
            mtrc.dispath_collect_jobs(nthreads, samtools_path, insert_size, derivation, clip_dist, folder)

            ##run reads collect for discordant
            left_reads=raw_reads[0]
            right_reads=raw_reads[1]
            if folder[-1]!="/":
                folder=folder+"/"
            sfout_discord_pos=folder+"discordant_reads_list"

            drc=DiscordantReadsCollector(sf_fai, sf_bam, folder, nthreads)
            drc.collect_discordant_regions_v2(sfout_discord_pos)
            drc.dispath_collect_jobs()

            drc.merge_dispatch_reads_for_gaps_v2(left_reads, right_reads)
            #dispatch_reads_for_gaps_to_validate_contigs(left_reads, right_reads)
            drc.dispatch_high_quality_reads_for_gaps(left_reads, right_reads)
            #####################################

        #then merge the reads
        rmerger=ReadsMerger()
        rmerger.merge_reads_v2(sf_fai, sf_gap_pos, folder_list, "gap_reads", working_folder+MERGE_FOLDER, nthreads)
        rmerger.merge_reads_v2(sf_fai, sf_gap_pos, folder_list, "gap_reads_alignment", working_folder+MERGE_FOLDER, nthreads)
        rmerger.merge_reads_v2(sf_fai, sf_gap_pos, folder_list, "gap_reads_high_quality", working_folder+MERGE_FOLDER, nthreads)

        ##here remove the temporary files???? ##########################################################################
        for folder in folder_list:
            clean_all(folder)

    if scommand=="Assembly" or scommand=="All":
        gap_assembler=GapAssembler(sf_fai, sf_gap_pos, nthreads, working_folder+MERGE_FOLDER)
        gap_assembler.assemble_pipeline()

    return


if __name__ == "__main__":
    if len(sys.argv) <= 2:
        usage()
        raise SystemExit

    sfconfig, scommand = get_args()
    main_func(scommand,sfconfig)
