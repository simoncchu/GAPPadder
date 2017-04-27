import sys
from subprocess import *
from multiprocessing.dummy import Pool as ThreadPool

def run_cmd_reads_collect(cmd1):
    print cmd1
    check_output(cmd1, shell=True)

class MultiThrdReadsCollector:
    def __init__(self, sf_fai, sf_bam, sf_gap_pos, anchor_mapq):
        self.sf_fai=sf_fai
        self.sf_bam=sf_bam
        self.sf_gap_pos=sf_gap_pos
        self.anchor_mapq=anchor_mapq

    def dispath_collect_jobs(self, nthreads, samtools_path, insert_size, derivation, clip_dist, working_folder):
        m_scaffold_has_gap_list={}
        with open(self.sf_gap_pos) as fin_gap_pos:
            for line in fin_gap_pos:
                fields=line.split()
                m_scaffold_has_gap_list[fields[3]]=1

        cmd_list=[]
        with open(self.sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                if m_scaffold_has_gap_list.has_key(fields[0])==False:
                    continue

                cmd="{0} view {1} \"{2}\" | python collect_reads_for_gaps.py {3} {4} {5} {6} {7} {8} -"\
                    .format(samtools_path, self.sf_bam, fields[0], self.sf_gap_pos, self.anchor_mapq, working_folder,
                            insert_size, derivation, clip_dist)
                cmd_list.append(cmd)

        pool = ThreadPool(nthreads)
        pool.map(run_cmd_reads_collect, cmd_list)
        pool.close()
        pool.join()

# if __name__ == "__main__":
#     sf_fai=sys.argv[1]
#     sf_bam=sys.argv[2]
#     sf_gap_pos=sys.argv[3]
#     anchor_mapq=int(sys.argv[4])
#
#     dispath_collect_jobs(int(sys.argv[5]))
