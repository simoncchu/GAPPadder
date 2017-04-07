import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# handle=open("test00.fa","w")
# for record in SeqIO.parse(sys.argv[1], "fasta"):
#     #record.seq[0]="N"
#     #seq=str(record.seq)
#     seq1="N"
#     new_r=SeqRecord(Seq(seq1)+record.seq[0:10], id = record.id, description="")
#     SeqIO.write(new_r,handle,"fasta")
# handle.close()

