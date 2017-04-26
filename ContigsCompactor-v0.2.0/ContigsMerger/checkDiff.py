import sys

dv1={}
with open(sys.argv[1]) as finV1:
    for line in finV1:
        dv1[line.rstrip()]=1

with open(sys.argv[2]) as finV0:
    for line in finV0:
        if dv1.has_key(line.rstrip())==False:
            print line.rstrip()



# s1="GAGGTCTGAATATCCACTTGCAGACTTTACAAACAGAGTGTTTCCTAACTGCTCTATGAAAAGAAAGGTTAAACTCTGTGAGTTGAACGCACACATCACAAAGGAGTTTCTGAGAATCATTCTGTCTAGTTT"
# s2="TCTCAGCCCAAAATCTCCTTAAGGTGATAAGCAACTTCAGCAAAGTCTCAGGAAACAAAATCAATGTGCAAAAATCACAAGCACTCTTATACACCAATAACAGACAAACAGAG"
#
# dS1={}
# k=10
# i=0
# while(i<len(s1)-k+1):
#     dS1[ s1[i:i+k] ]=1
#     i=i+1
#
# i=0
# while(i<len(s2)-k+1):
#     if dS1.has_key(s2[i:i+k]):
#         print s2[i:i+k]
#         print "hit"
#         break
#     i=i+1
#
