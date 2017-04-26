import sys
import os

def readContigFa(sffa,brstrip):
    dcontigs={}
    name=""
    seq=""
    last_name=""
    ffa=open(sffa)
    for line in ffa:
        if brstrip==True:
            line=line.rstrip()
        if line[0]=='>':
            fileds=line.split()

            name=fileds[0][1:]
            name=name.rstrip()
            if last_name!="":
                dcontigs[last_name]=seq
                seq=""
            last_name=name
        else:
            seq=seq+line
    dcontigs[last_name]=seq
    ffa.close()
    return dcontigs

dcontigs=readContigFa(sys.argv[1],True)

dqualified={}
with open(sys.argv[2]) as finBlast:
    for line in finBlast:
        fields=line.split()
        if fields[0]==fields[1]:
            continue
        contigLen_que=len(dcontigs[fields[0]])
        #contigLen_sub=len(dcontigs[fields[1]])
        hitLen=int(fields[7])-int(fields[6])

        key=fields[0]
        if float(hitLen)/float(contigLen_que) >= float(sys.argv[3]):
            dqualified[key]=1

for key in dqualified:
    print ">"+key
    print dcontigs[key].rstrip()
