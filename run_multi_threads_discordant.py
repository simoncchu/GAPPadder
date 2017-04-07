import sys
import os
from subprocess import *
from multiprocessing.dummy import Pool as ThreadPool

sf_bam=""
sf_fai=""
sf_discordant_pos=""


##In this version, allow one regions linked to more than 1 gaps
def collect_discordant_regions_v2(sf_fai,sf_out):
    m_scaffold_id={}
    l_scaffold_id=[]
    cnt=0##scaffold id, started from 0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1
                l_scaffold_id.append(scaffold_id)

    ##structure: scaffold-Id-Of-Discordant-Region_pos [linked-scaffold-id_gap-id]
    m_collected_pos={}

    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            scaffold_id=fields[0]

            i_id=m_scaffold_id[scaffold_id]

            sf_scaffold_right_list="./scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(scaffold_id)
            sf_scaffold_left_list="./scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(scaffold_id)

            if os.path.exists(sf_scaffold_left_list)==False or os.path.exists(sf_scaffold_right_list)==False:
                continue

            with open(sf_scaffold_left_list) as fin_left_list:
                for record in fin_left_list:
                    info_fields=record.split()
                    if info_fields[3]=="discordant":
                        #mate_scaffold  pos_mate_scaffold  current_scaffold  current_gap_id
                        mate_scaffold_id=info_fields[5]
                        if mate_scaffold_id=="=":
                            mate_scaffold_id = scaffold_id
                        i_mate_id=m_scaffold_id[mate_scaffold_id]
                        mate_pos=info_fields[6]

                        key="{0}_{1}".format(i_mate_id, mate_pos)
                        if m_collected_pos.has_key(key)==False:
                            m_collected_pos[key]=[]

                        gap_id=info_fields[1]
                        linked_gap="{0}_{1}".format(i_id,gap_id)
                        m_collected_pos[key].append(linked_gap)

                        # fout_discordant.write(mate_scaffold_id+" "+ mate_pos +" "
                        #                       +scaffold_id+" "+info_fields[1]+"\n")

            with open(sf_scaffold_right_list) as fin_right_list:
                for record in fin_right_list:
                    info_fields=record.split()
                    if info_fields[3]=="discordant":
                        mate_scaffold_id=info_fields[5]
                        if mate_scaffold_id=="=":
                            mate_scaffold_id = scaffold_id
                        i_mate_id=m_scaffold_id[mate_scaffold_id]
                        mate_pos=info_fields[6]

                        key="{0}_{1}".format(i_mate_id, mate_pos)
                        if m_collected_pos.has_key(key)==False:
                            m_collected_pos[key]=[]

                        gap_id=info_fields[1]
                        linked_gap="{0}_{1}".format(i_id,gap_id)
                        m_collected_pos[key].append(linked_gap)
                        #mate_scaffold  pos_mate_scaffold  current_scaffold  current_gap_id


    with open(sf_out,"w") as fout_discordant:
        for key in m_collected_pos:
            fields=key.split("_")
            m_id=fields[0]
            m_pos=fields[1]

            for content in m_collected_pos[key]:
                fld=content.split("_")
                s_id=fld[0]
                g_id=fld[1]
                fout_discordant.write(m_id+" "+ m_pos +" "+s_id+" "+g_id+"\n")

    m_collected_pos.clear()

    # ##sort by gap first, and then filter out records with threshold
    # sorted_by_gap=sf_out+".sorted_by_gap.txt"
    # cmd="sort -k3n -k4n -k1n -k2n {0} > {1}".format(sf_out,sorted_by_gap)
    # Popen(cmd, shell = True, stdout = PIPE).communicate()
    #
    # with open(sorted_by_gap) as fin_sorted_by_gap:
    #     for line in fin_sorted_by_gap:
    #         fields_gap=line.split()

    ##sort by map-to-scaffold, and split
    sorted=sf_out+".sorted.txt"
    cmd="sort -k1n -k2n -k3n -k4n {0} > {1}".format(sf_out, sorted)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #split
    pre=""
    ftemp=None
    with open(sorted) as fin_sorted:
        for line in fin_sorted:
            fields=line.split()
            id=int(fields[0])
            name=l_scaffold_id[id]
            if pre=="":
                ftemp=open("./discordant_temp/"+name+".list","w")
            elif name!=pre:
                ftemp.close()
                ftemp=open("./discordant_temp/"+name+".list","w")
            ftemp.write(line)
            pre=name
        ftemp.close()


def split_by_gap(sf_fai, sf_sorted):
    l_scaffold_id=[]
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                l_scaffold_id.append(scaffold_id)
    pre=""
    ftemp=None
    with open(sf_sorted) as fin_sorted:
        for line in fin_sorted:
            fields=line.split()
            id=int(fields[2])
            name=l_scaffold_id[id]
            if pre=="":
                ftemp=open("./discordant_temp_by_gap/"+name+".by_gap.list","w")
            elif name!=pre:
                ftemp.close()
                ftemp=open("./discordant_temp_by_gap/"+name+".by_gap.list","w")
            ftemp.write(line)
            pre=name
    ftemp.close()

def run_cmd(cmd1):
    print cmd1
    check_output(cmd1, shell=True)


def dispath_collect_jobs(n):
    global sf_fai
    global sf_bam
    cmd_list=[]
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            cmd="samtools view {0} \"{1}\" | python collect_discordant_low_mapq_reads.py -"\
                .format(sf_bam, fields[0])
            cmd_list.append(cmd)

    pool = ThreadPool(n)
    pool.map(run_cmd, cmd_list)
    pool.close()
    pool.join()


#merge the discordant, clip, 1-0 mapped reads, and dispatch to each gap
#generate final list
def merge_reads_for_gaps(sf_fai, sf_raw_left, sf_raw_right):
    m_scaffold_id={}
    l_scaffold_id=[]
    cnt=0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1
                l_scaffold_id.append(scaffold_id)

    m_left_reads={}
    m_right_reads={}

    s_reads_left="left_reads.list"
    f_left=open(s_reads_left,"w")
    s_reads_right="right_reads.list"
    f_right=open(s_reads_right,"w")

    #get all the discordant read list
    for key in m_scaffold_id:
        sf_discord_left="discordant_reads_list/{0}_cluster_by_discordant_reads_left.list".format(key)
        if os.path.exists(sf_discord_left)==False:
            continue
        f_discord_left=open(sf_discord_left)
        for line in f_discord_left:
            fields=line.split()
            m_left_reads[fields[0]]=1
        f_discord_left.close()

        sf_discord_right="discordant_reads_list/{0}_cluster_by_discordant_reads_right.list".format(key)
        f_discord_right=open(sf_discord_right)
        for line in f_discord_right:
            fields=line.split()
            m_right_reads[fields[0]]=1
        f_discord_right.close()

    for key in m_left_reads:
        f_left.write(key+"\n")
    for key in m_right_reads:
        f_right.write(key+"\n")
    del m_left_reads
    del m_right_reads

    #get all the clip and one-map-mate-unmap reads
    for key in m_scaffold_id:
        sf_clip_10_left="scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(key)
        f_clip_10_left=open(sf_clip_10_left)
        for line in f_clip_10_left:
            fields=line.split()
            if fields[3]!="discordant":
                f_left.write(fields[0]+"\n")
        f_clip_10_left.close()

        sf_clip_10_right="scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(key)
        f_clip_10_right=open(sf_clip_10_right)
        for line in f_clip_10_right:
            fields=line.split()
            if fields[3]!="discordant":
                f_right.write(fields[0]+"\n")
        f_clip_10_right.close()

    f_left.close()
    f_right.close()

#get all the reads by the list
    s_final_left_reads="left_reads.fastq"
    s_final_right_reads="right_reads.fastq"

    cmd_list=[]
    #cmd="/data2/chongchu/ParseFastqById/seqtk {0} {1} {2}".format(sf_raw_left, s_reads_left, s_final_left_reads)
    cmd="/gpfs/ddntempfs/chc12015/tools/ParseFastqById/seqtk {0} {1} {2}".format(sf_raw_left, s_reads_left, s_final_left_reads)
    cmd_list.append(cmd)
    #cmd="/data2/chongchu/ParseFastqById/seqtk  {0} {1} {2}".format(sf_raw_right, s_reads_right ,s_final_right_reads)
    cmd="/gpfs/ddntempfs/chc12015/tools/ParseFastqById/seqtk  {0} {1} {2}".format(sf_raw_right, s_reads_right ,s_final_right_reads)
    cmd_list.append(cmd)

    pool = ThreadPool(2)
    pool.map(run_cmd, cmd_list)
    pool.close()
    pool.join()


#merge the discordant, clip, 1-0 mapped reads, and dispatch to each gap
#generate final list
def merge_dispatch_reads_for_gaps(sf_fai, sf_raw_left, sf_raw_right):
    m_scaffold_id={}
    l_scaffold_id=[]
    cnt=0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1
                l_scaffold_id.append(scaffold_id)

    m_left_reads={}
    m_right_reads={}

    #get all the discordant read list
    for key in m_scaffold_id:
        sf_discord_left="discordant_reads_list/{0}_cluster_by_discordant_reads_left.list".format(key)
        if os.path.exists(sf_discord_left)==False:
            continue
        f_discord_left=open(sf_discord_left)
        for line in f_discord_left:
            fields=line.split()
            if m_left_reads.has_key(fields[0])==False:
                m_left_reads[fields[0]]={}
            m_left_reads[fields[0]][fields[1]]=1
        f_discord_left.close()

        sf_discord_right="discordant_reads_list/{0}_cluster_by_discordant_reads_right.list".format(key)
        f_discord_right=open(sf_discord_right)
        for line in f_discord_right:
            fields=line.split()
            if m_right_reads.has_key(fields[0])==False:
                m_right_reads[fields[0]]={}
            m_right_reads[fields[0]][fields[1]]=1
        f_discord_right.close()

    # for key in m_left_reads:
    #     f_left.write(key+"\n")
    # for key in m_right_reads:
    #     f_right.write(key+"\n")
    # del m_left_reads
    # del m_right_reads

    #get all the clip and one-map-mate-unmap reads
    for key in m_scaffold_id:
        sf_clip_10_left="scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(key)
        if os.path.exists(sf_clip_10_left)==False:
            continue
        f_clip_10_left=open(sf_clip_10_left)
        key_id=m_scaffold_id[key]

        for line in f_clip_10_left:
            fields=line.split()
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_left.write(fields[0]+"\n")
            if m_left_reads.has_key(fields[0])==False:
                m_left_reads[fields[0]]={}
            m_left_reads[fields[0]][sbelong]=1
        f_clip_10_left.close()

        sf_clip_10_right="scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(key)
        f_clip_10_right=open(sf_clip_10_right)
        for line in f_clip_10_right:
            fields=line.split()
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_right.write(fields[0]+"\n")
            if m_right_reads.has_key(fields[0])==False:
                m_right_reads[fields[0]]={}
            m_right_reads[fields[0]][sbelong]=1
        f_clip_10_right.close()

    s_reads_left="left_reads.list"
    f_left=open(s_reads_left,"w")
    s_reads_right="right_reads.list"
    f_right=open(s_reads_right,"w")

    for key_read in m_left_reads:
        for key_gap in m_left_reads[key_read]:
            f_left.write(key_gap+" "+key_read+"\n")

    for key_read in m_right_reads:
        for key_gap in m_right_reads[key_read]:
            f_right.write(key_gap+" "+key_read+"\n")
    f_left.close()
    f_right.close()

    #dispatch reads
    #clear old results
    if os.path.exists("./empty_dir")==False:
        cmd="mkdir ./empty_dir"
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rsync -a --delete empty_dir/ gap_reads/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()


    #deal with left reads
    cnt=0
    read_id=""
    read_head=""
    seq=""
    quality=""
    m_dispatched_reads={}
    with open(sf_raw_left) as fin_left_raw:
        for line in fin_left_raw:
            if cnt%4==0:
                id_fields=line.split()
                read_id=id_fields[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_left_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_1",seq,"+",quality)
                    for key_belong in m_left_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>500:
                            sf_gap_reads="gap_reads/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into file
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()

    ##deal with right reads
    cnt=0
    with open(sf_raw_right) as fin_right_raw:
        for line in fin_right_raw:
            if cnt%4==0:
                id_fields=line.split()
                read_id=id_fields[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_right_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_2",seq,"+",quality)
                    for key_belong in m_right_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>300:
                            sf_gap_reads="gap_reads/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into file
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()


def merge_dispatch_reads_for_gaps_v2(sf_fai, sf_raw_left, sf_raw_right):
    m_scaffold_id={}
    l_scaffold_id=[]
    cnt=0
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            scaffold_id=fields[0]
            m_scaffold_id[scaffold_id]=cnt
            cnt=cnt+1
            l_scaffold_id.append(scaffold_id)

    m_left_reads={}

    #get all the discordant read list
    for key in m_scaffold_id:
        sf_discord_left="discordant_reads_list/{0}_cluster_by_discordant_reads_left.list".format(key)
        if os.path.exists(sf_discord_left)==False:
            continue
        f_discord_left=open(sf_discord_left)
        for line in f_discord_left:
            fields=line.split()
            if m_left_reads.has_key(fields[0])==False:
                m_left_reads[fields[0]]={}
            m_left_reads[fields[0]][fields[1]]=1
        f_discord_left.close()

    #get all the clip and one-map-mate-unmap reads
    for key in m_scaffold_id:
        sf_clip_10_left="scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(key)
        if os.path.exists(sf_clip_10_left)==False:
            continue
        f_clip_10_left=open(sf_clip_10_left)
        key_id=m_scaffold_id[key]

        for line in f_clip_10_left:
            fields=line.split()
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_left.write(fields[0]+"\n")
            if m_left_reads.has_key(fields[0])==False:
                m_left_reads[fields[0]]={}
            m_left_reads[fields[0]][sbelong]=1
        f_clip_10_left.close()

    s_reads_left="left_reads.list"
    f_left=open(s_reads_left,"w")
    for key_read in m_left_reads:
        for key_gap in m_left_reads[key_read]:
            f_left.write(key_gap+" "+key_read+"\n")
    f_left.close()

    #dispatch reads
    #clear old results
    if os.path.exists("./empty_dir")==False:
        cmd="mkdir ./empty_dir"
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rsync -a --delete empty_dir/ gap_reads/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #deal with left reads
    cnt=0
    read_id=""
    read_head=""
    seq=""
    quality=""
    m_dispatched_reads={}
    with open(sf_raw_left) as fin_left_raw:
        for line in fin_left_raw:
            if cnt%4==0:
                id_fields=line.split()
                temp_field=id_fields[0].split("/")
                read_id=temp_field[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_left_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_1",seq,"+",quality)
                    for key_belong in m_left_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>500:
                            sf_gap_reads="gap_reads/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into file
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()
    m_left_reads.clear()

    ###################################################################################################################
    ##deal with right reads
    m_right_reads={}
    for key in m_scaffold_id:
        sf_discord_right="discordant_reads_list/{0}_cluster_by_discordant_reads_right.list".format(key)
        if os.path.exists(sf_discord_right)==False:
            continue
        f_discord_right=open(sf_discord_right)
        for line in f_discord_right:
            fields=line.split()
            if m_right_reads.has_key(fields[0])==False:
                m_right_reads[fields[0]]={}
            m_right_reads[fields[0]][fields[1]]=1
        f_discord_right.close()

    for key in m_scaffold_id:
        sf_clip_10_right="scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(key)
        if os.path.exists(sf_clip_10_right)==False:
            continue
        f_clip_10_right=open(sf_clip_10_right)
        key_id=m_scaffold_id[key]
        for line in f_clip_10_right:
            fields=line.split()
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_right.write(fields[0]+"\n")
            if m_right_reads.has_key(fields[0])==False:
                m_right_reads[fields[0]]={}
            m_right_reads[fields[0]][sbelong]=1
        f_clip_10_right.close()

    s_reads_right="right_reads.list"
    f_right=open(s_reads_right,"w")
    for key_read in m_right_reads:
        for key_gap in m_right_reads[key_read]:
            f_right.write(key_gap+" "+key_read+"\n")
    f_right.close()

    cnt=0
    with open(sf_raw_right) as fin_right_raw:
        for line in fin_right_raw:
            if cnt%4==0:
                id_fields=line.split()
                temp_field=id_fields[0].split("/")
                read_id=temp_field[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_right_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_2",seq,"+",quality)
                    for key_belong in m_right_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>500:
                            sf_gap_reads="gap_reads/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into file
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()
    m_right_reads.clear()


#merge the discordant, clip, 1-0 mapped reads, and dispatch to each gap
#generate final list
def dispatch_reads_for_gaps_to_validate_contigs(sf_fai, sf_raw_left, sf_raw_right):
    m_scaffold_id={}
    l_scaffold_id=[]
    cnt=0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1
                l_scaffold_id.append(scaffold_id)

    m_left_reads={}
    m_right_reads={}

    #get all the clip, one-map-mate-unmap, and only discordant ones
    for key in m_scaffold_id:
        sf_clip_10_left="scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(key)
        if os.path.exists(sf_clip_10_left)==False:
            continue
        f_clip_10_left=open(sf_clip_10_left)
        key_id=m_scaffold_id[key]

        for line in f_clip_10_left:
            fields=line.split()
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_left.write(fields[0]+"\n")
            if m_left_reads.has_key(fields[0])==False:
                m_left_reads[fields[0]]={}
            m_left_reads[fields[0]][sbelong]=1
        f_clip_10_left.close()

        sf_clip_10_right="scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(key)
        f_clip_10_right=open(sf_clip_10_right)
        for line in f_clip_10_right:
            fields=line.split()
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_right.write(fields[0]+"\n")
            if m_right_reads.has_key(fields[0])==False:
                m_right_reads[fields[0]]={}
            m_right_reads[fields[0]][sbelong]=1
        f_clip_10_right.close()

    #clear old results
    if os.path.exists("./empty_dir")==False:
        cmd="mkdir ./empty_dir"
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rsync -a --delete empty_dir/ gap_reads_for_alignment/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #deal with left reads
    cnt=0
    read_id=""
    read_head=""
    seq=""
    quality=""
    m_dispatched_reads={}
    with open(sf_raw_left) as fin_left_raw:
        for line in fin_left_raw:
            if cnt%4==0:
                id_fields=line.split()
                read_id=id_fields[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_left_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_1",seq,"+",quality)
                    for key_belong in m_left_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>500:
                            sf_gap_reads="gap_reads_for_alignment/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into file
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads_for_alignment/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()

    ##deal with right reads
    cnt=0
    with open(sf_raw_right) as fin_right_raw:
        for line in fin_right_raw:
            if cnt%4==0:
                id_fields=line.split()
                read_id=id_fields[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_right_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_2",seq,"+",quality)
                    for key_belong in m_right_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>300:
                            sf_gap_reads="gap_reads_for_alignment/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into files
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads_for_alignment/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()


def dispatch_high_quality_reads_for_gaps(sf_fai, sf_raw_left, sf_raw_right):
    m_scaffold_id={}
    l_scaffold_id=[]
    cnt=0
    with open(sf_fai) as fin_fai:
            for line in fin_fai:
                fields=line.split()
                scaffold_id=fields[0]
                m_scaffold_id[scaffold_id]=cnt
                cnt=cnt+1
                l_scaffold_id.append(scaffold_id)

    m_left_reads={}
    #get all the clip, one-map-mate-unmap, and only discordant ones
    for key in m_scaffold_id:
        sf_clip_10_left="scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(key)
        if os.path.exists(sf_clip_10_left)==False:
            continue
        f_clip_10_left=open(sf_clip_10_left)
        key_id=m_scaffold_id[key]

        for line in f_clip_10_left:
            fields=line.split()
            mapq=int(fields[2])
            if mapq!=60:
                continue
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_left.write(fields[0]+"\n")
            if m_left_reads.has_key(fields[0])==False:
                m_left_reads[fields[0]]={}
            m_left_reads[fields[0]][sbelong]=1
        f_clip_10_left.close()

    #clear old results
    if os.path.exists("./empty_dir")==False:
        cmd="mkdir ./empty_dir"
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rsync -a --delete empty_dir/ gap_reads_high_quality/"
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #deal with left reads
    cnt=0
    read_id=""
    read_head=""
    seq=""
    quality=""
    m_dispatched_reads={}
    with open(sf_raw_left) as fin_left_raw:
        for line in fin_left_raw:
            if cnt%4==0:
                id_fields=line.split()
                temp_field=id_fields[0].split("/")
                read_id=temp_field[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_left_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_1",seq,"+",quality)
                    for key_belong in m_left_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>500:
                            sf_gap_reads="gap_reads_high_quality/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into file
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads_high_quality/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()
    m_left_reads.clear()


    ############################################################################################################
    m_right_reads={}
    for key in m_scaffold_id:
        key_id=m_scaffold_id[key]
        sf_clip_10_right="scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(key)
        if os.path.exists(sf_clip_10_right)==False:
            continue
        f_clip_10_right=open(sf_clip_10_right)
        for line in f_clip_10_right:
            fields=line.split()
            mapq=int(fields[2])
            if mapq!=60:
                continue
            gap_id=fields[1]
            sbelong="{0}_{1}".format(key_id,gap_id)
            #if fields[3]!="discordant":
                #f_right.write(fields[0]+"\n")
            if m_right_reads.has_key(fields[0])==False:
                m_right_reads[fields[0]]={}
            m_right_reads[fields[0]][sbelong]=1
        f_clip_10_right.close()

    ##deal with right reads
    cnt=0
    with open(sf_raw_right) as fin_right_raw:
        for line in fin_right_raw:
            if cnt%4==0:
                id_fields=line.split()
                temp_field=id_fields[0].split("/")
                read_id=temp_field[0][1:].rstrip()
                read_head=line.rstrip()
            elif cnt%4==1:
                seq=line.rstrip()
            elif cnt%4==3:
                quality=line.rstrip()
                if m_right_reads.has_key(read_id)==True:
                    record="{0}\n{1}\n{2}\n{3}\n".format("@"+read_id+"_2",seq,"+",quality)
                    for key_belong in m_right_reads[read_id]:
                        if m_dispatched_reads.has_key(key_belong)==False:
                            m_dispatched_reads[key_belong]=[]
                        m_dispatched_reads[key_belong].append(record)
                        if len(m_dispatched_reads[key_belong])>300:
                            sf_gap_reads="gap_reads_high_quality/{0}.fastq".format(key_belong)
                            with open(sf_gap_reads,"a+") as fin_gap_reads:
                                for info in m_dispatched_reads[key_belong]:
                                    fin_gap_reads.write(info)
                            del m_dispatched_reads[key_belong][:]
            cnt=cnt+1

        #write the rest into files
        for sbl in m_dispatched_reads:
            sf_gap_reads="gap_reads_high_quality/{0}.fastq".format(sbl)
            with open(sf_gap_reads,"a+") as fin_gap_reads:
                for info in m_dispatched_reads[sbl]:
                    fin_gap_reads.write(info)
            del m_dispatched_reads[sbl][:]
    m_dispatched_reads.clear()
    m_right_reads.clear()


# #each gap has a buffer and corresponding file
# #dispatch the reads for each gap
# def dispatch_reads_for_gaps(sf_fai, sf_all_left_reads, sf_all_right_reads):
#     with open(sf_all_left_reads) as fin_left_reads:
if __name__ == "__main__":
    sf_fai=sys.argv[1]
    sf_bam=sys.argv[2]
    sfout_discord_pos=sys.argv[3]

    collect_discordant_regions_v2(sf_fai,sfout_discord_pos)
    dispath_collect_jobs(int(sys.argv[4]))

    merge_dispatch_reads_for_gaps_v2(sf_fai, sys.argv[5], sys.argv[6])
    #dispatch_reads_for_gaps_to_validate_contigs(sf_fai, sys.argv[5], sys.argv[6])
    dispatch_high_quality_reads_for_gaps(sf_fai, sys.argv[5], sys.argv[6])
