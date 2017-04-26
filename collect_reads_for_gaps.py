import sys

class GapReadsCollector:
    def __init__(self, insert_size, derivation, dist_clip):
        self.dist1=insert_size-3*derivation
        self.dist2=insert_size+3*derivation
        self.dist_clip=dist_clip #by default 250

    #0: unclipped
    #1: left-clipped
    #2: right-clipped
    #3: both-clipped
    def is_clipped(self, cigar):
        cnt=0
        l=len(cigar)
        if cigar[l-1]=="S" or cigar[l-1]=="H":#right clipped
            cnt=2

        for i in range(l):
            if cigar[i]>="0" and cigar[i]<="9":
                continue
            else:
                if cigar[i]=="S" or cigar[i]=="H": ##left clipped
                    cnt=cnt+1
                break
        return cnt


    #need to consider the situation that:
    #nearby gaps may share some reads
    def get_focal_region_of_scaffold_v2(self, sf_gap_pos, s_scaffod):
        focal_region={}
        gap_pos={}
        cnt=1
        with open(sf_gap_pos) as fin_gap_pos:
            for line in fin_gap_pos:
                fields=line.split()
                ref=fields[3]
                start=int(fields[0])
                end=int(fields[1])

                if ref==s_scaffod:
                    gap_pos[cnt]=[]
                    gap_pos[cnt].append(start)
                    gap_pos[cnt].append(end)

                    for i in range(self.dist2):
                        if start-i>=0:
                            if focal_region.has_key(start-i)==False:
                                    focal_region[start-i]={}

                            if i <= self.dist_clip:
                                focal_region[start-i][str(cnt)+"_0c"]=1
                            else:
                                focal_region[start-i][str(cnt)+"_0d"]=1

                        if focal_region.has_key(end+i)==False:
                            focal_region[end+i]={}
                        if i <= self.dist_clip:
                            focal_region[end+i][str(cnt)+"_1c"]=1
                        else:
                            focal_region[end+i][str(cnt)+"_1d"]=1
                    cnt=cnt+1

        return focal_region, gap_pos


    def parse_reads_fall_in_gaps_one_scaffold(self, sf_gap_pos, anchor_mapq, working_folder):
        pre_ref=""
        f_cluster_left=open(working_folder+"cluster_by_gap_reads_left.list","w")
        f_cluster_right=open(working_folder+"cluster_by_gap_reads_right.list","w")

        focal_region={}
        gap_pos={}

        for sam_record in sys.stdin:
            sam_fields=sam_record.split()
            #first get first or second in pair
            bfirst=True
            flag=int(sam_fields[1])
            if (flag & 0x40)==0:
                bfirst=False

            cigar=sam_fields[5]
            #read_seq=sam_fields[9]
            qname=sam_fields[0]
            map_pos=int(sam_fields[3])
            map_quality=sam_fields[4]
            ref=sam_fields[2]
            mate_ref=sam_fields[6]
            mate_pos=int(sam_fields[7])

            if pre_ref!=ref:
                #focal_region, gap_pos=get_focal_region_of_scaffold(sf_gap_pos,ref) alignments
                focal_region, gap_pos=self.get_focal_region_of_scaffold_v2(sf_gap_pos,ref)

                f_cluster_left.close()
                f_cluster_right.close()
                sf_cluster_left=working_folder+"scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(ref)
                sf_cluster_right=working_folder+"scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(ref)
                f_cluster_left=open(sf_cluster_left,"w")
                f_cluster_right=open(sf_cluster_right,"w")

            if focal_region.has_key(map_pos)==True:
                for pos_info in focal_region[map_pos]:
                    s_content=str(pos_info)
                    content_fields=s_content.split("_")
                    scnt=content_fields[0]
                    sflag=content_fields[1] ##0c 0d 1c 1d

                    gap_start=int(gap_pos[int(scnt)][0])
                    gap_end=int(gap_pos[int(scnt)][1])
                    gap_len=gap_end-gap_start

                    clip_flag=self.is_clipped(cigar)

                    #left-edge-of-gap and clipped at right  or
                    #right-edge-of-gap and clipped at left
                    if (sflag=="0c" and clip_flag>=2) or (sflag=="1c" and (clip_flag==1 or clip_flag==3)):##clip region
                        if bfirst==True:
                            f_cluster_left.write(qname+" "+scnt+" "+map_quality+" clip\n")
                        else:
                            f_cluster_right.write(qname+" "+scnt+" "+map_quality+" clip\n")

                    ##read and its mate both are mapped and the anchor maqp is satisfied
                    if (flag & 0x4)==0 and (flag & 0x8)==0 and int(map_quality)>=anchor_mapq:
                        if mate_ref != "=":
                            if bfirst==True:#mate is second
                                f_cluster_right.write(qname+" "+scnt + " "+map_quality+
                                    " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos) + " "+str(gap_len) +"\n")
                            else:
                                f_cluster_left.write(qname+" "+scnt + " "+map_quality+
                                    " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos)+ " "+str(gap_len)  +"\n")
        ##note when get region reads, first check the mapq first, if > 30, then discard
        ##as if both reads mapq>30, then both are included
                        else: #same ref, but distance pass the threshold
                            # f1r2=True
                            # if (bfirst==True and map_pos>mate_pos) or (bfirst==False and map_pos<mate_pos):
                            #     f1r2=False
                            temp_insert=int(sam_fields[8])
                            if bfirst==True:
                                #if (abs(temp_insert)-gap_len)>=dist2:
                                if (abs(temp_insert))>=self.dist2:
                                    f_cluster_right.write(qname+" "+scnt + " "+map_quality+
                                         " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos) + " "+str(gap_len) +"\n")
                            else:
                                if (abs(temp_insert))>=self.dist2:
                                #if (abs(temp_insert)-gap_len)>=dist2:
                                    f_cluster_left.write(qname+" "+scnt + " "+map_quality+
                                         " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos)+ " "+str(gap_len) +"\n")

                    ###one map, mate unmapped reads
                    elif (flag & 0x4)==0 and (flag & 0x8)!=0:
                        if bfirst==True:#mate is second
                            f_cluster_right.write(qname+" "+scnt+
                                                  " "+map_quality+" unmap\n") ##note: we are going to collect the mate read, but the mapq is the current read
                        else:
                            f_cluster_left.write(qname+" "+scnt+
                                                 " "+map_quality+" unmap\n")
            pre_ref=ref

        f_cluster_left.close()
        f_cluster_right.close()


    def parse_reads_fall_in_gaps_one_scaffold_short_is(self, sf_gap_pos, anchor_mapq, working_folder): ##for short insert size data
        pre_ref=""
        f_cluster_left=open(working_folder+"cluster_by_gap_reads_left.list","w")
        f_cluster_right=open(working_folder+"cluster_by_gap_reads_right.list","w")

        focal_region={}
        gap_pos={}

        for sam_record in sys.stdin:
            sam_fields=sam_record.split()
            #first get first or second in pair
            bfirst=True
            flag=int(sam_fields[1])
            if (flag & 0x40)==0:
                bfirst=False

            cigar=sam_fields[5]
            #read_seq=sam_fields[9]
            qname=sam_fields[0]
            map_pos=int(sam_fields[3])
            map_quality=sam_fields[4]
            ref=sam_fields[2]
            mate_ref=sam_fields[6]
            mate_pos=int(sam_fields[7])

            if pre_ref!=ref:
                #focal_region, gap_pos=get_focal_region_of_scaffold(sf_gap_pos,ref)
                focal_region, gap_pos=self.get_focal_region_of_scaffold_v2(sf_gap_pos,ref)

                f_cluster_left.close()
                f_cluster_right.close()
                sf_cluster_left=working_folder+"scaffold_reads_list_all/{0}_cluster_by_gap_reads_left.list".format(ref)
                sf_cluster_right=working_folder+"scaffold_reads_list_all/{0}_cluster_by_gap_reads_right.list".format(ref)
                f_cluster_left=open(sf_cluster_left,"w")
                f_cluster_right=open(sf_cluster_right,"w")


            if focal_region.has_key(map_pos)==True:
                for pos_info in focal_region[map_pos]:
                    s_content=str(pos_info)
                    content_fields=s_content.split("_")
                    scnt=content_fields[0]
                    sflag=content_fields[1] ##0c 0d 1c 1d


                    gap_start=int(gap_pos[int(scnt)][0])
                    gap_end=int(gap_pos[int(scnt)][1])
                    gap_len=gap_end-gap_start

                    clip_flag=self.is_clipped(cigar)

                    #left-edge-of-gap and clipped at right  or
                    #right-edge-of-gap and clipped at left
                    if (sflag=="0c" and clip_flag>=2) or (sflag=="1c" and (clip_flag==1 or clip_flag==3)):##clip region
                        if bfirst==True:
                            f_cluster_left.write(qname+" "+scnt+" "+map_quality+" clip\n")
                        else:
                            f_cluster_right.write(qname+" "+scnt+" "+map_quality+" clip\n")

                    ##read and its mate both are mapped and the anchor maqp is satisfied
                    if (flag & 0x4)==0 and (flag & 0x8)==0 and int(map_quality)>=anchor_mapq:
                        if mate_ref != "=":
                            if bfirst==True:#mate is second
                                f_cluster_right.write(qname+" "+scnt + " "+map_quality+
                                    " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos) + " "+str(gap_len) +"\n")
                            else:
                                f_cluster_left.write(qname+" "+scnt + " "+map_quality+
                                    " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos)+ " "+str(gap_len)  +"\n")
        ##note when get region reads, first check the mapq first, if > 30, then discard
        ##as if both reads mapq>30, then both are included
                        else: #same ref, but distance pass the threshold
                            # f1r2=True
                            # if (bfirst==True and map_pos>mate_pos) or (bfirst==False and map_pos<mate_pos):
                            #     f1r2=False
                            temp_insert=int(sam_fields[8])
                            if bfirst==True:
                                #if (abs(temp_insert)-gap_len)>=dist2:
                                if (abs(temp_insert))>=self.dist2 or (abs(temp_insert))<=self.dist1:
                                    f_cluster_right.write(qname+" "+scnt + " "+map_quality+
                                         " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos) + " "+str(gap_len) +"\n")
                            else:
                                if (abs(temp_insert))>=self.dist2 or (abs(temp_insert))<=self.dist1:
                                #if (abs(temp_insert)-gap_len)>=dist2:
                                    f_cluster_left.write(qname+" "+scnt + " "+map_quality+
                                         " discordant "+str(map_pos)+" "+ mate_ref+" "+str(mate_pos)+ " "+str(gap_len) +"\n")

                    ###one map, mate unmapped reads
                    elif (flag & 0x4)==0 and (flag & 0x8)!=0:
                        if bfirst==True:#mate is second
                            f_cluster_right.write(qname+" "+scnt+
                                                  " "+map_quality+" unmap\n") ##note: we are going to collect the mate read, but the mapq is the current read
                        else:
                            f_cluster_left.write(qname+" "+scnt+
                                                 " "+map_quality+" unmap\n")
            pre_ref=ref

        f_cluster_left.close()
        f_cluster_right.close()


if __name__ == "__main__":
    sf_gap_pos=sys.argv[1]
    anchor_mapq=int(sys.argv[2])
    working_folder=sys.argv[3]
    insert_size=int(sys.argv[4])
    derivation=int(sys.argv[5])
    dist_clip=int(sys.argv[6])

    grc=GapReadsCollector(insert_size, derivation, dist_clip)
    if insert_size>=750:
        grc.parse_reads_fall_in_gaps_one_scaffold(sf_gap_pos, anchor_mapq, working_folder)
    else:
        grc.parse_reads_fall_in_gaps_one_scaffold_short_is(self, sf_gap_pos, anchor_mapq, working_folder)
