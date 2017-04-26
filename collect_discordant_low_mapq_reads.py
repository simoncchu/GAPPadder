import sys
import os

def get_focal_region(sf_scaffod):
    m_pos_gaps={}
    focal_region={}

    with open(sf_scaffod) as fin_discordant_pos: ###when no this file happen???
        pre_pos=-1
        for line in fin_discordant_pos:
            fields=line.split()

            map_pos=int(fields[1])

            if map_pos!=pre_pos:
                m_pos_gaps[map_pos]=[]
            belong_to=fields[2]+"_"+fields[3]

            m_pos_gaps[map_pos].append(belong_to)

            for i in range(200):##
                if map_pos-i>=0:
                    focal_region[map_pos-i]=map_pos
            for i in range(300):
                focal_region[map_pos+i]=map_pos

            pre_pos=map_pos
    return focal_region, m_pos_gaps


def parse_discordant_reads_one_scaffold(working_folder):
    pre_ref=""
    f_cluster_left=open(working_folder+"cluster_by_discordant_reads_left.list","w")
    f_cluster_right=open(working_folder+"cluster_by_discordant_reads_right.list","w")

    focal_region={}
    m_pos_gaps={}

    for sam_record in sys.stdin:
        sam_fields=sam_record.split()
        #first get first or second in pair
        bfirst=True
        flag=int(sam_fields[1])
        if (flag & 0x40)==0:
            bfirst=False

        qname=sam_fields[0]
        map_pos=int(sam_fields[3])
        map_quality=int(sam_fields[4])
        ref=sam_fields[2]

        if map_quality>0:
            continue

        if pre_ref!=ref:
            sf_scaffold=working_folder+"discordant_temp/"+ref+".list"
            if os.path.exists(sf_scaffold)==False:
                continue
            focal_region, m_pos_gaps=get_focal_region(sf_scaffold)
            f_cluster_left.close()
            f_cluster_right.close()
            sf_cluster_left=working_folder+"discordant_reads_list/{0}_cluster_by_discordant_reads_left.list".format(ref)
            sf_cluster_right=working_folder+"discordant_reads_list/{0}_cluster_by_discordant_reads_right.list".format(ref)
            f_cluster_left=open(sf_cluster_left,"w")
            f_cluster_right=open(sf_cluster_right,"w")

        if focal_region.has_key(map_pos)==False:
            continue

        source_map_pos=int(focal_region[map_pos])
        if m_pos_gaps.has_key(source_map_pos)==False:
            continue

        if focal_region.has_key(map_pos)==True:
            for sbelong_to in m_pos_gaps[source_map_pos]:
                if bfirst==True:
                    f_cluster_left.write(qname+" "+sbelong_to+" "+str(map_quality)+"\n")
                else:
                    f_cluster_right.write(qname+" "+sbelong_to+" "+str(map_quality)+"\n")

        pre_ref=ref

    f_cluster_left.close()
    f_cluster_right.close()


if __name__ == "__main__":
    working_folder=sys.argv[1]
    parse_discordant_reads_one_scaffold(working_folder)
