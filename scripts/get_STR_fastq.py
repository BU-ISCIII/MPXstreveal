#!/usr/bin/env python
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: get_STR_fastq.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch

from __future__ import division
import os,re
from pathlib import Path
from multiprocessing import Pool

try:
    import configparser
except:
    from six.moves import configparser


def conf_parse(conf_file):
    cf = configparser.ConfigParser()
    cf.read(conf_file)
    conf_dict = {}
    for section in cf.sections():
        for (key,value) in cf.items(section):
            conf_dict[key] = cf.get(section, key)
    return conf_dict



def get_STR_fq_from_bam(args_list):
    '''
    (1)Get STR bam forward stand and reverse stand bam from  bam file
    (2)then transfor to fastq file using bamToFastq
    (3)get the inverse complementary sequence from reverse fastq file
    (4)combined forward and reverse strand fastq files

    '''
    sample,pos,marker_name,bam_file,result_dir,stand,mysamtools, mybamToFastq,myseqtk,myusearch,type,assemble_pairs= args_list
    COMMAND_index="{0} index -b {1}"
    if not os.path.exists(bam_file+".bai"):
        os.system(COMMAND_index.format(mysamtools,bam_file))

    STR_bam = os.path.join(result_dir,marker_name+"_reads_"+sample+".bam")
    #STR_bam_sort = os.path.join(result_dir,marker_name+"_reads_"+sample+"_sortByname.bam")
    STR_bam_fw = os.path.join(result_dir,marker_name+"_reads_"+sample+"_fw.bam")
    STR_bam_rev = os.path.join(result_dir,marker_name+"_reads_"+sample+"_rev.bam")
    STR_bam_sort_fw = os.path.join(result_dir,marker_name+"_reads_"+sample+"_sortByname_fw.bam")
    STR_bam_sort_rev = os.path.join(result_dir,marker_name+"_reads_"+sample+"_sortByname_rev.bam")
    #fq_R1 = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R1.fastq")
    #fq_R2 = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R2.fastq")
    fq_rev_tmp = os.path.join(result_dir,marker_name+"_reads_"+sample+"_rev_tmp.fastq")
    fq_rev = os.path.join(result_dir,marker_name+"_reads_"+sample+"_rev.fastq")
    fq_fw = os.path.join(result_dir,marker_name+"_reads_"+sample+"_fw.fastq")
    fq_R1_all = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R1_all.fastq")
    fq_R2_all = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R2_all.fastq")
    merge_fq = os.path.join(result_dir,marker_name+"_reads_"+sample+"_merge.fastq")
    merge_fq_forward =os.path.join(result_dir,marker_name+"_reads_"+sample+"_merge_forward.fastq")
    COMMAND_samtools = "{0} view -b1 {1} {2}>{3}"
    COMMAND_sort = "{0} sort -n {1}>{2}"
    COMMAND_samtoolsReverse = "{0} view -b1 -f 0x10 {1}>{2}"
    COMMAND_samtoolsForward = "{0} view -b1 -F 0x10 {1}>{2}"
    #COMMAND_bamToFastq = "{0} -i {1} -fq {2} -fq2 {3}>/dev/null 2>&1"
    if type == "single":
        COMMAND_bamToFastq = "{0} fastq {1}>{2}"
    else:
        COMMAND_bamToFastq = "{0} fastq -0 /dev/null {1}>{2}"
    COMMAND_reverse = "{0} seq -r {1}>{2}"
    COMMAND_merge = "{0} -fastq_mergepairs {1} -reverse {2} -fastqout {3}"

    os.system(COMMAND_samtools.format(mysamtools,bam_file,pos,STR_bam))
    try:
        os.system(COMMAND_samtoolsReverse.format(mysamtools,STR_bam,STR_bam_rev))
        os.system(COMMAND_sort.format(mysamtools,STR_bam_rev,STR_bam_sort_rev))
        os.system(COMMAND_bamToFastq.format(mysamtools,STR_bam_sort_rev,fq_rev_tmp))

        os.system(COMMAND_reverse.format(myseqtk,fq_rev_tmp,fq_rev))

        os.system(COMMAND_samtoolsForward.format(mysamtools,STR_bam,STR_bam_fw))
        os.system(COMMAND_sort.format(mysamtools,STR_bam_fw,STR_bam_sort_fw))
        os.system(COMMAND_bamToFastq.format(mysamtools,STR_bam_sort_fw,fq_fw))
        #os.system("grep -s -A 3 '/1' {0} | grep -v '^--'>{1}".format(fq_fw,fq_R1_all))
        #os.system("grep -s -A 3 '/2' {0} | grep -v '^--'>{1}".format(fq_rev,fq_R2_all))
    except Exception as e:
        print(e)
    os.system("cat {0} {1}>{2}".format(fq_fw,fq_rev,merge_fq))
#    os.system(COMMAND_merge.format(myusearch,fq_R1_all,fq_R2_all,merge_fq)) #TODO if assembly pairs true


#        if stand == "+":
#            fq_R2_reverse = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R2_trans.fastq")
#            os.system(COMMAND_reverse.format(myseqtk,fq_R2,fq_R2_reverse))
#            if assemble_pairs ==False:
#                os.system("cat {0} {1}>{2}".format(fq_R1,fq_R2_reverse,merge_fq))
#            else:
#                os.system(COMMAND_merge.format(myusearch,fq_R1,fq_R2,merge_fq))
#        else:
#            fq_R1_reverse = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R1_trans.fastq")
#            os.system(COMMAND_reverse.format(myseqtk,fq_R1,fq_R1_reverse))
#            if assemble_pairs ==False:
#                os.system("cat {0} {1}>{2}".format(fq_R1_reverse,fq_R2,merge_fq))
#            else:
#                os.system(COMMAND_merge.format(myusearch,fq_R1,fq_R2,merge_fq_forward))
#                os.system(COMMAND_reverse.format(myseqtk,merge_fq_forward,merge_fq))
    # else:
    #     print("Entra en single")
    #     COMMAND_bamToFastq = "{0} -i {1} -fq {2}"
    #     try:
    #         os.system(COMMAND_bamToFastq.format(mybamToFastq,STR_bam_sort,fq_R1))
    #         if stand == "+":
    #             os.system("mv {0} {1}".format(fq_R1,merge_fq))
    #         else:
    #             os.system(COMMAND_reverse.format(myseqtk,fq_R1,fq_R1_reverse))
    #             os.system("mv {0} {1}".format(fq_R1_reverse,merge_fq))
    #     except Exception as e:
    #         print(e)
    print("{}: finished getting STR fq from bam!".format(marker_name))


def check_bam_file(bam_file):
    my_file = Path(bam_file)
    try:
        bam = my_file.resolve()
    except FileNotFoundError:
        print("Bam file doesn't exist!!")


def main(sample,bam_file,working_path,ref_bed,type,assemble_pairs,num_processors):
    ## get path of linux tools
    conf_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    config = conf_parse(os.path.join(conf_path,"conf.py"))
    mysamtools=eval(config["mysamtools"])
    mybamToFastq=eval(config["mybamtofastq"])
    myseqtk=eval(config["myseqtk"])
    myusearch=eval(config["myusearch"])

    ## check Bam file
    check_bam_file(bam_file)
    ## create STR fq dir
    STR_fastq_dir=os.path.join(working_path,"STRfq")
    if not os.path.exists(STR_fastq_dir):
        os.makedirs(STR_fastq_dir)

    info_list = []
    bed = open(ref_bed,"r")
    next(bed)
    N = 0
    for line in bed:
        line = line.strip()
        mylist = line.split("\t")
        pos = mylist[0]+":"+mylist[1]+"-"+mylist[2]
        marker_name,stand = mylist[5],mylist[8]
        info_list.append([sample,pos,marker_name,bam_file,STR_fastq_dir,stand,mysamtools, mybamToFastq,myseqtk,myusearch,type,assemble_pairs])
        N+=1

    pool = Pool(num_processors)
    pool.imap(get_STR_fq_from_bam, info_list)
    pool.close()
    pool.join()
    print('Finished getting STR fq from bam for total {} markers!'.format(N))


if __name__=='__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],sys.argv[7]))
