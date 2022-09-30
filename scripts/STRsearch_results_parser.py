#!/usr/bin/env python

####IMPORTS#####
import re
import sys
import argparse

####FUNCTIONS######
def parse_args(args=None):
    Description = "Filter STRsearch output"
    Epilog = """Example usage: python3 str_parse.py <merged_fastq> <flanking_regions> <number_flanking> <mismatches> <out_fasta> <out_long>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("merged_fastq", help="Merged fastq file from STRsearch")
    parser.add_argument("flanking_regions", help="File with reads containing flanking regions from custom STRsearch")
    parser.add_argument("number_flanking", help="Match both flanking regions or only one of them")
    parser.add_argument("mismatches", help="Number of mismatches allowed in the flanking regions")
    parser.add_argument("out_fasta", help="Output multifasta file")
    parser.add_argument("out_long", help="Output long table")
    return parser.parse_args(args)

def create_dict_multifa (merged_fastq, flanking_regions, number_flanking, mismatches, out_fasta, str_seq_dict, str_mark):
    with open(flanking_regions, 'r') as fr_in:
        next(fr_in)
        fout = open(out_fasta, "w")
        for line in fr_in:
            line = re.split("\t", line)
            sample_name = line[0]
            str_pattern = line[1]
            allele_number = line[2]
            distance_5 = line[3]
            distance_3 = line[4]
            read_name = line[5]
            flank_left_start = line[6]
            flank_left_end = line[7]
            flank_left_mistmatches = line[8]
            flank_right_start = line[9]
            flank_right_end = line[10]
            flank_right_mistmatches = line[11]
            if number_flanking == "both" and int(flank_left_mistmatches) <= int(mismatches) and int(flank_right_mistmatches) <= int(mismatches):
                with open(merged_fastq, 'r') as merge_in:
                    read = []
                    for lines in merge_in:
                        read.append(lines.rstrip())
                        if len(read) == 4:
                            header=read[0]
                            sequence=read[1]
                            if read_name == header :
                                if int(flank_left_end) < int(flank_right_start):
                                    out_header='>'+header+'\n'
                                    fout.write(out_header)
                                    out_sequence=sequence[int(flank_left_start):int(flank_right_end)]+'\n'
                                    fout.write(out_sequence)
                                    if str_pattern not in str_seq_dict:
                                        new_key = { str_pattern:
                                                   { sequence[int(flank_left_start):int(flank_right_end)] : {
                                                           'sample_name': sample_name,
                                                           'str_mark' : str_mark,
                                                           'num_reps': allele_number,
                                                           'supporting_reads': 1,
                                                           'allele_frequency':''
                                                        }
                                                      }
                                                  }
                                        str_seq_dict.update(new_key)
                                    else:
                                        if sequence[int(flank_left_start):int(flank_right_end)] in str_seq_dict[str_pattern]:

                                            old_support=str_seq_dict[str_pattern][sequence[int(flank_left_start):int(flank_right_end)]]['supporting_reads']
                                            new_support=old_support+1
                                            str_seq_dict[str_pattern][sequence[int(flank_left_start):int(flank_right_end)]]['supporting_reads']=new_support


                                        else :
                                            new_sequence={ sequence[int(flank_left_start):int(flank_right_end)] : {
                                                           'sample_name': sample_name,
                                                           'str_mark' : str_mark,
                                                           'num_reps': allele_number,
                                                           'supporting_reads': 1,
                                                           'allele_frequency':''
                                                        }
                                                      }
                                            str_seq_dict[str_pattern].update(new_sequence)
                                else:
                                    print("Reverse read:")
                                    print(read_name+'\t'+flank_left_start+'\t'+flank_left_end+'\t'+flank_right_start+'\t'+flank_right_end)
                                    print(header)
                                    print(sequence)
                                    print(sequence[int(flank_right_start):int(flank_left_end)])
                            read = []

    #         else: #TODO
    #             print("Solo un flanking")
    return(str_seq_dict)

def calculate_total_suppport(str_seq_dict,all_support):
    for str_id in str_seq_dict:
        for sequence in str_seq_dict[str_id]:
            support_value = str_seq_dict[str_id][sequence]['supporting_reads']
            all_support=all_support+support_value
    return(all_support)

def calculate_frequency(str_seq_dict,all_support):
    for str_id in str_seq_dict:
        for sequence in str_seq_dict[str_id]:
            support_value = str_seq_dict[str_id][sequence]['supporting_reads']
            allele_freq = support_value/all_support
            str_seq_dict[str_id][sequence]['allele_frequency'] = allele_freq
    return(str_seq_dict)

def create_long_table(str_seq_dict,out_long):
    lout = open(out_long, "w")
    header="Sample_name\tSTR_mark\tSTR_structure\tNumer_reps\tSupporting_reads\tAlleleFrequency\tSequence\n"
    lout.write(header)
    for str_id in str_seq_dict:
        for sequence in str_seq_dict[str_id]:
            oline=str_seq_dict[str_id][sequence]['sample_name']+'\t'+str_seq_dict[str_id][sequence]['str_mark']+'\t'+str_id+'\t'+str(str_seq_dict[str_id][sequence]['num_reps'])+'\t'+str(str_seq_dict[str_id][sequence]['supporting_reads'])+'\t'+str(str_seq_dict[str_id][sequence]['allele_frequency'])+'\t'+sequence+'\n'
            lout.write(oline)

def main(args=None):
    args = parse_args(args)
    str_mark=re.split('/',re.split("_reads", args.merged_fastq)[0])[1]
    str_seq_dict={}
    initial_dict=create_dict_multifa(args.merged_fastq, args.flanking_regions, args.number_flanking, args.mismatches, args.out_fasta, str_seq_dict, str_mark)
    all_support=0
    all_support=calculate_total_suppport(initial_dict,all_support)
    final_dict=calculate_frequency(initial_dict,all_support)
    create_long_table(final_dict,args.out_long)


if __name__ == "__main__":
    sys.exit(main())
