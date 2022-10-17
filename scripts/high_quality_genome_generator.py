################################################################################################
### Create hight quality genome including curated STR repeats in genome consensus with ivar  ###
################################################################################################

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import sys
import argparse

def parse_args(args=None):
    Description = "Create hight quality genome including curated STR repeats in genome consensus with ivar"
    Epilog = """Example usage: python high_quality_genome_generator.py <in_file> <in_tsv> <file_out>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("in_file", help="Path to input fasta file for STR introduction.")
    parser.add_argument("in_tsv", help="Path to TSV file with STR anotation information.")
    parser.add_argument("file_out", help="Path to print output High Quality genome fasta file.")
    return parser.parse_args(args)

def expand_str(pattern):
  pattern = pattern.split(" ")
  print(pattern)
  expanded_seq = ""
  for pt in pattern:
      x = re.search(r"\[(.*)\](\d+)",pt)
      if x:
        pat = x.group(1)
        repeat = x.group(2)
        expanded_seq += pat*int(repeat)
      else:
        expanded_seq += pt
  return expanded_seq

args = parse_args()

#### csv with str info
str_info_path = args.in_tsv
#### consensus genome
genome_path = args.in_file

out_file = args.file_out

## Read input files
str_info_file = open(str_info_path,"r")
str_info = {}
next(str_info_file)
for line in str_info_file:
    line = line.rstrip()
    line = line.split("\t")
    str_region = {}
    str_region["flanking_left"] = line[9]
    str_region["flanking_right"] = line[10]
    str_region["pattern"] = line[11]
    str_info[line[5]] = str_region
record = list(SeqIO.parse(genome_path, "fasta"))[0]

output_seq = record
len_seq = len(record.seq)
for k in str_info.keys():
    print(k)
    begin_seq_index = output_seq.seq.find(str_info[k]["flanking_left"])
    begin_seq = output_seq[0:begin_seq_index + len(str_info[k]["flanking_left"])]
    str_seq = expand_str(str_info[k]["pattern"])
    end_seq_index = output_seq.seq.find(str_info[k]["flanking_right"])
    print(end_seq_index)
    end_seq = output_seq[end_seq_index:len_seq]
    output_seq = begin_seq + str_seq + end_seq
    len_seq = len(output_seq)

with open(out_file, "w") as output_handle:
    SeqIO.write(output_seq, output_handle, "fasta")
