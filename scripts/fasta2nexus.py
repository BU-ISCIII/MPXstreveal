from Bio import SeqIO
import argparse

def parse_args(args=None):
    Description = "Convert fasta to nexus"
    Epilog = """Example usage: python fasta2nexus.py <in_file> <file_out>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("in_file", help="Path to input fasta file to convert.")
    parser.add_argument("file_out", help="Path to print output High Quality genome fasta file.")
    return parser.parse_args(args)

args = parse_args()

fasta_file = args.in_file
nexus_file = args.file_out

count = SeqIO.convert(fasta_file, "fasta", nexus_file, "nexus", "DNA")
print("Converted %i records" % count)
