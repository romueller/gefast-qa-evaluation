#!/usr/bin/python

# Convenience script / method for converting a FASTQ file to a FASTA file.

from Bio import SeqIO


# Generate a FASTA file (two-line format) from the given FASTQ file by dropping the quality scores.
def convert(fastq_file, fasta_file):
    with open(fastq_file, 'r') as in_file, open(fasta_file, 'w') as out_file:
        SeqIO.write((record for record in SeqIO.parse(in_file, 'fastq')), out_file, 'fasta-2line')


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('fastq_file', help = 'input FASTQ file')
    argparser.add_argument('fasta_file', help = 'output FASTA file')
    args = argparser.parse_args()

    convert(args.fastq_file, args.fasta_file)