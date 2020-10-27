#!/usr/bin/python

# Different methods for the assignment of taxonomic information to sequences.
#
# Expanded and adapted from:
# https://github.com/torognes/vsearch-eval
#
# File: https://github.com/torognes/vsearch-eval/blob/master/cluster/scripts/tax.sh

import subprocess as sp

from Bio import SeqIO

from . import reformat_silva_taxonomy as rst
from .. import util

DEF_THREADS = 8


# Determine the ground truth by matching against the reference sequences.
def match_against_references(sequences, references, threshold, output, vsearch, threads = DEF_THREADS):
    sp.call('%s --usearch_global %s --db %s --id %f --blast6out %s --strand plus --threads %i'
            % (vsearch, sequences, references, threshold, output, threads), shell = True)


# Replace reference identifiers (obtained from FASTA headers) in a taxonomic assignment
# with the taxonomy corresponding to the identifier.
def relabel(taxonomic_assignments, taxonomy, output):
    # parse taxonomy file
    tax_of_rep = dict()
    with open(taxonomy, 'r') as in_file:
        for line in in_file:
            id, taxon = line.rstrip('\n').split('\t')
            taxon = ''.join(taxon.split(' '))
            tax_of_rep[id] = taxon

    # parse taxonomic assignments
    assignments = []
    with open(taxonomic_assignments, 'r') as in_file:
        for line in in_file:
            assignments.append(line.rstrip('\n').split('\t'))

    # look up the taxonomy and replace the identifier in the assignment with the taxonomy
    for elem in assignments:
        elem[1] = tax_of_rep[elem[1]]

    # print relabelled assignments
    with open(output, 'w') as out_file:
        for elem in assignments:
            out_file.write('%s\n' % '\t'.join(elem))


# Produce a reference FASTA file in which the entries are labelled with taxonomic information
# by reading a taxonomy file and a reference FASTA file and keeping only those FASTA entries for which
# the taxonomy file has clean, and sufficiently "complete" taxonomic information
def reduce_and_label(rep_set, taxonomy, reference_file, new_rep_set = None, new_taxonomy = None, level = rst.DEF_LEVEL):
    # read taxonomy information and keep only those entries with clean, "complete" taxonomic information
    tax_info = dict()
    new_taxonomy_lines = []
    with open(taxonomy, 'r') as in_file:
        for line in in_file:
            seq_id, tax = rst.remove_levels(line.strip().split('\t'), level).split('\t')
            if tax.find('Ambiguous_taxa') == -1 and tax.find('D_%i__uncultured' % level) == -1 and tax.find('D_%i__' % level) != -1:
                tax_info[seq_id] = tax
                new_taxonomy_lines.append('%s\t%s\n' % (seq_id, tax))

    # read FASTA file of representatives
    # keep and relabel only those sequences for which we have kept taxonomic information in the previous step
    with open(rep_set, 'r') as in_file, open(reference_file, 'w') as out_file:
        for record in SeqIO.parse(in_file, 'fasta'):
            if record.id in tax_info:
                out_file.write('>%s\n%s\n' % (tax_info[record.id], record.seq))

    # optionally, create reduced representative and taxonomy files
    if new_rep_set is not None:
        with open(rep_set, 'r') as in_file, open(new_rep_set, 'w') as out_file:
            SeqIO.write((record for record in SeqIO.parse(in_file, 'fasta') if record.id in tax_info), out_file, 'fasta-2line')
    if new_taxonomy is not None:
        with open(new_taxonomy, 'w') as out_file:
            for line in new_taxonomy_lines:
                out_file.write(line)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers(help = 'task to be performed', dest = 'command')

    match_parser = subparsers.add_parser('match', help = 'match sequences against references')
    match_parser.add_argument('sequences', help = 'sequences to match against the references')
    match_parser.add_argument('references', help = 'reference sequences')
    match_parser.add_argument('threshold', type = float, help = 'matching threshold')
    match_parser.add_argument('output', help = 'output file')
    match_parser.add_argument('--vsearch', default = util.DEF_VSEARCH_PATH, help = 'path to VSEARCH binary')
    match_parser.add_argument('-t', '--threads', type = int, default = DEF_THREADS, help = 'number of threads')

    relabel_parser = subparsers.add_parser('relabel', help = 'relabel reference identifiers in a taxonomic assignment with the taxonomy corresponding to the identifier')
    relabel_parser.add_argument('taxonomic_assignments', help = 'taxonomic assignments to relabel')
    relabel_parser.add_argument('taxonomy', help = 'taxonomy of references')
    relabel_parser.add_argument('output', help = 'output file')

    reduce_parser = subparsers.add_parser('reduce', help = 'reduce and label FASTA reference file to those entries with clean, and sufficiently "complete" taxonomic information')
    reduce_parser.add_argument('rep_set', help = 'FASTA file of representatives')
    reduce_parser.add_argument('taxonomy', help = 'taxonomy of references')
    reduce_parser.add_argument('reference_file', help = 'reduced and relabelled reference sequences')
    reduce_parser.add_argument('--new_rep_set', help = 'reduced representatives file')
    reduce_parser.add_argument('--new_taxonomy', help = 'reduced taxonomy')
    reduce_parser.add_argument('-l', '--level', type = int, default = rst.DEF_LEVEL, help = 'highest level of taxonomic information that has to be complete')

    args = argparser.parse_args()


    if args.command == 'match':
        match_against_references(args.sequences, args.references, args.threshold, args.output, args.vsearch, args.threads)

    if args.command == 'relabel':
        relabel(args.taxonomic_assignments, args.taxonomy, args.output)

    if args.command == 'reduce':
        reduce_and_label(args.rep_set, args.taxonomy, args.new_rep_set, new_rep_set = args.new_rep_set, new_taxonomy = args.new_taxonomy, level = args.level)
