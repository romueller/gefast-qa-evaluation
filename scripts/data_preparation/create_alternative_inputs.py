#!/usr/bin/python

# Creates input files for USEARCH and VSEARCH (with VSEARCH to avoid memory limit of 32-bit USEARCH).

import subprocess as sp

from .. import util

DEF_OLD_SEP = '_'
DEF_NEW_SEP = ';size='
DEF_MIN_SIZE = 1


# Create alternative input files.
def create_file(binary, in_file, out_file, derep = False, sort_criterion = None, min_size = DEF_MIN_SIZE,
                old_sep = DEF_OLD_SEP, new_sep = DEF_NEW_SEP):
    sp.call("sed 's/%s/%s/' %s > %s" % (old_sep, new_sep, in_file, out_file), shell = True)

    if derep:
        sp.call('mv %s %s.tmp; %s --derep_fulllength %s.tmp --sizein --fasta_width 0 --sizeout --output %s --minuniquesize %i; rm %s.tmp'
                % (out_file, out_file, binary, out_file, out_file, min_size, out_file), shell = True)

    if sort_criterion == 'length':
        sp.call('mv %s %s.tmp; %s --sortbylength %s.tmp --output %s --sizeout; rm %s.tmp'
                % (out_file, out_file, binary, out_file, out_file, out_file), shell = True)

    if sort_criterion == 'size':
        sp.call('mv %s %s.tmp; %s --sortbysize %s.tmp --output %s --minsize %i --sizeout; rm %s.tmp'
                % (out_file, out_file, binary, out_file, out_file, min_size, out_file), shell = True)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('in_file', help = 'input FASTA file')
    argparser.add_argument('out_file', help = 'output file')
    argparser.add_argument('--old_sep', default = DEF_OLD_SEP, help = 'old abundance separator to be replaced')
    argparser.add_argument('--new_sep', default = DEF_NEW_SEP, help = 'new abundance separator')
    argparser.add_argument('--minsize', type = int, default = DEF_MIN_SIZE, help = 'minimum sequence abundance')
    argparser.add_argument('--dereplicate', action = 'store_true', help = 'dereplicate the sequences before sorting them')
    argparser.add_argument('--vsearch', default = util.DEF_VSEARCH_PATH, help = 'path to VSEARCH binary')
    sorting = argparser.add_mutually_exclusive_group(required = False)
    sorting.add_argument('--sortbylength', action = 'store_true', help = 'sort the sequences by their length')
    sorting.add_argument('--sortbysize', action = 'store_true', help = 'sort the sequences by their abundance')

    args = argparser.parse_args()


    sort_criterion = None
    if args.sortbylength:
        sort_criterion = 'length'
    if args.sortbysize:
        sort_criterion = 'size'

    create_file(args.vsearch, args.in_file, args.out_file, args.derep, sort_criterion, args.minsize, args.old_sep, args.new_sep)
