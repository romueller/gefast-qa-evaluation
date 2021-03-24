#!/usr/bin/python

# Organises the measured (time & memory) execution of UPARSE-OTU algorithm with different parameters / options.

import subprocess as sp

from .. import util

DEF_LOG_CMD = ''


# Run UPARSE after obtaining further arguments from a string (and log runtime and memory consumption).
def run_uparse_str_args(binary, in_file, out_file, args_string, args_sep = util.MISC_ARGS_SEP, log_cmd = DEF_LOG_CMD):
    op_c = False  # clean

    if 'clean' in args_string:
        op_c = True

    run_usearch(binary, in_file, out_file, clean = op_c, log_cmd = log_cmd)


# Run UPARSE (and log runtime and memory consumption).
def run_usearch(binary, in_file, out_file, clean = False, log_cmd = DEF_LOG_CMD):
    sp.call('%s %s -cluster_otus %s -uparseout %s -fulldp -minsize 1' % (log_cmd, binary, in_file, out_file), shell = True)
    if clean:
        sp.call('rm %s' % out_file, shell = True)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('in_file', help = 'input FASTA file')
    argparser.add_argument('out_file', help = 'output file (uparseout)')
    argparser.add_argument('-c', '--clean', action = 'store_true', help = 'remove output file')
    argparser.add_argument('--usearch', default = util.DEF_USEARCH_PATH, help = 'path to USEARCH binary')
    argparser.add_argument('--log_cmd', default = DEF_LOG_CMD, help = 'log command')
    args = argparser.parse_args()

    run_usearch(args.usearch, args.in_file, args.out_file, args.clean, args.log_cmd)
