#!/usr/bin/python

# Organises the measured (time & memory) execution of VSEARCH with different parameters / options.

import subprocess as sp

from .. import util

DEF_LOG_CMD = ''
DEF_NUM_THREADS = 1


# Run VSEARCH after obtaining further arguments from a string (and log runtime and memory consumption).
def run_vsearch_str_args(binary, in_file, threshold, cluster_opt, out_file,
                         args_string, args_sep = util.MISC_ARGS_SEP, log_cmd = DEF_LOG_CMD):
    op_u = False
    op_t = DEF_NUM_THREADS
    op_c = False

    if 'usersort' in args_string:
        op_u = True

    val = util.get_value(args_string, 'threads', args_sep)
    if val is not None:
        op_t = int(val)

    if 'clean' in args_string:
        op_c = True

    run_vsearch(binary, in_file, threshold, cluster_opt, out_file, usersort = op_u,
                threads = op_t, clean = op_c, log_cmd = log_cmd)


# Run VSEARCH (and log runtime and memory consumption).
def run_vsearch(binary, in_file, threshold, cluster_opt, out_file, usersort = False,
                threads = DEF_NUM_THREADS, clean = False, log_cmd = DEF_LOG_CMD):
    op_us = '--usersort' if usersort else ''

    sp.call('%s %s --threads %i --%s %s --id %f --uc %s %s'
            % (log_cmd, binary, threads, cluster_opt, in_file, threshold, out_file, op_us), shell = True)
    if clean:
        sp.call('rm %s' % out_file, shell = True)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('in_file', help = 'input FASTA file')
    argparser.add_argument('threshold', type = float, help = 'clustering threshold')
    argparser.add_argument('cluster_opt', help = 'clustering option (cluster_fast, cluster_size, cluster_smallmem)')
    argparser.add_argument('out_file', help = 'output file')
    argparser.add_argument('-u', '--usersort', action = 'store_true', help = 'indicate that the sequences are sorted by another criterion')
    argparser.add_argument('-t', '--threads', default = DEF_NUM_THREADS, help = 'number of threads')
    argparser.add_argument('-c', '--clean', action = 'store_true', help = 'remove output file')
    argparser.add_argument('--vsearch', default = util.DEF_VSEARCH_PATH, help = 'path to VSEARCH binary')
    argparser.add_argument('--log_cmd', default = DEF_LOG_CMD, help = 'log command')
    args = argparser.parse_args()

    run_vsearch(args.vsearch, args.in_file, args.threshold, args.cluster_opt, args.out_file, args.usersort,
                args.threads, args.clean, args.log_cmd)
