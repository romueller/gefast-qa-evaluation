#!/usr/bin/python

# Organises the measured (time & memory) execution of USEARCH with different parameters / options.

import subprocess as sp

from .. import util

DEF_LOG_CMD = ''
DEF_NUM_THREADS = 1


# Run USEARCH after obtaining further arguments from a string (and log runtime and memory consumption).
def run_usearch_str_args(binary, in_file, threshold, cluster_opt, order, out_file,
                         args_string, args_sep = util.MISC_ARGS_SEP, log_cmd = DEF_LOG_CMD):
    op_t = DEF_NUM_THREADS  # threads
    op_c = False  # clean

    val = util.get_value(args_string, 'threads', args_sep)
    if val is not None:
        op_t = int(val)

    if 'clean' in args_string:
        op_c = True

    run_usearch(binary, in_file, threshold, cluster_opt, order, out_file, threads = op_t, clean = op_c, log_cmd = log_cmd)


# Run USEARCH (and log runtime and memory consumption).
def run_usearch(binary, in_file, threshold, cluster_opt, order, out_file,
                threads = DEF_NUM_THREADS, clean = False, log_cmd = DEF_LOG_CMD):
    op_t = ('-threads %i' % threads) if cluster_opt == 'cluster_fast' else ''
    op_s = '-sort%s %s' % ('' if cluster_opt == 'cluster_fast' else 'edby', order)

    sp.call('%s %s %s -%s %s -id %f -uc %s %s -fulldp'
            % (log_cmd, binary, op_t, cluster_opt, in_file, threshold, out_file, op_s), shell = True)
    if clean:
        sp.call('rm %s' % out_file, shell = True)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('in_file', help = 'input FASTA file')
    argparser.add_argument('threshold', type = float, help = 'clustering threshold')
    argparser.add_argument('cluster_opt', help = 'clustering option (cluster_fast, cluster_smallmem)')
    argparser.add_argument('order', help = 'sorting criterion (length, size)')
    argparser.add_argument('out_file', help = 'output file')
    argparser.add_argument('-t', '--threads', default = DEF_NUM_THREADS, help = 'number of threads')
    argparser.add_argument('-c', '--clean', action = 'store_true', help = 'remove output file')
    argparser.add_argument('--usearch', default = util.DEF_USEARCH_PATH, help = 'path to USEARCH binary')
    argparser.add_argument('--log_cmd', default = DEF_LOG_CMD, help = 'log command')
    args = argparser.parse_args()

    run_usearch(args.usearch, args.in_file, args.threshold, args.cluster_opt, args.order, args.out_file,
                args.threads, args.clean, args.log_cmd)
