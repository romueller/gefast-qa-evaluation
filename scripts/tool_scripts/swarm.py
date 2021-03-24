#!/usr/bin/python

# Organises the measured (time & memory) execution of Swarm (v2 and higher) with different parameters / options.

import subprocess as sp

from .. import util

DEF_LOG_CMD = ''
DEF_NUM_THREADS = 1


# Run Swarm after obtaining further arguments from a string (and log runtime and memory consumption).
def run_swarm_str_args(binary, in_file, threshold, out_file, args_string, args_sep = util.MISC_ARGS_SEP, log_cmd = DEF_LOG_CMD):
    op_f = False  # fastidious
    op_n = False  # no_breaking
    op_z = False  # usearch_abundance
    op_t = DEF_NUM_THREADS  # threads
    op_c = False  # clean

    if 'fastidious' in args_string:
        op_f = True

    if 'no_breaking' in args_string:
        op_n = True

    if 'usearch_abundance' in args_string:
        op_z = True

    val = util.get_value(args_string, 'threads', args_sep)
    if val is not None:
        op_t = int(val)

    if 'clean' in args_string:
        op_c = True

    run_swarm(binary, in_file, threshold, out_file, fastidious = op_f, no_breaking = op_n,
              usearch_separator = op_z, threads = op_t, clean = op_c, log_cmd = log_cmd)


# Run Swarm (and log runtime and memory consumption).
def run_swarm(binary, in_file, threshold, out_file, fastidious = False, no_breaking = False,
              usearch_separator = False, threads = DEF_NUM_THREADS, clean = False, log_cmd = DEF_LOG_CMD):
    op_f = '-f' if fastidious and threshold == 1 else ''
    op_n = '-n' if no_breaking else ''
    op_z = '-z' if usearch_separator else ''

    sp.call('%s %s %s %s %s -d %i -o %s -a 1 -t %i %s'
            % (log_cmd, binary, op_f, op_n, op_z, threshold, out_file, threads, in_file), shell = True)
    if clean:
        sp.call('rm %s' % out_file, shell = True)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('in_file', help = 'input FASTA file')
    argparser.add_argument('threshold', type = int, help = 'clustering threshold')
    argparser.add_argument('out_file', help = 'output file')
    argparser.add_argument('-f', '--fastidious', action = 'store_true', help = 'activates fastidious clustering phase')
    argparser.add_argument('-n', '--no_breaking', action = 'store_true', help = 'deactivate breaking mechanism')
    argparser.add_argument('-z', '--usearch_abundance', action = 'store_true', help = 'accept ";size=" as abundance separator')
    argparser.add_argument('-t', '--threads', default = DEF_NUM_THREADS, help = 'number of threads')
    argparser.add_argument('-c', '--clean', action = 'store_true', help = 'remove output file')
    argparser.add_argument('--swarm', default = util.DEF_SWARM_PATH, help = 'path to Swarm binary')
    argparser.add_argument('--log_cmd', default = DEF_LOG_CMD, help = 'log command')
    args = argparser.parse_args()

    run_swarm(args.swarm, args.in_file, args.threshold, args.out_file, args.fastidious, args.no_breaking,
              args.usearch_abundance, args.threads, args.clean, args.log_cmd)
