#!/usr/bin/python

# Runs a clustering-quality and performance analysis on the given data set.
#
# Analyses the precision, recall and adjusted Rand index of the computed clusters (when a taxonomy file is provided)
# and records the runtime and memory consumption (when logging is activated).

import os
import subprocess as sp

from .. import util
from ..cluster_evaluation import compute_confusion_table as cct
from ..cluster_evaluation import evaluate_confusion_table as ect
from ..data_preparation import create_alternative_inputs as cai
from ..data_preparation import fastq2fasta as f2f
from ..tool_scripts import gefast as ex_gefast
from ..tool_scripts import swarm as ex_swarm
from ..tool_scripts import uparse as ex_uparse
from ..tool_scripts import usearch as ex_usearch
from ..tool_scripts import vsearch as ex_vsearch

DEF_R_CODE_PATH = 'scripts/tool_scripts/dada2.R'


# Compute clusters according to the tasks file and, if a taxonomy file is provided, analyse their quality.
def run_analysis(analysis_name, reads, tasks, out_dir, binaries, taxonomy_file = None,
                 log = False, log_cmd = util.DEF_LOG_CMD, log_header = util.DEF_LOG_HEADER, clear = False):

    util.check_task_sanity(reads, tasks)

    sp.call('mkdir -p %s' % out_dir, shell = True)

    if log:
        log_file = '%s/%s__performance.log' % (out_dir, analysis_name)
        sp.call('mkdir -p %s' % out_dir, shell = True)
        with open(log_file, 'w') as out_file:
            out_file.write('%s\n' % log_header)
        log_cmd = log_cmd % log_file
    else:
        log_cmd = ''


    reads_name, reads_path, reads_type = reads
    print('Processing reads %s...' % reads_name)

    # Preparation of metric files and ground truths
    metrics_files = dict()
    confusion_table = '%s/%s__conf_table' % (out_dir, analysis_name)
    if taxonomy_file is not None:
        tax_name, tax_file = taxonomy_file
        metrics_files[tax_name] = '%s/%s__metrics.csv' % (out_dir, analysis_name)
        with open(metrics_files[tax_name], 'w') as out_file:
            out_file.write('reads;gt;task;tool;mode;threshold;precision;recall;adjrandindex\n')

    # Preparation of FASTA version of input file
    reads_path_fasta = reads_path
    if reads_type == 'fastq':
        reads_path_fasta = '%s/%s.fasta.tmp' % (out_dir, reads_name)
        f2f.convert(reads_path, reads_path_fasta)

    # Clustering and evaluation
    for task, tool, mode, config, thresholds, misc in tasks:
        print('Executing task %s...' % task)

        if tool == 'gefast':
            # mode = GeFaST mode, config = source for base configuration file
            # misc: deactivate breaking, remove output file, change abundance separator
            refinement_thresholds = util.get_value(misc, 'refinement_threshold', util.MISC_ARGS_SEP)
            if refinement_thresholds is None:
                refinement_thresholds = [None] * len(thresholds)
            else:
                refinement_thresholds = refinement_thresholds.split(',,')  # two commas to distinguish it from the list notation of a single run
                if len(refinement_thresholds) == 1:
                    refinement_thresholds = [refinement_thresholds[0]] * len(thresholds)

            for i, t in enumerate(thresholds):
                if log:
                    with open(log_file, 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%f;' % (reads_name, task, tool, mode, t))

                otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
                ex_gefast.run_gefast_str_args(binaries[tool], mode, reads_path, config, t, otus_file, misc,
                                              args_sep = util.MISC_ARGS_SEP, log_cmd = log_cmd,
                                              refinement_threshold = refinement_thresholds[i])

                if taxonomy_file is not None:
                    tax_name, tax_file = taxonomy_file
                    cct.compute_table(tax_file, otus_file, confusion_table)
                    metric_values = ect.evaluate(confusion_table)
                    with open(metrics_files[tax_name], 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%s;%f;%s\n' % (
                            reads_name, tax_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

                if clear:
                    os.remove(otus_file)

        if tool == 'swarm':
            # mode and config do not apply
            # misc: activate fastidious, deactivate breaking, remove output file, change abundance separator
            derep = ('derep' in misc)
            if derep:
                alt_reads_file = '%s/%s_%s_alt_reads_dereplicated.fasta' % (out_dir, analysis_name, task)
                cai.create_file(binaries['vsearch'], reads_path_fasta, '%s.tmp' % alt_reads_file, derep = derep)
                cai.create_file(binaries['vsearch'], '%s.tmp' % alt_reads_file, alt_reads_file, old_sep = ';size=', new_sep = '_')
                os.remove('%s.tmp' % alt_reads_file)
                reads_path_swarm = alt_reads_file
            else:
                reads_path_swarm = reads_path_fasta

            for t in thresholds:
                if log:
                    with open(log_file, 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%f;' % (reads_name, task, tool, mode, t))

                otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
                ex_swarm.run_swarm_str_args(binaries[tool], reads_path_swarm, t, otus_file, misc,
                                            args_sep = util.MISC_ARGS_SEP, log_cmd = log_cmd)

                if taxonomy_file is not None:
                    tax_name, tax_file = taxonomy_file
                    cct.compute_table(tax_file, otus_file, confusion_table)
                    metric_values = ect.evaluate(confusion_table)
                    with open(metrics_files[tax_name], 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%s;%f;%s\n' % (
                            reads_name, tax_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

                if clear:
                    os.remove(otus_file)

            if derep:
                os.remove(alt_reads_file)

        if tool == 'usearch':
            # mode is used for cluster_opt
            # config does not apply
            # misc: order, number of threads, remove output file; for alternative input: minsize, old / new abundance separator
            cluster_opt = mode
            order = util.get_value(misc, 'order', util.MISC_ARGS_SEP)
            min_size = util.get_value(misc, 'minsize', util.MISC_ARGS_SEP)
            min_size = cai.DEF_MIN_SIZE if min_size is None else int(min_size)
            old_sep = util.get_value(misc, 'old_sep', util.MISC_ARGS_SEP)
            old_sep = cai.DEF_OLD_SEP if old_sep is None else old_sep
            new_sep = util.get_value(misc, 'new_sep', util.MISC_ARGS_SEP)
            new_sep = cai.DEF_NEW_SEP if new_sep is None else new_sep

            alt_reads_file = '%s/%s_%s_alt_reads_usearch.fasta' % (out_dir, analysis_name, task)
            cai.create_file(binaries['vsearch'], reads_path_fasta, alt_reads_file, sort_criterion = order,
                            min_size = min_size, old_sep = old_sep, new_sep = new_sep)

            for t in thresholds:
                if log:
                    with open(log_file, 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%f;' % (reads_name, task, tool, mode, t))

                otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
                ex_usearch.run_usearch_str_args(binaries[tool], alt_reads_file, t, cluster_opt, order, otus_file,
                                                misc, args_sep = util.MISC_ARGS_SEP, log_cmd = log_cmd)

                # transform output file from UCLUST format to list-style output of Swarm / GeFaST
                clusters = []
                with open(otus_file, 'r') as in_file:
                    for line in in_file:
                        row_type, cluster_num, _, _, _, _, _, _, seq_label, _ = line.strip().split('\t')
                        if row_type == 'S':
                            clusters.append([seq_label.replace(';size=', '_').rstrip(';')])
                        elif row_type == 'H':
                            clusters[int(cluster_num)].append(seq_label.replace(';size=', '_').rstrip(';'))
                with open(otus_file, 'w') as out_file:
                    for c in clusters:
                        out_file.write('%s\n' % ' '.join(c))

                if taxonomy_file is not None:
                    tax_name, tax_file = taxonomy_file
                    cct.compute_table(tax_file, otus_file, confusion_table)
                    metric_values = ect.evaluate(confusion_table)
                    with open(metrics_files[tax_name], 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%s;%f;%s\n' % (
                            reads_name, tax_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

                if clear:
                    os.remove(otus_file)

            os.remove(alt_reads_file)

        if tool == 'vsearch':
            # mode is used for cluster_opt
            # config does not apply
            # misc: order, number of threads, minsize, remove output files
            cluster_opt = mode
            order = util.get_value(misc, 'order', util.MISC_ARGS_SEP)
            min_size = util.get_value(misc, 'minsize', util.MISC_ARGS_SEP)
            min_size = cai.DEF_MIN_SIZE if min_size is None else int(min_size)
            old_sep = util.get_value(misc, 'old_sep', util.MISC_ARGS_SEP)
            old_sep = cai.DEF_OLD_SEP if old_sep is None else old_sep
            new_sep = util.get_value(misc, 'new_sep', util.MISC_ARGS_SEP)
            new_sep = cai.DEF_NEW_SEP if new_sep is None else new_sep

            alt_reads_file = '%s/%s_%s_alt_reads_vsearch.fasta' % (out_dir, analysis_name, task)
            cai.create_file(binaries['vsearch'], reads_path_fasta, alt_reads_file, sort_criterion = order,
                            min_size = min_size, old_sep = old_sep, new_sep = new_sep)

            for t in thresholds:
                if log:
                    with open(log_file, 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%f;' % (reads_name, task, tool, mode, t))

                otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
                ex_vsearch.run_vsearch_str_args(binaries[tool], alt_reads_file, t, cluster_opt, otus_file, misc,
                                                args_sep = util.MISC_ARGS_SEP, log_cmd = log_cmd)

                # transform output file from UCLUST format to list-style output of Swarm / GeFaST
                clusters = []
                with open(otus_file, 'r') as in_file:
                    for line in in_file:
                        row_type, cluster_num, _, _, _, _, _, _, seq_label, _ = line.strip().split('\t')
                        if row_type == 'S':
                            clusters.append([seq_label.replace(';size=', '_')])
                        elif row_type == 'H':
                            clusters[int(cluster_num)].append(seq_label.replace(';size=', '_'))
                with open(otus_file, 'w') as out_file:
                    for c in clusters:
                        out_file.write('%s\n' % ' '.join(c))

                if taxonomy_file is not None:
                    tax_name, tax_file = taxonomy_file
                    cct.compute_table(tax_file, otus_file, confusion_table)
                    metric_values = ect.evaluate(confusion_table)
                    with open(metrics_files[tax_name], 'a') as out_file:
                        out_file.write('%s;%s;%s;%s;%s;%f;%s\n' % (
                            reads_name, tax_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

                if clear:
                    os.remove(otus_file)

            os.remove(alt_reads_file)

        if tool == 'uparse':
            # mode and config do not apply
            # misc: remove output file; for alternative input: minsize, old / new abundance separator
            order = util.get_value(misc, 'order', util.MISC_ARGS_SEP)
            min_size = util.get_value(misc, 'minsize', util.MISC_ARGS_SEP)
            min_size = cai.DEF_MIN_SIZE if min_size is None else int(min_size)
            old_sep = util.get_value(misc, 'old_sep', util.MISC_ARGS_SEP)
            old_sep = cai.DEF_OLD_SEP if old_sep is None else old_sep
            new_sep = util.get_value(misc, 'new_sep', util.MISC_ARGS_SEP)
            new_sep = cai.DEF_NEW_SEP if new_sep is None else new_sep

            alt_reads_file = '%s/%s_%s_alt_reads_uparse.fasta' % (out_dir, analysis_name, task)
            cai.create_file(binaries['vsearch'], reads_path_fasta, alt_reads_file, sort_criterion = order,
                            min_size = min_size, old_sep = old_sep, new_sep = new_sep, derep = ('derep' in misc))

            t = 0.97 # UPARSE-OTU uses a fixed threshold
            if log:
                with open(log_file, 'a') as out_file:
                    out_file.write('%s;%s;%s;%s;%f;' % (reads_name, task, tool, mode, t))

            otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
            ex_uparse.run_uparse_str_args(binaries[tool], alt_reads_file, otus_file, misc,
                                          args_sep = util.MISC_ARGS_SEP, log_cmd = log_cmd)

            # transform uparseout output file from tabbed classification format to list-style output of Swarm / GeFaST
            clusters = []
            cluster_id_map = dict()
            chimeras = set()
            with open(otus_file, 'r') as in_file:
                for line in in_file:
                    seq_label, classification, info = line.strip().split('\t')
                    seq_id = seq_label.split(';size=')[0]
                    if classification.startswith('otu'): # representative of new OTU
                        clusters.append([seq_label.replace(';size=', '_')])
                        cluster_id_map[seq_id] = len(cluster_id_map)

                    elif classification in ['noisy_chimera', 'perfect_chimera']: # chimeric entry
                        chimeras.add(seq_id)

                    elif classification in ['match', 'perfect']: # similar to previous OTU (or chimera)
                        # determine top-entry (if any) and the contained identifier (w/o abundance and 'top=')
                        target = [token for token in info.split(';') if token.startswith('top=')]
                        target_id = None if len(target) == 0 else target[0][4:]

                        if target_id in cluster_id_map: # add to cluster of similar representative
                            clusters[cluster_id_map[target_id]].append(seq_label.replace(';size=', '_'))
                        elif target_id in chimeras: # record as chimera too
                            chimeras.add(seq_id)
                            # chimeras.add(seq_label.replace(';size=', '_'))
                        else:
                            print('Unexpected entry without top: %s' % line)

                    else:
                        print('Unexpected classification: %s' % line)

            with open(otus_file, 'w') as out_file:
                for c in clusters:
                    out_file.write('%s\n' % ' '.join(c))

            if taxonomy_file is not None:
                tax_name, tax_file = taxonomy_file
                cct.compute_table(tax_file, otus_file, confusion_table)
                metric_values = ect.evaluate(confusion_table)
                with open(metrics_files[tax_name], 'a') as out_file:
                    out_file.write('%s;%s;%s;%s;%s;%f;%s\n' % (
                        reads_name, tax_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

            if clear:
                os.remove(otus_file)

            os.remove(alt_reads_file)

    if reads_type == 'fastq':
        os.remove(reads_path_fasta)

    if taxonomy_file is not None:
        os.remove(confusion_table)


# Compute clusters using DADA2 and, if taxonomic assignments are provided, analyse their quality.
# Should only be executed after run_analysis(...).
def run_dada2_analysis(analysis_name, forward_reads, reverse_reads, out_dir, paths, taxonomy_file = None,
                       log = False, log_cmd = util.DEF_LOG_CMD, clear = False):

    sp.call('mkdir -p %s' % out_dir, shell = True)

    if log:
        log_file = '%s/%s__performance.log' % (out_dir, analysis_name)
        log_cmd = log_cmd % log_file
    else:
        log_cmd = ''

    out_prefix = '%s/%s' % (out_dir, analysis_name)
    metrics_files = dict()
    confusion_table = '%s/%s__conf_table' % (out_dir, analysis_name)
    if taxonomy_file is not None:
        tax_name, tax_file = taxonomy_file
        metrics_files[tax_name] = '%s/%s_%s__metrics.csv' % (out_dir, analysis_name, tax_name)


    if log:
        with open(log_file, 'a') as out_file:
            out_file.write('%s;dada2;dada2;;0;' % analysis_name)

    if reverse_reads is None:
        sp.call('%s R -e \'source("%s"); dada2_cluster_single("%s", "%s")\''
                % (log_cmd, paths['r_code'], forward_reads, out_prefix), shell = True)

    else:
        sp.call('%s R -e \'source("%s"); dada2_cluster_paired("%s", "%s", "%s")\''
                % (log_cmd, paths['r_code'], forward_reads, reverse_reads, out_prefix), shell = True)

    if taxonomy_file is not None:
        # with chimeras
        otus_file = '%s_dada2_otus.txt' % out_prefix
        tax_name, tax_file = taxonomy_file
        cct.compute_table(tax_file, otus_file, confusion_table)
        metric_values = ect.evaluate(confusion_table)
        with open(metrics_files[tax_name], 'a') as out_file:
            out_file.write('%s;%s;dada2;dada2;;0;%s\n' % (analysis_name, tax_name, ';'.join(map(str, metric_values))))
        if clear:
            os.remove(otus_file)
        # without chimeras
        otus_file = '%s_dada2_nc_otus.txt' % out_prefix
        tax_name, tax_file = taxonomy_file
        cct.compute_table(tax_file, otus_file, confusion_table)
        metric_values = ect.evaluate(confusion_table)
        with open(metrics_files[tax_name], 'a') as out_file:
            out_file.write('%s;%s;dada2-nc;dada2;;0;%s\n' % (analysis_name, tax_name, ';'.join(map(str, metric_values))))
        if clear:
            os.remove(otus_file)
        os.remove(confusion_table)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers(help = 'task to be performed', dest = 'command')

    run_parser = subparsers.add_parser('run', help = 'execute analysis')
    run_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_parser.add_argument('reads_file', help = 'name of read collection and path to file, separated by a colon (e.g. file1:path/to/file)')
    run_parser.add_argument('tasks_file', help = 'list file of tasks')
    run_parser.add_argument('out_dir', help = 'output directory')
    run_parser.add_argument('--tax_file', help = 'pair consisting of name of and path to taxonomy file (components of pair separated by colon; activates analysis of clustering quality)')
    run_parser.add_argument('--gefast', default = util.DEF_GEFAST_PATH, help = 'path to GeFaST binary')
    run_parser.add_argument('--swarm', default = util.DEF_SWARM_PATH, help = 'path to Swarm (v2 and higher) binary')
    run_parser.add_argument('--usearch', default = util.DEF_USEARCH_PATH, help = 'path to USEARCH binary')
    run_parser.add_argument('--vsearch', default = util.DEF_VSEARCH_PATH, help = 'path to VSEARCH binary')
    run_parser.add_argument('--log', action = 'store_true', help = 'activate performance logging')
    run_parser.add_argument('--log_cmd', default = util.DEF_LOG_CMD, help = 'logging command')
    run_parser.add_argument('--log_file_header', default = util.DEF_LOG_HEADER, help = 'header of log file')
    run_parser.add_argument('--clear', action = 'store_true', help = 'remove clustering outputs')

    run_dada2_parser = subparsers.add_parser('run_dada2', help = 'execute analysis with DADA2')
    run_dada2_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_dada2_parser.add_argument('forward_reads', help = 'FASTQ file with (preprocessed) forward reads')
    run_dada2_parser.add_argument('out_dir', help = 'output directory')
    run_dada2_parser.add_argument('-r', '--reverse_reads', help = 'FASTQ file with (preprocessed) reverse reads')
    run_dada2_parser.add_argument('--tax_file', help = 'pair consisting of name of and path to taxonomy file (components of pair separated by colon; activates analysis of clustering quality)')
    run_dada2_parser.add_argument('--r_code', default = DEF_R_CODE_PATH, help = 'path to underlying R code')
    run_dada2_parser.add_argument('--log', action = 'store_true', help = 'activate performance logging')
    run_dada2_parser.add_argument('--log_cmd', default = util.DEF_LOG_CMD, help = 'logging command')
    run_dada2_parser.add_argument('--clear', action='store_true', help='remove clustering outputs')

    args = argparser.parse_args()


    if args.command == 'run':
        reads = util.parse_reads_list(args.reads_file)[0]
        tasks = util.parse_tasks_list_file(args.tasks_file)
        tax_file = args.tax_file.split(':') if args.tax_file is not None else None

        binaries = {'gefast': args.gefast, 'swarm': args.swarm, 'uparse': args.usearch, 'usearch': args.usearch, 'vsearch': args.vsearch}

        run_analysis(args.analysis_name, reads, tasks, args.out_dir, binaries, taxonomy_file = tax_file,
                     log = args.log, log_cmd = args.log_cmd, log_header = args.log_file_header, clear = args.clear)


    if args.command == 'run_dada2':
        tax_file = args.tax_file.split(':') if args.tax_file is not None else None

        paths = dict()
        paths['r_code'] = args.r_code

        run_dada2_analysis(args.analysis_name, args.forward_reads, args.reverse_reads, args.out_dir, paths, taxonomy_file = tax_file,
                           log = args.log, log_cmd = args.log_cmd, clear = args.clear)
