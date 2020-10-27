#!/usr/bin/python

# Prepares and runs a clustering-quality analysis on different mock-community data sets
# similar to the analysis in the following paper:
#
# Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. (2016)
# DADA2: High-resolution sample inference from Illumina amplicon data.
# https://doi.org/10.1038/nmeth.3869
#
# Analyses the precision, recall and adjusted Rand index of the computed clusters.
#
# Offers three kinds of analysis runs:
# (1) GeFaST with its diverse modes to evaluate the impact of quality-weighted alignments.
# (2) DADA2 (used as a clustering tool, for comparison with DADA2-inspired GeFaST clustering / refinement approaches)
# (3) USEARCH and VSEARCH (for comparison with DADA2-inspired GeFaST clustering / refinement approaches)

import os
import subprocess as sp

from Bio import SeqIO

from .. import util
from ..cluster_evaluation import compute_confusion_table as cct
from ..cluster_evaluation import evaluate_confusion_table as ect
from ..data_preparation import assign_taxa as at
from ..data_preparation import create_alternative_inputs as cai
from ..data_preparation import fastq2fasta as f2f
from ..data_preparation import reformat_silva_taxonomy as rst
from ..tool_scripts import gefast as ex_gefast
from ..tool_scripts import usearch as ex_usearch
from ..tool_scripts import vsearch as ex_vsearch


DEF_R_CODE_PATH = 'scripts/tool_scripts/dada2.R'
DEF_USEARCH8_PATH = 'usearch8'
DEF_USEARCH10_PATH = 'usearch10'


def prepare_balanced(out_dir, paths, single = False, max_n = None, max_ee = None, trunc_q = None, trunc_len = None,
                     trim_left = None, min_size = None, fastq_min_ovlen = None, fastq_max_diffs = None, clean = False):

    sp.call('mkdir -p %s' % out_dir, shell = True)

    # get default values when for unspecified arguments
    max_n = 0 if max_n is None else max_n
    max_ee = 2 if max_ee is None else max_ee
    trunc_q = 2 if trunc_q is None else trunc_q
    trunc_len = "240,220" if trunc_len is None else trunc_len
    trim_left = "20,20" if trim_left is None else trim_left
    min_size = 2 if min_size is None else min_size
    fastq_min_ovlen = 20 if fastq_min_ovlen is None else fastq_min_ovlen
    fastq_max_diffs = 1 if fastq_max_diffs is None else fastq_max_diffs

    # download
    sp.call('wget -P %s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR777/ERR777695/ERR777695_1.fastq.gz' % out_dir, shell = True)
    if not single:
        sp.call('wget -P %s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR777/ERR777695/ERR777695_2.fastq.gz' % out_dir, shell = True)

    forward_reads = '%s/ERR777695_1.fastq.gz' % out_dir
    reverse_reads = '%s/ERR777695_2.fastq.gz' % out_dir if not single else ''

    # filtering (e.g. length, quality)
    if single:
        sp.call('R -e \'source("%s"); filter_single("%s", %i, %i, %i, %s, %s)\''
                % (paths['r_code'], forward_reads, max_n, max_ee, trunc_q, trunc_len.split(',')[0], trim_left.split(',')[0]), shell = True)
    else:
        sp.call('R -e \'source("%s"); filter_paired("%s", "%s", %i, %i, %i, c(%s), c(%s))\''
                % (paths['r_code'], forward_reads, reverse_reads, max_n, max_ee, trunc_q, trunc_len, trim_left), shell = True)

    sp.call('rm %s %s' % (forward_reads, reverse_reads), shell = True)

    # dereplicate and change size annotation
    if single:
        sp.call('gzip -d %s/ERR777695_1_sf.fastq.gz' % out_dir, shell = True)

        filtered_reads = '%s/ERR777695_1_sf.fastq' % out_dir
        dereplicated_reads = '%s/ERR777695_1_sfd.fastq' % out_dir

        # USEARCH 9+ modifies the quality scores (of singleton sequences) during dereplication in an inconsistent way
        # -> back to USEARCH 8 (as in the analysis files of DADA2)
        sp.call('%s -derep_fulllength %s -fastqout %s.tmp -sizeout -minuniquesize %i'
                % (paths['usearch8'], filtered_reads, dereplicated_reads, min_size), shell = True)

        sp.call("cat %s.tmp | sed '/;size=/s/;$//' | sed 's/;size=/_/' > %s" % (dereplicated_reads, dereplicated_reads), shell = True)
        sp.call('rm %s %s.tmp' % (filtered_reads if clean else '', dereplicated_reads), shell = True)

    else:
        sp.call('gzip -d %s/ERR777695_1_pf.fastq.gz %s/ERR777695_2_pf.fastq.gz' % (out_dir, out_dir), shell = True)

        filt_forward_reads = '%s/ERR777695_1_pf.fastq' % out_dir
        filt_reverse_reads = '%s/ERR777695_2_pf.fastq' % out_dir
        merged_reads = '%s/ERR777695_pfm.fastq' % out_dir
        dereplicated_reads = '%s/ERR777695_pfmd.fastq' % out_dir

        # USEARCH 8 showed merging problem on at least one data set -> USEARCH 10
        sp.call('%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastq_minovlen %i -fastq_maxdiffs %i' %
                (paths['usearch10'], filt_forward_reads, filt_reverse_reads, merged_reads, fastq_min_ovlen, fastq_max_diffs), shell = True)
        # USEARCH 9+ modifies the quality scores (of singleton sequences) during dereplication in an inconsistent way
        # -> back to USEARCH 8 (as in the analysis files of DADA2)
        sp.call('%s -derep_fulllength %s -fastqout %s.tmp -sizeout -minuniquesize %i'
                % (paths['usearch8'], merged_reads, dereplicated_reads, min_size), shell = True)

        sp.call("cat %s.tmp | sed '/;size=/s/;$//' | sed 's/;size=/_/' > %s" % (dereplicated_reads, dereplicated_reads), shell = True)
        sp.call('rm %s %s %s %s.tmp' % (filt_forward_reads if clean else '', filt_reverse_reads if clean else '',
                                        merged_reads if clean else '', dereplicated_reads), shell = True)


def prepare_hmp(out_dir, paths, single = False, max_n = None, max_ee = None, trunc_q = None, trunc_len = None,
                     trim_left = None, min_size = None, fastq_min_ovlen = None, fastq_max_diffs = None, clean = False):

    sp.call('mkdir -p %s' % out_dir, shell = True)

    # get default values when for unspecified arguments
    max_n = 0 if max_n is None else max_n
    max_ee = 2 if max_ee is None else max_ee
    trunc_q = 2 if trunc_q is None else trunc_q
    trunc_len = "240,200" if trunc_len is None else trunc_len
    trim_left = "20,20" if trim_left is None else trim_left
    min_size = 2 if min_size is None else min_size
    fastq_min_ovlen = 20 if fastq_min_ovlen is None else fastq_min_ovlen
    fastq_max_diffs = 1 if fastq_max_diffs is None else fastq_max_diffs

    # download
    sp.call('wget -P %s https://mothur.s3.us-east-2.amazonaws.com/data/MiSeqDevelopmentData/130403.tar' % out_dir, shell = True)
    sp.call('tar -xvf %s/130403.tar -C %s Mock1_S1_L001_R1_001.fastq.bz2; ' % (out_dir, out_dir)
            + 'bunzip2 -c < %s/Mock1_S1_L001_R1_001.fastq.bz2 | gzip -c > %s/Mock1_S1_L001_R1_001.fastq.gz'
            % (out_dir, out_dir), shell = True)
    if not single:
        sp.call('tar -xvf %s/130403.tar -C %s Mock1_S1_L001_R2_001.fastq.bz2; ' % (out_dir, out_dir)
                + 'bunzip2 -c < %s/Mock1_S1_L001_R2_001.fastq.bz2 | gzip -c > %s/Mock1_S1_L001_R2_001.fastq.gz'
                % (out_dir, out_dir), shell = True)

    sp.call('rm %s/130403.tar %s/Mock1_S1_L001_R1_001.fastq.bz2 %s/Mock1_S1_L001_R2_001.fastq.bz2'
            % (out_dir, out_dir, out_dir), shell=True)

    forward_reads = '%s/Mock1_S1_L001_R1_001.fastq.gz' % out_dir
    reverse_reads = '%s/Mock1_S1_L001_R2_001.fastq.gz' % out_dir if not single else ''

    # filtering (e.g. length, quality)
    if single:
        sp.call('R -e \'source("%s"); filter_single("%s", %i, %i, %i, %s, %s)\''
                % (paths['r_code'], forward_reads, max_n, max_ee, trunc_q, trunc_len.split(',')[0], trim_left.split(',')[0]), shell = True)
    else:
        sp.call('R -e \'source("%s"); filter_paired("%s", "%s", %i, %i, %i, c(%s), c(%s))\''
                % (paths['r_code'], forward_reads, reverse_reads, max_n, max_ee, trunc_q, trunc_len, trim_left), shell = True)

    sp.call('rm %s %s' % (forward_reads, reverse_reads), shell = True)

    # dereplicate and change size annotation
    if single:
        sp.call('gzip -d %s/Mock1_S1_L001_R1_001_sf.fastq.gz' % out_dir, shell = True)

        filtered_reads = '%s/Mock1_S1_L001_R1_001_sf.fastq' % out_dir
        dereplicated_reads = '%s/Mock1_S1_L001_R1_001_sfd.fastq' % out_dir

        # USEARCH 9+ modifies the quality scores (of singleton sequences) during dereplication in an inconsistent way
        # -> back to USEARCH 8 (as in the analysis files of DADA2)
        sp.call('%s -derep_fulllength %s -fastqout %s.tmp -sizeout -minuniquesize %i'
                % (paths['usearch8'], filtered_reads, dereplicated_reads, min_size), shell = True)

        sp.call("cat %s.tmp | sed '/;size=/s/;$//' | sed 's/;size=/_/' > %s" % (dereplicated_reads, dereplicated_reads), shell = True)
        sp.call('rm %s %s.tmp' % (filtered_reads if clean else '', dereplicated_reads), shell = True)

    else:
        sp.call('gzip -d %s/Mock1_S1_L001_R1_001_pf.fastq.gz %s/Mock1_S1_L001_R2_001_pf.fastq.gz' % (out_dir, out_dir), shell = True)

        filt_forward_reads = '%s/Mock1_S1_L001_R1_001_pf.fastq' % out_dir
        filt_reverse_reads = '%s/Mock1_S1_L001_R2_001_pf.fastq' % out_dir
        merged_reads = '%s/Mock1_S1_L001_pfm.fastq' % out_dir
        dereplicated_reads = '%s/Mock1_S1_L001_pfmd.fastq' % out_dir

        # USEARCH 8 showed merging problem on at least one data set -> USEARCH 10
        sp.call('%s -fastq_mergepairs %s -reverse %s -fastqout %s -fastq_minovlen %i -fastq_maxdiffs %i' %
                (paths['usearch10'], filt_forward_reads, filt_reverse_reads, merged_reads, fastq_min_ovlen, fastq_max_diffs), shell = True)
        # USEARCH 9+ modifies the quality scores (of singleton sequences) during dereplication in an inconsistent way
        # -> back to USEARCH 8 (as in the analysis files of DADA2)
        sp.call('%s -derep_fulllength %s -fastqout %s.tmp -sizeout -minuniquesize %i'
                % (paths['usearch8'], merged_reads, dereplicated_reads, min_size), shell = True)

        sp.call("cat %s.tmp | sed '/;size=/s/;$//' | sed 's/;size=/_/' > %s" % (dereplicated_reads, dereplicated_reads), shell = True)
        sp.call('rm %s %s %s %s.tmp' % (filt_forward_reads if clean else '', filt_reverse_reads if clean else '',
                                        merged_reads if clean else '', dereplicated_reads), shell = True)


# Compute clusters using GeFaST according to the tasks file and,
# if taxonomic assignments are provided, analyse their quality.
def run_analysis(analysis_name, reads, tasks, out_dir, binaries, tax_files = None):

    util.check_task_sanity(reads, tasks)

    sp.call('mkdir -p %s' % out_dir, shell = True)
    reads_name, reads_path, reads_type = reads

    print('Processing reads %s...' % reads_name)

    # Preparation of metric files and ground truths
    metrics_files = dict()
    confusion_table = '%s/%s__conf_table' % (out_dir, analysis_name)
    if tax_files is not None:
        for tax_name, tax_file in tax_files:
            metrics_files[tax_name] = '%s/%s_%s__metrics.csv' % (out_dir, reads_name, tax_name)
            with open(metrics_files[tax_name], 'w') as out_file:
                out_file.write('reads;gt;task;tool;mode;threshold;precision;recall;adjrandindex\n')

    for task, tool, mode, config, thresholds, misc in tasks:
        print('Executing task %s...' % task)

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
            otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
            ex_gefast.run_gefast_str_args(binaries[tool], mode, reads_path, config, t, otus_file, misc,
                                          args_sep = util.MISC_ARGS_SEP, refinement_threshold = refinement_thresholds[i])
            if tax_files is not None:
                for tax_name, tax_file in tax_files:
                    cct.compute_table(tax_file, otus_file, confusion_table)
                    metric_values = ect.evaluate(confusion_table)
                    with open(metrics_files[tax_name], 'a') as out_file:
                        out_file.write('%s;0;%s;%s;%s;%f;%s\n' % (reads_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

    if tax_files is not None:
        os.remove(confusion_table)


# Compute clusters using DADA2 and, if taxonomic assignments are provided, analyse their quality.
def run_dada2_analysis(analysis_name, forward_reads, reverse_reads, out_dir, paths, tax_files = None):

    sp.call('mkdir -p %s' % out_dir, shell = True)

    out_prefix = '%s/%s' % (out_dir, analysis_name)
    metrics_files = dict()
    confusion_table = '%s/%s__conf_table' % (out_dir, analysis_name)
    if tax_files is not None:
        for tax_name, tax_file in tax_files:
            metrics_files[tax_name] = '%s/%s_%s__metrics.csv' % (out_dir, analysis_name, tax_name)
            with open(metrics_files[tax_name], 'w') as out_file:
                out_file.write('reads;gt;task;tool;mode;threshold;precision;recall;adjrandindex\n')

    if reverse_reads is None:
        sp.call('R -e \'source("%s"); dada2_cluster_single("%s", "%s")\''
                % (paths['r_code'], forward_reads, out_prefix), shell = True)

    else:
        sp.call('R -e \'source("%s"); dada2_cluster_paired("%s", "%s", "%s")\''
                % (paths['r_code'], forward_reads, reverse_reads, out_prefix), shell = True)

    if tax_files is not None:
        # with chimeras
        otus_file = '%s_dada2_otus.txt' % out_prefix
        for tax_name, tax_file in tax_files:
            cct.compute_table(tax_file, otus_file, confusion_table)
            metric_values = ect.evaluate(confusion_table)
            with open(metrics_files[tax_name], 'a') as out_file:
                out_file.write('%s;0;dada2;dada2;;0;%s\n' % (analysis_name, ';'.join(map(str, metric_values))))
        # without chimeras
        otus_file = '%s_dada2_nc_otus.txt' % out_prefix
        for tax_name, tax_file in tax_files:
            cct.compute_table(tax_file, otus_file, confusion_table)
            metric_values = ect.evaluate(confusion_table)
            with open(metrics_files[tax_name], 'a') as out_file:
                out_file.write('%s;0;dada2-nc;dada2;;0;%s\n' % (analysis_name, ';'.join(map(str, metric_values))))
        os.remove(confusion_table)


# Compute clusters using USEARCH and / or VSEARCH and, if taxonomic assignments are provided, analyse their quality.
def run_uvsearch_analysis(analysis_name, reads, tasks, out_dir, binaries, tax_files = None):

    util.check_task_sanity(reads, tasks)

    sp.call('mkdir -p %s' % out_dir, shell = True)
    reads_name, reads_path, reads_type = reads

    print('Processing reads %s...' % reads_name)

    # Preparation of metric files and ground truths
    metrics_files = dict()
    confusion_table = '%s/%s__conf_table' % (out_dir, analysis_name)
    if tax_files is not None:
        for tax_name, tax_file in tax_files:
            metrics_files[tax_name] = '%s/%s_%s__metrics.csv' % (out_dir, reads_name, tax_name)
            with open(metrics_files[tax_name], 'w') as out_file:
                out_file.write('reads;gt;task;tool;mode;threshold;precision;recall;adjrandindex\n')

    # Preparation of FASTA version of innput file
    reads_path_fasta = reads_path
    if reads_type == 'fastq':
        reads_path_fasta = '%s/%s.fasta.tmp' % (out_dir, reads_name)
        f2f.convert(reads_path, reads_path_fasta)

    for task, tool, mode, config, thresholds, misc in tasks:
        print('Executing task %s...' % task)

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

            alt_reads_file = '%s/%s_%s_alt_reads_usearch.fasta.tmp' % (out_dir, analysis_name, task)
            cai.create_file(binaries['vsearch'], reads_path_fasta, alt_reads_file, sort_criterion = order,
                            min_size = min_size, old_sep = old_sep, new_sep = new_sep)

            for t in thresholds:
                otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
                ex_usearch.run_usearch_str_args(binaries[tool], alt_reads_file, t, cluster_opt, order, otus_file,
                                                misc, args_sep = util.MISC_ARGS_SEP)

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

                if tax_files is not None:
                    for tax_name, tax_file in tax_files:
                        cct.compute_table(tax_file, otus_file, confusion_table)
                        metric_values = ect.evaluate(confusion_table)
                        with open(metrics_files[tax_name], 'a') as out_file:
                            out_file.write('%s;0;%s;%s;%s;%f;%s\n' % (
                            reads_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

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
                otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
                ex_vsearch.run_vsearch_str_args(binaries[tool], alt_reads_file, t, cluster_opt, otus_file, misc,
                                                args_sep = util.MISC_ARGS_SEP)

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

                if tax_files is not None:
                    for tax_name, tax_file in tax_files:
                        cct.compute_table(tax_file, otus_file, confusion_table)
                        metric_values = ect.evaluate(confusion_table)
                        with open(metrics_files[tax_name], 'a') as out_file:
                            out_file.write('%s;0;%s;%s;%s;%f;%s\n' % (
                            reads_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

            os.remove(alt_reads_file)

    if reads_type == 'fastq':
        os.remove(reads_path_fasta)

    if tax_files is not None:
        os.remove(confusion_table)



if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers(help = 'task to be performed', dest = 'command')

    prepare_parser = subparsers.add_parser('prepare', help = 'prepare read data and references')
    prepare_parser.add_argument('data_set', help = 'data set to prepare (balanced, hmp)')
    prepare_parser.add_argument('out_dir', help = 'output directory')
    prepare_parser.add_argument('--max_n', type = int, help = 'maximum number of ambiguous bases')
    prepare_parser.add_argument('--max_ee', type = int, help = 'maximum number of "expected errors"')
    prepare_parser.add_argument('--trunc_q', type = int, help = 'truncate reads at the first instance of a quality score less than or equal')
    prepare_parser.add_argument('--trunc_len', help = 'truncate reads after specified number of bases and discard shorter ones (two comma-separated values for paired reads)')
    prepare_parser.add_argument('--trim_left', help = 'number of nucleotides to remove from the start of each read (two comma-separated values for paired reads)')
    prepare_parser.add_argument('--min_size', type = int, help = 'minimum cluster size (dereplication)')
    prepare_parser.add_argument('--fastq_min_ovlen', type = int, help = 'discard pair if alignment is shorter (only paired reads)')
    prepare_parser.add_argument('--fastq_max_diffs', type = int, help = 'maximum number of mismatches in the alignment (only paired reads)')
    prepare_parser.add_argument('--single', action = 'store_true', help = 'prepare only the forward reads')
    prepare_parser.add_argument('--clean', action = 'store_true', help = 'discard intermediate read files')
    prepare_parser.add_argument('--r_code', default = DEF_R_CODE_PATH, help = 'path to underlying R code')
    prepare_parser.add_argument('--usearch8', default = DEF_USEARCH8_PATH, help = 'path to USEACH v8.0 binary')
    prepare_parser.add_argument('--usearch10', default = DEF_USEARCH10_PATH, help = 'path to USEARCH v10 binary')

    ref_parser = subparsers.add_parser('reference', help = 'compute reference')
    ref_parser.add_argument('ref_name', help = 'reference to prepare (balanced, hmp, silva)')
    ref_parser.add_argument('rep_set', help = 'FASTA file of representatives')
    ref_parser.add_argument('reference_file', help = '(reduced and relabelled) reference sequences')
    ref_parser.add_argument('-t', '--taxonomy', help = 'taxonomy of references (only needed for silva)')
    ref_parser.add_argument('-l', '--level', type = int, default = rst.DEF_LEVEL, help = 'highest level of taxonomic information that has to be complete (only used by silva)')

    tax_parser = subparsers.add_parser('taxonomy', help = 'compute taxonomy files')
    tax_parser.add_argument('reads_file', help = 'reads to which taxonomic information is to be assigned')
    tax_parser.add_argument('reference_file', help = 'reference sequences with taxonomic information')
    tax_parser.add_argument('tax_file', help = 'output taxonomy file in blast6out format')
    tax_parser.add_argument('gt_threshold', type = float, help = 'ground truth threshold')
    tax_parser.add_argument('--reads_type', default = 'fastq', help = 'type of reads file (fasta, fastq)')
    tax_parser.add_argument('--vsearch', default = util.DEF_VSEARCH_PATH, help = 'path to VSEARCH binary')

    run_parser = subparsers.add_parser('run', help = 'execute analysis')
    run_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_parser.add_argument('reads_file', help = 'name of read collection and path to file, separated by a colon (e.g. file1:path/to/file)')
    run_parser.add_argument('tasks_file', help = 'list file of tasks')
    run_parser.add_argument('out_dir', help = 'output directory')
    run_parser.add_argument('--tax_files', help = 'comma-separated list of pairs consisting of name of and path to taxonomy files (components of each pair separated by colon; activates analysis of clustering quality)')
    run_parser.add_argument('--gefast', default = util.DEF_GEFAST_PATH, help = 'path to GeFaST binary')

    run_dada2_parser = subparsers.add_parser('run_dada2', help = 'execute analysis with DADA2')
    run_dada2_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_dada2_parser.add_argument('forward_reads', help = 'FASTQ file with (preprocessed) forward reads')
    run_dada2_parser.add_argument('out_dir', help = 'output directory')
    run_dada2_parser.add_argument('-r', '--reverse_reads', help = 'FASTQ file with (preprocessed) reverse reads')
    run_dada2_parser.add_argument('--tax_files', help = 'comma-separated list of pairs consisting of name of and path to taxonomy files (components of each pair separated by colon; activates analysis of clustering quality)')
    run_dada2_parser.add_argument('--r_code', default = DEF_R_CODE_PATH, help = 'path to underlying R code')

    run_uv_parser = subparsers.add_parser('run_uvsearch', help = 'execute analysis with USEARCH and VSEARCH')
    run_uv_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_uv_parser.add_argument('reads_file', help = 'name of read collection and path to file, separated by a colon (e.g. file1:path/to/file)')
    run_uv_parser.add_argument('tasks_file', help = 'list file of tasks')
    run_uv_parser.add_argument('out_dir', help = 'output directory')
    run_uv_parser.add_argument('--tax_files', help = 'comma-separated list of pairs consisting of name of and path to taxonomy files (components of each pair separated by colon; activates analysis of clustering quality)')
    run_uv_parser.add_argument('--usearch', default = util.DEF_USEARCH_PATH, help = 'path to USEARCH binary')
    run_uv_parser.add_argument('--vsearch', default = util.DEF_VSEARCH_PATH, help = 'path to VSEARCH binary')


    args = argparser.parse_args()


    if args.command == 'prepare':
        paths = dict()
        paths['r_code'] = args.r_code
        paths['usearch8'] = args.usearch8
        paths['usearch10'] = args.usearch10

        if args.data_set == 'balanced':
            prepare_balanced(args.out_dir, paths, single = args.single, max_n = args.max_n, max_ee = args.max_ee,
                             trunc_q = args.trunc_q, trunc_len = args.trunc_len, trim_left = args.trim_left,
                             min_size = args.min_size, fastq_min_ovlen = args.fastq_min_ovlen,
                             fastq_max_diffs = args.fastq_max_diffs, clean = args.clean)

        elif args.data_set == 'hmp':
            prepare_hmp(args.out_dir, paths, single = args.single, max_n = args.max_n, max_ee = args.max_ee,
                        trunc_q = args.trunc_q, trunc_len = args.trunc_len, trim_left = args.trim_left,
                        min_size = args.min_size, fastq_min_ovlen = args.fastq_min_ovlen,
                        fastq_max_diffs = args.fastq_max_diffs, clean = args.clean)

        else:
            raise ValueError('ERROR: "data_set" has to be "balanced" or "hmp", but was "%s".' % args.data_set)


    if args.command == 'reference':
        if args.ref_name == 'balanced':
            print('Nothing to do. Provided reference data set for "%s" can be used directly.' % args.ref_name)
            sp.call('cp %s %s' % (args.rep_set, args.reference_file), shell = True)

        elif args.ref_name == 'hmp':
            # sequences belonging together are numbered consecutively
            with open(args.rep_set, 'r') as in_file, open(args.reference_file, 'w') as out_file:
                for record in SeqIO.parse(in_file, 'fasta'):
                    id = '.'.join(record.id.split('.')[:-1])
                    out_file.write('>%s\n%s\n' % (id, record.seq))

        elif args.ref_name == 'silva':
            if args.taxonomy is None:
                raise ValueError('ERROR: reference "silva" needs both a set of representatives and taxonomy information.')
            at.reduce_and_label(args.rep_set, args.taxonomy, args.reference_file, level = args.level)

        else:
            raise ValueError('ERROR: "ref_name" has to be "balanced", "hmp" or "silva", but was "%s".' % args.data_set)


    if args.command == 'taxonomy':
        fasta_reads = args.reads_file
        if args.reads_type == 'fastq':
            fasta_reads = '%s_tmp_fasta' % args.tax_file
            with open(args.reads_file, 'r') as in_file, open(fasta_reads, 'w') as out_file:
                for record in SeqIO.parse(in_file, 'fastq'):
                    out_file.write('>%s\n%s\n' % (record.id, str(record.seq)))

        at.match_against_references(fasta_reads, args.reference_file, args.gt_threshold, args.tax_file, args.vsearch)

        if args.reads_type == 'fastq':
            os.remove(fasta_reads)


    if args.command == 'run':
        reads = util.parse_reads_list(args.reads_file)[0]
        tasks = util.parse_tasks_list_file(args.tasks_file)
        tax_files = [t.split(':') for t in args.tax_files.split(',')] if args.tax_files is not None else None

        binary_paths = dict()
        binary_paths['gefast'] = args.gefast

        run_analysis(args.analysis_name, reads, tasks, args.out_dir, binary_paths, tax_files = tax_files)


    if args.command == 'run_dada2':
        tax_files = [t.split(':') for t in args.tax_files.split(',')] if args.tax_files is not None else None

        paths = dict()
        paths['r_code'] = args.r_code

        run_dada2_analysis(args.analysis_name, args.forward_reads, args.reverse_reads, args.out_dir, paths, tax_files)


    if args.command == 'run_uvsearch':
        reads = util.parse_reads_list(args.reads_file)[0]
        tasks = util.parse_tasks_list_file(args.tasks_file)
        tax_files = [t.split(':') for t in args.tax_files.split(',')] if args.tax_files is not None else None

        binary_paths = dict()
        binary_paths['usearch'] = args.usearch
        binary_paths['vsearch'] = args.vsearch

        run_uvsearch_analysis(args.analysis_name, reads, tasks, args.out_dir, binary_paths, tax_files = tax_files)
