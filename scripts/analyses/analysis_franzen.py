#!/usr/bin/python

# Prepares and runs a clustering-quality analysis on simulated amplicon reads
# similar to the analysis in the following paper:
#
# Franzén O, Hu J, Bao X, Itzkowitz SH, Peter I, Bashir A. (2015)
# Improved OTU-picking using long-read 16S rRNA gene amplicon sequencing and generic hierarchical clustering.
# https://doi.org/10.1186/s40168-015-0105-6
#
# Analyses the precision, recall and adjusted Rand index of the computed clusters.
#
# Offers three kinds of analysis runs:
# (1) GeFaST with its diverse modes to evaluate the impact of quality-weighted alignments.
# (2) DADA2 (used as a clustering tool, for comparison with DADA2-inspired GeFaST clustering / refinement approaches)
# (3) USEARCH and VSEARCH (for comparison with DADA2-inspired GeFaST clustering / refinement approaches)

import os.path
import pandas as pd
import random
import subprocess as sp

from Bio import AlignIO, SeqIO

from .. import util
from ..cluster_evaluation import compute_confusion_table as cct
from ..cluster_evaluation import evaluate_confusion_table as ect
from ..data_preparation import create_alternative_inputs as cai
from ..data_preparation import fastq2fasta as f2f
from ..tool_scripts import gefast as ex_gefast
from ..tool_scripts import usearch as ex_usearch
from ..tool_scripts import vsearch as ex_vsearch
from ..tool_scripts import uparse as ex_uparse
from ..tool_scripts import swarm as ex_swarm


DEF_R_CODE_PATH = 'scripts/tool_scripts/dada2.R'


# Download / extract / build the necessary data files.
# Format also the 16S rRNA MSA using Infernal (-> <out_dir>/bacteria16S_508_mod5.cmfile).
def prepare_data(out_dir, cmbuild_path):
    sp.call('mkdir -p %s' % out_dir, shell = True)

    # GreenGenes database
    sp.call('cd %s; ' % out_dir +
            'wget https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5.fasta.gz; ' +
            'gunzip -d gg_13_5.fasta.gz', shell = True)

    # hand-curated bacterial 16S rRNA MSA
    sp.call('cd %s; ' % out_dir +
            'wget http://rdp.cme.msu.edu/download/RDPinfernalTraindata.zip; ' +
            'unzip RDPinfernalTraindata.zip; ' +
            'rm -r RDPinfernalTraindata.zip __MACOSX/', shell = True)
    sp.call('%s --ere 1.4 %s/bacteria16S_508_mod5.cmfile %s/RDPinfernalTraindata/bacteria16S_508_mod5.stk' % (cmbuild_path, out_dir, out_dir),
            shell = True)



# Create a separate list file (of GreenGenes identifiers) per complexity and mock community.
def get_list_files(csv_file, out_dir):
    sp.call('mkdir -p %s' % out_dir, shell = True)
    df = pd.read_csv(csv_file, sep = ';')
    with open('%s/list.txt' % out_dir, 'w') as out_list:
        for (c, m), grp in df.groupby(['Complexity', 'Mock community']):
            with open('%s/%s_%i.txt' % (out_dir, c, m), 'w') as out_file:
                for acc in grp['GreenGenes identifier']:
                    out_file.write('%s\n' % acc)
            out_list.write('%s/%s_%i.txt\n' % (out_dir, c, m))


# Extract the FASTA entries of the specified sequences.
def get_references(gg_file, acc_file, ref_file, skip_ambiguous):
    accs = set()
    with open(acc_file, 'r') as in_file:
        for line in in_file:
            accs.add(line.strip())

    found = set()
    with open(gg_file, 'r') as in_file, open(ref_file, 'w') as out_file:
        for record in SeqIO.parse(in_file, 'fasta'):
            if record.id in accs:
                seq = str(record.seq)
                if skip_ambiguous:
                    seq = ''.join([c for c in seq if c in 'ACGTU'])
                out_file.write('>%s\n%s\n' % (record.id, seq))
                found.add(record.id)

    print('Sequences not found for mock community in %s (%i):\n%s' % (acc_file, len(accs - found), ', '.join(accs - found)))


# Extract amplicons from MSA in Stockholm format.
def extract(sto_file, first_pos, last_pos, fasta_file):
    msa = AlignIO.read(sto_file, 'stockholm')
    with open(fasta_file, 'w') as out_file:
        for i, record in enumerate(msa):
            gapped_subseq = str(record.seq)[first_pos - 1 : last_pos]
            out_file.write('>%s\n%s\n' % (record.id, gapped_subseq.replace('-', '').upper()))


# Simulate MiSeq reads (FASTA, FASTQ or both) using ART (art_illumina).
# The produced reads are a concatenation of the forward reads and the reverse-complement of the reverse reads,
# separated by 10 Ns (having the same quality as the last nucleotide of the forward read).
def simulate_miseq_reads(seq_file, out_prefix, length, coverage, type, num_n, art_path):
    sp.call('%s -amp -na -p -l %i -f %i -i %s -o %s' % (art_path, length, coverage, seq_file, out_prefix), shell = True)

    forward_reads = '%s1.fq' % out_prefix
    reverse_reads = '%s2.fq' % out_prefix
    seqs = dict()
    quals = dict()
    with open(forward_reads, 'r') as in_file:
        for record in SeqIO.parse(in_file, 'fastq'):
            seqs[record.id[:-2]] = '%s%s' % (str(record.seq), ''.join(['N'] * num_n))
            scores = SeqIO.QualityIO._get_sanger_quality_str(record)
            quals[record.id[:-2]] = '%s%s' % (scores, ''.join([str(scores[-1])] * num_n))
    with open(reverse_reads, 'r') as in_file:
        for record in SeqIO.parse(in_file, 'fastq'):
            seqs[record.id[:-2]] += str(record.seq.reverse_complement())
            quals[record.id[:-2]] += SeqIO.QualityIO._get_sanger_quality_str(record)[::-1]

    sp.call('rm %s %s' % (forward_reads, reverse_reads), shell = True)

    if type in ['fasta', 'both']:
        synthetic_reads = '%s_miseq.fasta' % out_prefix
        with open(synthetic_reads, 'w') as out_file:
            for s in seqs:
                out_file.write('>%s\n%s\n' % (s, seqs[s]))

    if type in ['fastq', 'both']:
        synthetic_reads = '%s_miseq.fastq' % out_prefix
        with open(synthetic_reads, 'w') as out_file:
            for s in seqs:
                out_file.write('@%s\n%s\n+\n%s\n' % (s, seqs[s], quals[s]))


# Simulate PacBio reads (FASTA, FASTQ or both) using pbsim.
def simulate_pacbio_reads(seq_file, out_prefix, model_qc, acc_mean, acc_sd, diff_ratio, data_type, depth, type, pbsim_path):
    sp.call('%s --accuracy-mean %f --accuracy-sd %f --difference-ratio %s --data-type %s --model_qc %s --depth %i --prefix %s %s'
            % (pbsim_path, acc_mean, acc_sd, diff_ratio, data_type, model_qc, depth, out_prefix, seq_file), shell = True)

    seq_ids = []
    with open(seq_file, 'r') as in_file:
        for record in SeqIO.parse(in_file, 'fasta'):
            seq_ids.append(record.id)

    synthetic_fasta = '%s_pacbio.fasta' % out_prefix
    synthetic_fastq = '%s_pacbio.fastq' % out_prefix

    with open(synthetic_fasta, 'w') as out_fasta, open(synthetic_fastq, 'w') as out_fastq:
        for i in range(1, len(seq_ids) + 1):
            file_prefix = '%s_%04d' % (out_prefix, i)
            with open('%s.fastq' % file_prefix, 'r') as in_file:
                # not possible because pbsim adds an additional quality value for every sequence (except the first and last one) -> corrupted FASTQ file that cannot be parsed
                #for record in SeqIO.parse(in_file, 'fastq'):
                #    out_file.write('>%s.%s\n%s\n' % (seq_ids[i], record.id.split('_')[1], str(record.seq)))
                lines = [l.strip() for l in in_file.readlines() if l != '\n']

                # pbsim FASTQ files sometimes have unexpected newline in nucleotide sequence
                k = 0
                j = 1
                while k < len(lines):
                    k += 1
                    seq = ''
                    lab = '+S%i_%i' % (i, j)
                    while not lines[k].startswith(lab):
                        seq += lines[k]
                        k += 1

                    k += 1
                    qual = ''
                    lab = '@S%i_%i' % (i, j + 1)
                    while k < len(lines) and not lines[k].startswith(lab):
                        qual += lines[k]
                        k += 1

                    # cap the quality scores at J
                    # (based on discussion on http://seqanswers.com/forums/showthread.php?t=48036)
                    qual = ''.join([(c if c <= 'J' else 'J') for c in qual])

                    # use hyphen instead of dot to separate id and repetition number (like MiSeq files)
                    if type in ['fasta', 'both']:
                        out_fasta.write('>%s-%s\n%s\n' % (seq_ids[i - 1], j, seq))
                    if type in ['fastq', 'both']:
                        out_fastq.write('@%s-%s\n%s\n+\n%s\n' % (seq_ids[i - 1], j, seq, qual))

                    j += 1

            sp.call('rm %s.*' % file_prefix, shell = True)

    if type != 'both':
        sp.call('rm %s' % (synthetic_fasta if type == 'fastq' else synthetic_fastq), shell = True)


# Subsample the reads (per genome), shuffle the resulting reads and output the sequences in FASTA format.
def create_shuffled_subset_fasta(seq_file, num_chosen, sub_file, rep_sep):
    seqs = dict()
    with open(seq_file, 'r') as in_file:
        for record in SeqIO.parse(in_file, 'fasta'):
            sid, rep = record.id.split(rep_sep)
            if sid not in seqs:
                seqs[sid] = []
            seqs[sid].append((rep, str(record.seq)))

    entries = []
    for s in seqs:
        indices = list(range(0, len(seqs[s])))
        random.shuffle(indices)
        choice = indices[:num_chosen]
        for c in choice:
            rep, seq = seqs[s][c]
            entries.append('>%s.%s\n%s\n' % (s, rep, seq))
    random.shuffle(entries)

    with open(sub_file, 'w') as out_file:
        for e in entries:
            out_file.write(e)


# Subsample the reads (per genome), shuffle the resulting reads and output the sequences in FASTQ (and FASTA) format.
# The sequences have the same order in the FASTQ and the FASTA file.
def create_shuffled_subset_fastq(seq_file, num_chosen, sub_file, rep_sep, both = False):
    seqs = dict()
    with open(seq_file, 'r') as in_file:
        for record in SeqIO.parse(in_file, 'fastq'):
            sid, _ = record.id.split(rep_sep)
            if sid not in seqs:
                seqs[sid] = []
            seqs[sid].append(record)

    entries = []
    for s in seqs:
        indices = list(range(0, len(seqs[s])))
        random.shuffle(indices)
        choice = indices[:num_chosen]
        for c in choice:
            entries.append(seqs[s][c])
    random.shuffle(entries)

    with open(sub_file, 'w') as out_file:
        SeqIO.write(entries, out_file, 'fastq')

    if both:
        fasta_file = '%s/%s.fasta' % (os.path.dirname(sub_file), os.path.splitext(os.path.basename(sub_file))[0])
        with open(fasta_file, 'w') as out_file:
            SeqIO.write(entries, out_file, 'fasta-2line')


# Create a taxonomy file in blast6out format.
#
# As only the first two columns (query and target labels) are needed for the later evaluation,
# the remaining columns are filled with dummy values.
def build_taxonomy_file(seq_file, tax_file, type):
    with open(seq_file, 'r') as in_file, open(tax_file, 'w') as out_file:
        for record in SeqIO.parse(in_file, type):
            ref_id = record.id.split('-')[0]
            out_file.write('%s\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n' % (record.id, ref_id))


# For each mock community, extract the reference sequences from the GreenGenes database
# and align them to the 16S rRNA MSA model.
# Subsequently, extract 4 sets of amplicons from the MSA and for each set of amplicons simulate
# MiSeq and / or PacBio reads.
def get_mock_community(name, gg_file, acc_file, msa_model, out_dir, cmalign_path, art_path, pbsim_path,
                       skip_ambiguous = True, num_n = 0, get_miseq = False, model_qc = None):

    sp.call('mkdir -p %s' % out_dir, shell = True)

    # extract reference sequences
    ref_file = '%s/%s_references.fasta' % (out_dir, name)
    get_references(gg_file, acc_file, ref_file, skip_ambiguous)

    # align sequences to 16S rRNA MSA model
    msa_file = '%s/%s_references.msa' % (out_dir, name)
    sp.call('%s %s %s > %s' % (cmalign_path, msa_model, ref_file, msa_file), shell = True)

    # extract 4 sets of amplicons covering different parts of the references
    # and, for each set of amplicons, simulate MiSeq and / or PacBio reads

    # interval 1 (MiSeq 2 x 150)
    int_name = 'V4'
    first_pos = 389
    last_pos = 801
    out_prefix = '%s/%s_%s' % (out_dir, name, int_name)
    amplicons_file = '%s_amplicons.fasta' % out_prefix
    if get_miseq:
        extract(msa_file, first_pos, last_pos, amplicons_file)

    if get_miseq:
        simulate_miseq_reads(amplicons_file, out_prefix, 150, 100, 'fastq', num_n, art_path)
        reads_file = '%s_miseq.fastq' % out_prefix
        reduced_shuffled_file = '%s_miseq_rs.fastq' % out_prefix
        tax_file = '%s_miseq_rs.tax' % out_prefix
        create_shuffled_subset_fastq(reads_file, 20, reduced_shuffled_file, '-', both = False)
        build_taxonomy_file(reduced_shuffled_file, tax_file, 'fastq')

    # interval 2 (MiSeq 2 x 250, PacBio 450)
    int_name = 'V3-V4'
    first_pos = 227
    last_pos = 801
    out_prefix = '%s/%s_%s' % (out_dir, name, int_name)
    amplicons_file = '%s_amplicons.fasta' % out_prefix
    if get_miseq or (model_qc is not None):
        extract(msa_file, first_pos, last_pos, amplicons_file)

    if get_miseq:
        simulate_miseq_reads(amplicons_file, out_prefix, 250, 100, 'fastq', num_n, art_path)
        reads_file = '%s_miseq.fastq' % out_prefix
        reduced_shuffled_file = '%s_miseq_rs.fastq' % out_prefix
        tax_file = '%s_miseq_rs.tax' % out_prefix
        create_shuffled_subset_fastq(reads_file, 20, reduced_shuffled_file, '-', both = False)
        build_taxonomy_file(reduced_shuffled_file, tax_file, 'fastq')

    if model_qc is not None:
        simulate_pacbio_reads(amplicons_file, out_prefix, model_qc, 0.99, 0.01, '6:21:73', 'CLR', 50, 'fastq', pbsim_path)
        reads_file = '%s_pacbio.fastq' % out_prefix
        reduced_shuffled_file = '%s_pacbio_rs.fastq' % out_prefix
        tax_file = '%s_pacbio_rs.tax' % out_prefix
        create_shuffled_subset_fastq(reads_file, 20, reduced_shuffled_file, '-', both = False)
        build_taxonomy_file(reduced_shuffled_file, tax_file, 'fastq')

    # interval 3 (PacBio 750)
    int_name = 'V1-V4'
    first_pos = 4
    last_pos = 801
    out_prefix = '%s/%s_%s' % (out_dir, name, int_name)
    amplicons_file = '%s_amplicons.fasta' % out_prefix
    if model_qc is not None:
        extract(msa_file, first_pos, last_pos, amplicons_file)

    if model_qc is not None:
        simulate_pacbio_reads(amplicons_file, out_prefix, model_qc, 0.99, 0.01, '6:21:73', 'CLR', 50, 'fastq', pbsim_path)
        reads_file = '%s_pacbio.fastq' % out_prefix
        reduced_shuffled_file = '%s_pacbio_rs.fastq' % out_prefix
        tax_file = '%s_pacbio_rs.tax' % out_prefix
        create_shuffled_subset_fastq(reads_file, 20, reduced_shuffled_file, '-', both = False)
        build_taxonomy_file(reduced_shuffled_file, tax_file, 'fastq')

    # interval 4 (PacBio 1450)
    int_name = 'V1-V6'
    first_pos = 4
    last_pos = 1506
    out_prefix = '%s/%s_%s' % (out_dir, name, int_name)
    amplicons_file = '%s_amplicons.fasta' % out_prefix
    if model_qc is not None:
        extract(msa_file, first_pos, last_pos, amplicons_file)

    if model_qc is not None:
        simulate_pacbio_reads(amplicons_file, out_prefix, model_qc, 0.99, 0.01, '6:21:73', 'CLR', 50, 'fastq', pbsim_path)
        reads_file = '%s_pacbio.fastq' % out_prefix
        reduced_shuffled_file = '%s_pacbio_rs.fastq' % out_prefix
        tax_file = '%s_pacbio_rs.tax' % out_prefix
        create_shuffled_subset_fastq(reads_file, 20, reduced_shuffled_file, '-', both = False)
        build_taxonomy_file(reduced_shuffled_file, tax_file, 'fastq')


# Compute clusters using GeFaST according to the tasks file and, if a taxonomy file is provided,
# analyse their quality.
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
def run_dada2_analysis(analysis_name, reads, out_dir, paths, tax_files = None):

    sp.call('mkdir -p %s' % out_dir, shell = True)

    out_prefix = '%s/%s' % (out_dir, analysis_name)
    metrics_files = dict()
    confusion_table = '%s/%s__conf_table' % (out_dir, analysis_name)
    if tax_files is not None:
        for tax_name, tax_file in tax_files:
            metrics_files[tax_name] = '%s/%s_%s__metrics.csv' % (out_dir, analysis_name, tax_name)
            with open(metrics_files[tax_name], 'w') as out_file:
                out_file.write('reads;gt;task;tool;mode;threshold;precision;recall;adjrandindex\n')

    sp.call('R -e \'source("%s"); dada2_cluster_single("%s", "%s")\''
            % (paths['r_code'], reads, out_prefix), shell = True)

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
            cai.create_file(binaries['vsearch'], reads_path_fasta, alt_reads_file, derep = True, sort_criterion = order,
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
            cai.create_file(binaries['vsearch'], reads_path_fasta, alt_reads_file, derep = True, sort_criterion = order,
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
            otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
            ex_uparse.run_uparse_str_args(binaries[tool], alt_reads_file, otus_file, misc,
                                         args_sep = util.MISC_ARGS_SEP)

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


# Compute clusters using Swarm and, if taxonomic assignments are provided, analyse their quality.
def run_swarm_analysis(analysis_name, reads, tasks, out_dir, binaries, tax_files = None):

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
            otus_file = '%s/%s_%f_otus.txt' % (out_dir, task, t)
            ex_swarm.run_swarm_str_args(binaries[tool], reads_path_swarm, t, otus_file, misc, args_sep = util.MISC_ARGS_SEP)

            if tax_files is not None:
                for tax_name, tax_file in tax_files:
                    cct.compute_table(tax_file, otus_file, confusion_table)
                    metric_values = ect.evaluate(confusion_table)
                    with open(metrics_files[tax_name], 'a') as out_file:
                        out_file.write('%s;0;%s;%s;%s;%f;%s\n' % (
                            reads_name, task, tool, mode, t, ';'.join(map(str, metric_values))))

        if derep:
            os.remove(alt_reads_file)

    if reads_type == 'fastq':
        os.remove(reads_path_fasta)

    if tax_files is not None:
        os.remove(confusion_table)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers(help = 'task to be performed', dest = 'command')

    prepare_parser = subparsers.add_parser('prepare', help = 'prepare tools and data')
    prepare_parser.add_argument('out_dir', help = 'output directory for list files')
    prepare_parser.add_argument('--cmbuild', help = 'path to cmbuild binary (part of Infernal)')

    lists_parser = subparsers.add_parser('lists', help = 'get mock-community list files')
    lists_parser.add_argument('csv_file', help = 'CSV file of genomes')
    lists_parser.add_argument('out_dir', help = 'output directory for list files')

    mock_parser = subparsers.add_parser('mock', help = 'simulate reads from mock community')
    mock_parser.add_argument('gg_file', help = 'Greengenes database FASTA file')
    mock_parser.add_argument('list_file', help = 'list file of mock communities')
    mock_parser.add_argument('msa_model', help = 'bacterial 16S rRNA MSA formatted with Infernal (cmbuild)')
    mock_parser.add_argument('out_dir', help = 'base output directory for samples')
    mock_parser.add_argument('--skip_ambiguous', action = 'store_true', help = 'discard ambiguous nucleotides (IUPAC code) while copying sequences from the Greengenes database')
    mock_parser.add_argument('--num_n', type = int, default = 0, help = 'number of Ns between forward and reverse portion in MiSeq reads')
    mock_parser.add_argument('--cmalign', help = 'path to cmalign binary (part of Infernal)')
    mock_parser.add_argument('--art', help = 'path to ART binary')
    mock_parser.add_argument('--pbsim', help = 'path to pbsim binary')
    mock_parser.add_argument('--miseq', action = 'store_true', help = 'create MiSeq-read mock communities')
    mock_parser.add_argument('--pacbio', help = 'model of quality code (pbsim; activates creating PacBio-read mock communities)')

    run_parser = subparsers.add_parser('run', help = 'execute analysis')
    run_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_parser.add_argument('reads_file', help = 'name of read collection and path to file, separated by a colon (e.g. file1:path/to/file)')
    run_parser.add_argument('tasks_file', help = 'list file of tasks')
    run_parser.add_argument('out_dir', help = 'output directory')
    run_parser.add_argument('--tax_files', help = 'comma-separated list of pairs consisting of name of and path to taxonomy files (components of each pair separated by colon; activates analysis of clustering quality)')
    run_parser.add_argument('--gefast', default = util.DEF_GEFAST_PATH, help = 'path to GeFaST binary')

    run_dada2_parser = subparsers.add_parser('run_dada2', help = 'execute analysis with DADA2')
    run_dada2_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_dada2_parser.add_argument('reads', help = 'FASTQ file with (preprocessed) reads')
    run_dada2_parser.add_argument('out_dir', help = 'output directory')
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

    run_swarm_parser = subparsers.add_parser('run_swarm', help = 'execute analysis with Swarm')
    run_swarm_parser.add_argument('analysis_name', help = 'name of analysis run')
    run_swarm_parser.add_argument('reads_file', help = 'name of read collection and path to file, separated by a colon (e.g. file1:path/to/file)')
    run_swarm_parser.add_argument('tasks_file', help = 'list file of tasks')
    run_swarm_parser.add_argument('out_dir', help = 'output directory')
    run_swarm_parser.add_argument('--tax_files', help = 'comma-separated list of pairs consisting of name of and path to taxonomy files (components of each pair separated by colon; activates analysis of clustering quality)')
    run_swarm_parser.add_argument('--swarm', default = util.DEF_SWARM_PATH, help = 'path to Swarm binary')

    args = argparser.parse_args()


    if args.command == 'prepare':
        prepare_data(args.out_dir, cmbuild_path = args.cmbuild)


    if args.command == 'lists':
        get_list_files(args.csv_file, args.out_dir)


    if args.command == 'mock':
        mock_files = []
        with open(args.list_file, 'r') as in_file:
            for line in in_file:
                mock_files.append(line.strip())

        for mf in mock_files:
            mock_name = os.path.splitext(os.path.basename(mf))[0]
            sample_dir = '%s/%s' % (args.out_dir, mock_name)

            print('Processing mock community %s...' % mock_name)
            get_mock_community(mock_name, args.gg_file, mf, args.msa_model, sample_dir, args.cmalign,
                               args.art, args.pbsim, skip_ambiguous = args.skip_ambiguous, num_n = args.num_n,
                               get_miseq = args.miseq, model_qc = args.pacbio)


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

        run_dada2_analysis(args.analysis_name, args.reads, args.out_dir, paths, tax_files)


    if args.command == 'run_uvsearch':
        reads = util.parse_reads_list(args.reads_file)[0]
        tasks = util.parse_tasks_list_file(args.tasks_file)
        tax_files = [t.split(':') for t in args.tax_files.split(',')] if args.tax_files is not None else None

        binary_paths = dict()
        binary_paths['usearch'] = args.usearch
        binary_paths['vsearch'] = args.vsearch

        run_uvsearch_analysis(args.analysis_name, reads, tasks, args.out_dir, binary_paths, tax_files = tax_files)


    if args.command == 'run_swarm':
        reads = util.parse_reads_list(args.reads_file)[0]
        tasks = util.parse_tasks_list_file(args.tasks_file)
        tax_files = [t.split(':') for t in args.tax_files.split(',')] if args.tax_files is not None else None

        binary_paths = dict()
        binary_paths['swarm'] = args.swarm

        run_swarm_analysis(args.analysis_name, reads, tasks, args.out_dir, binary_paths, tax_files = tax_files)
