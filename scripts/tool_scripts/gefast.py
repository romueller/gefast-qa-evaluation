#!/usr/bin/python

# Organises the measured (runtime & memory) execution of GeFaST with different parameters / options.

import subprocess as sp

from .. import configuration_handling as ch
from .. import util

DEF_LOG_CMD = ''


# Run GeFaST after obtaining further arguments from a string (and log runtime and memory consumption).
def run_gefast_str_args(binary, mode, in_file, config, threshold, out_file,
                        args_string, args_sep = util.MISC_ARGS_SEP, log_cmd = DEF_LOG_CMD, refinement_threshold = None):
    op_c = False  # clean
    op_m = False  # mothur
    op_l = False  # list_file
    op_tw = False  # use_score
    op_us = False  # use_score
    op_uq = False  # use_qgrams

    if 'clean' in args_string:
        op_c = True

    op_alphabet = util.get_value(args_string, 'alphabet', args_sep)
    op_out_i = util.get_value(args_string, 'output_internal', args_sep)
    op_out_s = util.get_value(args_string, 'output_statistics', args_sep)
    op_out_w = util.get_value(args_string, 'output_seeds', args_sep)
    op_out_u = util.get_value(args_string, 'output_uclust', args_sep)
    op_n = util.get_value(args_string, 'break_swarms', args_sep)
    op_boundary = util.get_value(args_string, 'boundary', args_sep)
    if 'mothur' in args_string:
        op_m = True
    op_min_len = util.get_value(args_string, 'min_length', args_sep)
    op_max_len = util.get_value(args_string, 'max_length', args_sep)
    op_min_abund = util.get_value(args_string, 'min_abundance', args_sep)
    op_max_abund = util.get_value(args_string, 'max_abundance', args_sep)
    op_sep = util.get_value(args_string, 'sep_abundance', args_sep)
    op_preprocessor = util.get_value(args_string, 'preprocessor', args_sep)
    op_clusterer = util.get_value(args_string, 'clusterer', args_sep)
    op_refiner = util.get_value(args_string, 'refiner', args_sep)
    op_out_gen = util.get_value(args_string, 'output_generator', args_sep)
    op_qual_enc = util.get_value(args_string, 'quality_encoding', args_sep)
    op_match = util.get_value(args_string, 'match_reward', args_sep)
    op_mismatch = util.get_value(args_string, 'mismatch_penalty', args_sep)
    op_gap_open = util.get_value(args_string, 'gap_opening_penalty', args_sep)
    op_gap_extend = util.get_value(args_string, 'gap_extension_penalty', args_sep)
    op_keep = util.get_value(args_string, 'keep_preprocessed', args_sep)
    op_misc = util.get_value(args_string, 'misc', args_sep)
    if 'list_file' in args_string:
        op_l = True

    # mode-specific parameters
    op_num_segs = util.get_value(args_string, 'num_extra_segments', args_sep)
    op_bands = util.get_value(args_string, 'bands_per_side', args_sep)
    op_distance = util.get_value(args_string, 'distance', args_sep)
    if 'two_way_segment_filter' in args_string:
        op_tw = True
    if 'use_score' in args_string:
        op_us = True
    if 'use_qgrams' in args_string:
        op_uq = True
    op_amplicon_storage = util.get_value(args_string, 'amplicon_storage', args_sep)

    op_representation = util.get_value(args_string, 'representation', args_sep)
    op_ampl_coll = util.get_value(args_string, 'amplicon_collection', args_sep)
    op_auxiliary = util.get_value(args_string, 'auxiliary', args_sep)

    run_gefast(binary, mode, in_file, config, threshold, out_file, refinement_threshold = refinement_threshold,
               clean = op_c, log_cmd = log_cmd, alphabet = op_alphabet, out_internal = op_out_i, out_statistics = op_out_s,
               out_seeds = op_out_w, out_uclust = op_out_u, break_swarms = op_n, boundary = op_boundary, mothur = op_m,
               min_length = op_min_len, max_length = op_max_len, min_abund = op_min_abund, max_abund = op_max_abund,
               separator = op_sep, preprocessor = op_preprocessor, clusterer = op_clusterer, refiner = op_refiner,
               output_generator = op_out_gen, quality_encoding = op_qual_enc, match_reward = op_match,
               mismatch_penalty = op_mismatch, gap_open_penalty = op_gap_open, gap_extend_penalty = op_gap_extend,
               keep_preprocessed = op_keep, misc = op_misc, list_file = op_l, num_extra_segments = op_num_segs, bands_per_side = op_bands,
               distance = op_distance, two_way_segment_filter = op_tw, use_score = op_us, use_qgrams = op_uq,
               amplicon_storage = op_amplicon_storage, representation = op_representation, amplicon_collection = op_ampl_coll,
               auxiliary = op_auxiliary)


# Run GeFaST and log time / memory.
def run_gefast(binary, mode, in_file, config, threshold, out_file, refinement_threshold = None, clean = False,
               log_cmd = DEF_LOG_CMD, alphabet = None, out_internal = None, out_statistics = None, out_seeds = None,
               out_uclust = None, break_swarms = None, boundary = None, mothur = False, min_length = None, max_length = None,
               min_abund = None, max_abund = None, separator = None, preprocessor = None, clusterer = None, refiner = None,
               output_generator = None, quality_encoding = None, match_reward = None, mismatch_penalty = None, gap_open_penalty = None,
               gap_extend_penalty = None, keep_preprocessed = None, misc = None, list_file = False, num_extra_segments = None, bands_per_side = None, distance = None,
               two_way_segment_filter = False, use_score = False, use_qgrams = False, amplicon_storage = None, representation = None, amplicon_collection = None,
               auxiliary = None):
    tmp_config = '%s_tmp_config' % out_file
    if config != '':
        ch.copy_config(config, tmp_config)
    else:
        ch.create_config(tmp_config)
    ch.add_to_config(tmp_config, 'threshold', threshold)
    ch.add_to_config(tmp_config, 'output_otus', out_file)
    if refinement_threshold is not None:
        ch.add_to_config(tmp_config, 'refinement_threshold', refinement_threshold)

    if alphabet is not None:
        ch.add_to_config(tmp_config, 'alphabet', alphabet)
    if out_internal is not None:
        ch.add_to_config(tmp_config, 'output_internal', out_internal)
    if out_statistics is not None:
        ch.add_to_config(tmp_config, 'output_statistics', out_statistics)
    if out_seeds is not None:
        ch.add_to_config(tmp_config, 'output_seeds', out_seeds)
    if out_uclust is not None:
        ch.add_to_config(tmp_config, 'output_uclust', out_uclust)
    if break_swarms is not None:
        ch.add_to_config(tmp_config, 'break_swarms', break_swarms)
    if boundary is not None:
        ch.add_to_config(tmp_config, 'boundary', boundary)
    if mothur:
        ch.add_to_config(tmp_config, 'mothur', '1')
    if min_length is not None:
        ch.add_to_config(tmp_config, 'min_length', min_length)
    if max_length is not None:
        ch.add_to_config(tmp_config, 'max_length', max_length)
    if min_abund is not None:
        ch.add_to_config(tmp_config, 'min_abundance', min_abund)
    if max_abund is not None:
        ch.add_to_config(tmp_config, 'max_abundance', max_abund)
    if separator is not None:
        ch.add_to_config(tmp_config, 'sep_abundance', separator)
    if preprocessor is not None:
        ch.add_to_config(tmp_config, 'preprocessor', preprocessor)
    if clusterer is not None:
        ch.add_to_config(tmp_config, 'clusterer', clusterer)
    if refiner is not None:
        ch.add_to_config(tmp_config, 'refiner', refiner)
    if output_generator is not None:
        ch.add_to_config(tmp_config, 'output_generator', output_generator)
    if quality_encoding is not None:
        ch.add_to_config(tmp_config, 'quality_encoding', quality_encoding)
    if match_reward is not None:
        ch.add_to_config(tmp_config, 'match_reward', match_reward)
    if mismatch_penalty is not None:
        ch.add_to_config(tmp_config, 'mismatch_penalty', mismatch_penalty)
    if gap_open_penalty is not None:
        ch.add_to_config(tmp_config, 'gap_opening_penalty', gap_open_penalty)
    if gap_extend_penalty is not None:
        ch.add_to_config(tmp_config, 'gap_extension_penalty', gap_extend_penalty)
    if keep_preprocessed is not None:
        ch.add_to_config(tmp_config, 'keep_preprocessed', keep_preprocessed)
    if misc is not None:
        ch.add_to_config(tmp_config, 'misc', misc)

    op_list_file = '--list_file' if list_file else ''

    if num_extra_segments is not None:
        ch.add_to_config(tmp_config, 'num_extra_segments', num_extra_segments)
    if bands_per_side is not None:
        ch.add_to_config(tmp_config, 'bands_per_side', bands_per_side)
    if distance is not None:
        ch.add_to_config(tmp_config, 'distance', distance)
    if two_way_segment_filter:
        ch.add_to_config(tmp_config, 'two_way_segment_filter', '1')
    if use_score:
        ch.add_to_config(tmp_config, 'use_score', '1')
    if use_qgrams:
        ch.add_to_config(tmp_config, 'use_qgrams', '1')
    if amplicon_storage is not None:
        ch.add_to_config(tmp_config, 'amplicon_storage', amplicon_storage)
    if representation is not None:
        ch.add_to_config(tmp_config, 'representation', representation)
    if amplicon_collection is not None:
        ch.add_to_config(tmp_config, 'amplicon_collection', amplicon_collection)
    if auxiliary is not None:
        ch.add_to_config(tmp_config, 'auxiliary', auxiliary)

    sp.call('%s %s %s %s %s %s' % (log_cmd, binary, mode, in_file, tmp_config, op_list_file), shell = True)

    ch.remove_file(tmp_config)

    if clean:
        sp.call('rm %s' % out_file, shell = True)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('mode', help = 'GeFaST mode')
    argparser.add_argument('in_file', help = 'input FASTA file')
    argparser.add_argument('config', help = 'configuration file')
    argparser.add_argument('threshold', type = int, help = 'clustering threshold')
    argparser.add_argument('out_file', help = 'output file')
    argparser.add_argument('-r', '--refinement_threshold', help = 'refinement threshold(s)')

    argparser.add_argument('-a', '--alphabet', help = 'eponymous GeFaST option')
    argparser.add_argument('-i', '--output_internal', help = 'eponymous GeFaST option')
    argparser.add_argument('-s', '--output_statistics', help = 'eponymous GeFaST option')
    argparser.add_argument('-w', '--output_seeds', help = 'eponymous GeFaST option')
    argparser.add_argument('-u', '--output_uclust', help = 'eponymous GeFaST option')
    argparser.add_argument('-n', '--break_swarms', help = 'eponymous GeFaST option')
    argparser.add_argument('-b', '--boundary', help = 'eponymous GeFaST option')
    argparser.add_argument('--mothur', action = 'store_true', help = 'eponymous GeFaST option')
    argparser.add_argument('--min_length', help = 'eponymous GeFaST option')
    argparser.add_argument('--max_length', help = 'eponymous GeFaST option')
    argparser.add_argument('--min_abundance', help = 'eponymous GeFaST option')
    argparser.add_argument('--max_abundance', help = 'eponymous GeFaST option')
    argparser.add_argument('--sep_abundance', help = 'eponymous GeFaST option')
    argparser.add_argument('--preprocessor', help = 'eponymous GeFaST option')
    argparser.add_argument('--clusterer', help = 'eponymous GeFaST option')
    argparser.add_argument('--refiner', help = 'eponymous GeFaST option')
    argparser.add_argument('--output_generator', help = 'eponymous GeFaST option')
    argparser.add_argument('--quality_encoding', help = 'eponymous GeFaST option')
    argparser.add_argument('-m', '--match_reward', help = 'eponymous GeFaST option')
    argparser.add_argument('-p', '--mismatch_penalty', help = 'eponymous GeFaST option')
    argparser.add_argument('-g', '--gap_opening_penalty', help = 'eponymous GeFaST option')
    argparser.add_argument('-e', '--gap_extension_penalty', help = 'eponymous GeFaST option')
    argparser.add_argument('-k', '--keep_preprocessed', help = 'eponymous GeFaST option')
    argparser.add_argument('--list_file', action = 'store_true', help = 'eponymous GeFaST option')
    argparser.add_argument('--num_extra_segments', help = 'GeFaST option: number of extra segments for segment filter')
    argparser.add_argument('--bands_per_side', help = 'GeFaST option: number of bands per side in DP-matrix of alignment computation')
    argparser.add_argument('--distance', help = 'GeFaST option: used function for distance calculation')
    argparser.add_argument('--two_way_segment_filter', action = 'store_true', help = 'GeFaST option: use two-way segment filter')
    argparser.add_argument('--use_score', action = 'store_true', help = 'GeFaST option: use scoring function')
    argparser.add_argument('--use_qgrams', action = 'store_true', help = 'GeFaST option: use q-gram filter')
    argparser.add_argument('--amplicon_storage', help = 'GeFaST option: select data structure for amplicon storage')
    argparser.add_argument('--representation', help = 'GeFaST option: select data structure for feature representation')
    argparser.add_argument('--amplicon_collection', help = 'GeFaST option: select data structure for amplicon collection')
    argparser.add_argument('--auxiliary', help = 'GeFaST option: select data structure for auxiliary data')
    argparser.add_argument('--misc', help = 'GeFaST option: miscellaneous options')

    argparser.add_argument('-c', '--clean', action = 'store_true', help = 'remove output file')
    argparser.add_argument('--gefast', default = util.DEF_GEFAST_PATH, help = 'path to GeFaST binary')
    argparser.add_argument('--log_cmd', default = DEF_LOG_CMD, help = 'log command')
    args = argparser.parse_args()

    run_gefast(args.gefast, args.mode, args.in_file, args.config, args.threshold, args.out_file,
               refinement_threshold = args.refinement_threshold, clean = args.clean, log_cmd = args.log_cmd,
               alphabet = args.alphabet, out_internal = args.output_internal, out_statistics = args.output_statistics,
               out_seeds = args.output_seeds, out_uclust = args.output_uclust, break_swarms = args.break_swarms,
               boundary = args.boundary, mothur = args.mothur, min_length = args.min_length, max_length = args.max_length,
               min_abund = args.min_abundance, max_abund = args.max_abundance, separator = args.sep_abundance,
               preprocessor = args.preprocessor, clusterer = args.clusterer, refiner = args.refiner,
               output_generator = args.output_generator, quality_encoding = args.quality_encoding,
               match_reward = args.match_reward, mismatch_penalty = args.mismatch_penalty,
               gap_open_penalty = args.gap_opening_penalty, gap_extend_penalty = args.gap_extension_penalty,
               keep_preprocessed = args.keep_preprocessed, misc = args.misc, list_file = args.list_file,
               num_extra_segments = args.num_extra_segments, bands_per_side = args.bands_per_side,
               distance = args.distance, two_way_segment_filter = args.two_way_segment_filter, use_score = args.use_score,
               use_qgrams = args.use_qgrams, amplicon_storage = args.amplicon_storage, representation = args.representation,
               amplicon_collection = args.amplicon_collection, auxiliary = args.auxiliary)
