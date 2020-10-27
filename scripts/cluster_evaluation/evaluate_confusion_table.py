#!/usr/bin/python

# Computes precision, recall and adjusted Rand index from a confusion table.
#
# Adapted from:
# Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014)
# Swarm: robust and fast clustering method for amplicon-based studies.
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
# Section: 2.6.4 Clustering metrics

import math
import numpy as np

DEF_TABLE_SEPARATOR = '\t'
DEF_OUTPUT_SEPARATOR = ','


# Read confusion table from file into a matrix.
def read_confusion_table(file, sep = DEF_TABLE_SEPARATOR):
    matrix = []
    with open(file, 'r') as in_file:
        next(in_file)  # skip header line
        for line in in_file:
            matrix.append(line.rstrip('\n').split(sep)[1:])

    return np.mat(matrix, 'float64')


# Compute precision, recall and adjusted Rand index from the confusion-table matrix.
def compute_metrics(conf_table):
    total = np.float64(np.sum(conf_table))

    if total > 0.0:
        recall = np.max(conf_table, axis = 0).sum() / total
        precision = np.max(conf_table, axis = 1).sum() / total

        row_totals = np.sum(conf_table, axis = 1)
        col_totals = np.sum(conf_table, axis = 0)

        row_combs = np.sum(np.multiply(row_totals, row_totals - 1) / 2)
        col_combs = np.sum(np.multiply(col_totals, col_totals - 1) / 2)
        combs = np.sum(np.multiply(conf_table, conf_table - 1) / 2)
        total_combs = total * (total - 1) / 2

        # check for pathological(?) cases (expected index = 1.0 -> index = 1.0 -> formula for corrected divides 0 by 0)
        # formula for expected index from Hubert & Arabie (1985), Comparing Partitions, Journal of Classification
        exp_index = 1 + 2 * row_combs * col_combs / (total_combs * total_combs) - (row_combs + col_combs) / total_combs
        if exp_index == 1.0:
            ari = 1.0
        else:
            tmp = (row_combs * col_combs) / total_combs
            ret = combs - tmp
            ari = ret / (0.5 * (row_combs + col_combs) - tmp)

    else:
        recall = 0.0
        precision = 0.0
        ari = 0.0

    return precision, recall, ari


# Run evaluation of confusion table and return metrics.
def evaluate(table_file, table_sep = DEF_TABLE_SEPARATOR):
    return compute_metrics(read_confusion_table(table_file, table_sep))


# Run evaluation of confusion table and print the metrics to the standard output device.
def print_evaluate(table_file, table_sep = DEF_TABLE_SEPARATOR, out_sep = DEF_OUTPUT_SEPARATOR):
    precision, recall, ari = evaluate(table_file, table_sep)

    print('%f%s%f%s%f' % (precision, out_sep, recall, out_sep, ari))


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('confusion_table', help = 'confusion table')
    argparser.add_argument('-t', '--table_separator', default = DEF_TABLE_SEPARATOR, help = 'separator between elements in the confusion-table file')
    argparser.add_argument('-o', '--output_separator', default = DEF_OUTPUT_SEPARATOR, help = 'separator between amplicon identifier and abundance')
    args = argparser.parse_args()

    print_evaluate(args.confusion_table, args.table_separator, args.output_separator)
