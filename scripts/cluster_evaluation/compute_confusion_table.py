#!/usr/bin/python

# Computes a confusion table from clusters and taxonomic assignments.
#
# Adapted from:
# Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014)
# Swarm: robust and fast clustering method for amplicon-based studies.
# PeerJ 2:e593 https://doi.org/10.7717/peerj.593
#
# Supplement 8: https://doi.org/10.7717/peerj.593/supp-8
# Section: 2.6.3 Confusion table

DEF_SEPARATOR = '_'


# Parse taxonomic assignments.
def parse_taxonomy(taxonomy_file, sep = DEF_SEPARATOR):
    amplicon_info = dict()
    taxa_list = []

    with open(taxonomy_file, 'r') as in_file:
        for line in in_file:
            amplicon, taxon = line.split()[0:2]
            try:
                id, abundance = amplicon.split(sep)
            except ValueError:
                id, abundance = amplicon, 1
            amplicon_info[id] = (abundance, taxon)
            taxa_list.append(taxon)

    return amplicon_info, taxa_list


# Dereplicate taxa and create a list of unique taxon names plus a dictionary for counting the taxa.
def dereplicate_taxa(taxa):
    taxa_list = list(set(taxa)) + ["Unassigned"]
    taxa_dict = dict(zip(taxa_list, [0] * len(taxa_list)))

    return taxa_list, taxa_dict


# Parse clusters and output one line of the confusion table.
def parse_cluster(amplicons, amplicon_info, taxa_dict, taxa_list, sep = DEF_SEPARATOR):
    for amplicon in amplicons:
        try:
            id, abundance = amplicon.split(sep)
        except ValueError:
            id, abundance = amplicon, 1

        try:
            abundance, taxon = amplicon_info[id]
        except KeyError:
            taxon = "Unassigned"

        taxa_dict[taxon] += int(abundance)

    return [str(taxa_dict[taxon]) for taxon in taxa_list]  # abundances per taxa in this cluster


# Compute the confusion table from the input files.
def compute_table(taxonomy_file, clusters_file, table_file, sep = DEF_SEPARATOR):
    amplicon_info, taxa_list = parse_taxonomy(taxonomy_file, sep)
    taxa_list, taxa_dict = dereplicate_taxa(taxa_list)

    with open(clusters_file, 'r') as in_file, open(table_file, 'w') as out_file:
        out_file.write('clusters_vs_taxa\t%s\n' % '\t'.join(taxa_list))

        for i, line in enumerate(in_file):
            amplicons = line.strip().split()
            abundances_per_taxa = parse_cluster(amplicons, amplicon_info, taxa_dict.copy(), taxa_list, sep)

            out_file.write('%i\t%s\n' % (i + 1, '\t'.join(abundances_per_taxa)))


if __name__ == '__main__':
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument('taxonomy', help = 'taxonomic assignments')
    argparser.add_argument('clusters', help = 'sequence clusters')
    argparser.add_argument('confusion_table', help = 'confusion table (output file)')
    argparser.add_argument('-s', '--separator', default = DEF_SEPARATOR, help = 'separator between amplicon identifier and abundance')
    args = argparser.parse_args()

    compute_table(args.taxonomy, args.clusters, args.confusion_table, args.separator)
