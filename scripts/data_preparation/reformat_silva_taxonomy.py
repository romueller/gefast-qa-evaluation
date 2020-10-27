#!/usr/bin/python

# Scans the taxonomy file and limit the taxonomic information to the given level (default: 5 = genus).
# Remove entries with ambiguous taxonomic information or an unspecific label.

DEF_LEVEL = 5


# Remove taxonomic information from level (level + 1) onwards.
def remove_levels(tax_info, level):
    pos_spec = tax_info[1].find(';D_%i__' % (level + 1))
    if pos_spec == -1 and tax_info[1].endswith(';Ambiguous_taxa'):
        pos_spec = tax_info[1].rfind(';Ambiguous_taxa')

    tax_info[1] = tax_info[1][:pos_spec]

    return '\t'.join(tax_info)


# Create a reformatted and, potentially, reduced taxonomy file.
def reformat(taxonomy_file, output, level = DEF_LEVEL):
    with open(taxonomy_file, 'r') as in_file, open(output, 'w') as out_file:
        for line in in_file:
            line = remove_levels(line.rstrip('\n').split('\t'), level)
            if line.find('Ambiguous_taxa') == -1 and line.find('D_%i__uncultured' % level) == -1 and line.find('D_%i__' % level) != -1:
                out_file.write('%s\n' % line)


if __name__ == '__main__':
    import argparse
    argparser = argparse.ArgumentParser()
    argparser.add_argument('taxonomy', help = 'taxonomy file')
    argparser.add_argument('output', help = 'output file')
    argparser.add_argument('-l', '--level', type = int, default = DEF_LEVEL, help = 'highest level that has to be complete')
    args = argparser.parse_args()

    reformat(args.taxonomy, args.output, args.level)
