#!/usr/bin/python

# Constants and general helper methods

import numpy as np

# default paths to binaries of reoccurring tools
DEF_GEFAST_PATH = 'GeFaST'
DEF_USEARCH_PATH = 'usearch'
DEF_VSEARCH_PATH = 'vsearch'
DEF_SWARM_PATH = 'swarm'

# default performance-logging command and the associated output file header
DEF_LOG_CMD = '/usr/bin/time -f "%%e;%%M;%%C" -a -o %s'
DEF_LOG_HEADER = 'reads;task;tool;mode;threshold;time;memory;cmd'

# argument separator of the 'misc' column in a task list
MISC_ARGS_SEP = '|'


# ===== Handling of list files =====

# Detect the file type by checking the first character of the first line in the file.
# Thus, the method will fail when there are empty or comment lines at the beginning of the file.
# Currently supported file types: FASTA (>), FASTQ (@)
def detect_file_type(read_file):
    file_type = None
    with open(read_file, 'r') as in_file:
        first_line = in_file.readline()
        if first_line[0] == '>':
            file_type = 'fasta'
        elif first_line[0] == '@':
            file_type = 'fastq'
        else:
            raise ValueError('ERROR: Only FASTA and FASTQ files can be used as input. Please check file %s.' % read_file)

    return file_type


# Parse a reads list file into a list of tuples (<name>, <path>, <type>).
#  - <name>: some (useful) name for the reads file
#  - <path>: file path of the reads file
#  - <type>: file type as detected by detected_file_type(...), e.g. fasta
#
# The list file is expected to have the following structure:
#  - one entry per line
#  - entry syntax: <name>:<path>
#  - no colons in <name> or <path>
def parse_reads_list_file(list_file):
    reads = []
    with open(list_file, 'r') as in_file:
        for line in in_file:
            tokens = line.strip().split(':')
            tokens.append(detect_file_type(tokens[1]))
            reads.append(tokens)

    return reads


# Parse a reads list given as a comma-separated list of entries in a single string.
# The entry syntax and the returned list are as described for parse_reads_list_file(...).
def parse_reads_list(reads_list):
    reads = []
    for pair in reads_list.split(','):
        tokens = pair.split(':')
        tokens.append(detect_file_type(tokens[1]))
        reads.append(tokens)

    return reads


# Parse a references list file into a list of file paths.
#
# If the same reference is used with multiple read files,
# it is enough to specify the file path once together with the number of read files.
# The file path is then replicated to have equally long read file and reference lists.
def parse_references_list_file(list_file, num_read_files = None):
    references = []
    with open(list_file, 'r') as in_file:
        for line in in_file:
            references.append(line.strip())
    if len(references) == 1 and num_read_files is not None:
        references = [references[0]] * num_read_files

    return references


# Parse a references list given as a comma-separated list of entries in a single string.
# The returned list is as described for parse_references_list_file(...).
def parse_references_list(ref_list, num_read_files = None):
    references = ref_list.split(',')
    if len(references) == 1 and num_read_files is not None:
        references = [references[0]] * num_read_files

    return references


# Parse a tasks list file into a list of tuples (<name>, <tool>, <mode>, <config>, <thresholds>, <misc>).
#  - <name>: some (useful) name for the task
#  - <tool>: name of the tool executing the task
#  - <mode>: mode / main option of the executing tool (empty if not supported by a tool)
#  - <config>: file path to configuration file of the executing tool (empty if not supported / used by a tool)
#  - <thresholds>: threshold specification, a comma-separated list or a range notation of the form <from>:<to>:<by>
#  - <misc>: any other (tool-specific) argument to the task-executing tool
#
# A task consists of one or multiple 'runs', depending on the number of specified thresholds.
#
# If different refinement thresholds should be used for the different thresholds, then the different specifications
# have to be separated by double commas (,,).
#
# The list file is expected to have the following structure:
#  - one entry per line
#  - entry syntax: <name>;<tool>;<mode>;<config>;<thresholds>;<misc>
#  - no semicolons in <name>, <tool> etc.
def parse_tasks_list_file(list_file):
    tasks = []
    with open(list_file, 'r') as in_file:
        for line in in_file:
            if line[0] != '#':
                task, tool, mode, config, thresholds, misc = line.strip().split(';')

                if ':' in thresholds:
                    first, last, inc = list(map(float, thresholds.split(':')))
                    thresholds = list(np.arange(first, last + inc / 2.0, inc))
                else:
                    thresholds = list(map(float, thresholds.split(',')))

                tasks.append([task, tool, mode, config, thresholds, misc])

    return tasks


# Parse a tasks list given as a ||-separated list of entries in a single string.
# The entry syntax and the returned list are as described for parse_tasks_list_file(...).
def parse_tasks_list(tasks_list):
    tasks = []
    for task_string in tasks_list.split('||'):
        task, tool, mode, config, thresholds, misc = task_string.strip().split(';')

        if ':' in thresholds:
            first, last, inc = list(map(float, thresholds.split(':')))
            thresholds = list(np.arange(first, last + inc / 2.0, inc))
        else:
            thresholds = list(map(float, thresholds.split(',')))

        tasks.append([task, tool, mode, config, thresholds, misc])

    return tasks


# Extract the value of key-value pair in the misc column of task.
# Key-value pairs have the form <key>=<value> and are separated, by default, by MISC_ARGS_SEP.
# If there is no pair with the given key, None is returned.
def get_value(args_str, key, sep):
    pos = args_str.find('%s=' % key)
    if pos != -1:
        end_pos = args_str.find(sep, pos)
        end_pos = len(args_str) if end_pos == -1 else end_pos
        return args_str[pos + len(key) + 1 : end_pos]
    else:
        return None


# Perform a simple sanity check of the specified reads and tasks.
#
# Current checks:
# (a) If refinement thresholds are specified, there can only be one (used with all thresholds) or their number
#     has to match the number of thresholds.
# (b) If the task is to be executed by USEARCH, a correct mode and order information have to be given.
# (c) If the task is to be executed by VSEARCH, a correct mode and, for some modes, order information have to be given.
def check_task_sanity(reads, tasks):
    for task, tool, mode, config, thresholds, misc in tasks:

        # statement of refinement thresholds (if any)
        pos = misc.find('refinement_threshold=')
        if pos != -1:
            end_pos = misc.find(MISC_ARGS_SEP, pos)
            end_pos = len(misc) if end_pos == -1 else end_pos
            refinement_entries = misc[pos : end_pos].split(',')
            if len(refinement_entries) != 1 and len(refinement_entries) != len(thresholds):
                raise ValueError('ERROR: Invalid specification of refinement thresholds for task %s (options: none, one or the same number as thresholds).' % task)

        # usearch: cluster_opt and order specification
        if tool == 'usearch':
            if mode not in ['cluster_fast', 'cluster_smallmem']:
                raise ValueError('ERROR: Invalid USEARCH mode for task %s (options: cluster_fast, cluster_smallmem).' % task)
            if 'order' not in misc:
                raise ValueError('ERROR: No order was specified for USEARCH for task %s (options: length, size).' % task)

        # vsearch: cluster_opt and order specification
        if tool == 'vsearch':
            if mode not in ['cluster_fast', 'cluster_size', 'cluster_smallmem']:
                raise ValueError('ERROR: Invalid VSEARCH mode for task %s (options: cluster_fast, cluster_size, cluster_smallmem).' % task)
            if mode in ['cluster_smallmem']:
                order = get_value(misc, 'order', MISC_ARGS_SEP)
                if order is None:
                    raise ValueError('ERROR: No order was specified to create the appropriate input for VSEARCH (cluster_smallmem) for task %s (options: length, size).' % task)
                elif order == 'size' and 'usersort' not in misc:
                    raise ValueError('ERROR: "usersort" has to be specified when expecting abundance-sorted input (option "size") for VSEARCH (cluster_smallmem) for task %s.' % task)
