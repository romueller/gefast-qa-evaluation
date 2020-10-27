#!/usr/bin/python

# Helper methods for adjusting configuration files of GeFaST

import os


# Create an empty configuration file.
def create_config(path):
    with open(path, 'w') as out_file:
        pass


# Copy an existing configuration file.
def copy_config(path, duplicate):
    os.system('cp %s %s' % (path, duplicate))


# Add a configuration entry (key-value pair).
# If the configuration already contains an entry with the key,
# the value is not overwritten in the configuration.
def add_to_config(path, key, value):
    with open(path, 'a') as out_file:
        out_file.write('%s=%s\n' % (key, str(value)))


# Remove configuration entries (by key).
def remove_from_config(path, keys):
    with open(path, 'r') as in_file:
        lines = []
        for line in in_file:
            if not line.startswith('#'):
                key = line[:line.find('=')]
                if key not in keys:
                    lines.append(line)
    with open(path, 'w') as out_file:
        for line in lines:
            out_file.write(line)


# Change the value of a configuration entry.
# If the configuration does not contain an entry with the key,
# the configuration is not changed.
def change_value(path, key, value):
    with open(path, 'r') as in_file:
        lines = []
        for line in in_file:
            k = line[:line.find('=')]
            if k == key:
                lines.append('%s=%s\n' % (key, str(value)))
            else:
                lines.append(line)
    with open(path, 'w') as out_file:
        for line in lines:
            out_file.write(line)


# Remove an configuration file.
def remove_file(path):
    if os.path.exists(path):
        os.remove(path)
