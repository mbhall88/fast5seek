#!/usr/bin/env python

import argparse
import os
import re
import sys
import h5py


def get_fast5_dirs(fast5_dir, fastq, out_file):
    """Returns the paths to all fast5 files contained in the fastq file.

    Args:
            fast5_dir (str): The path to the directory where the fast5 files are.
            fastq (str): Path to the fastq file.
            out_file (File): File to write the fast5 paths contained
            in the fastq file to. 

    """
    # get a set of the read ids in the fastq file
    read_ids = extract_read_ids(fastq)

    # get the absolute path for the fast5 directory
    abs_path = os.path.abspath(os.path.realpath(fast5_dir))

    filepaths = set()
    for root, dirs, files in os.walk(abs_path):  # recursively walk fast5_dir
        for file_ in files:

            if file_.endswith(".fast5"):
                filepath = os.path.join(root, file_)
                
                try:
                    read_id = get_read_id(filepath)
                except IOError:  # file cannot be opened
                    continue

                # if the file is in the fastq file and it has not 
                # already been found...
                if read_id in read_ids and filepath not in filepaths:
                    filepaths.add(filepath)
                    out_file.write(filepath + '\n')

    sys.stderr.write('{} files found.\n'.format(len(filepaths)))


def get_read_group(list_of_names):
    """Extracts the correct group name for the group containing the read_id"""
    for name in list_of_names:
        if re.search(r'Read_\d+$', name):
            return name


def get_read_id(filepath):
    """Extracts the read_id from a given fast5 file.
    If the file cannot be opened, it will be skipped.

    Args:
            filepath (str): Path to the fast5 file.

    Returns:
            read_id (str): The read_id

    """
    group = 'Raw/Reads/'
    try:
        with h5py.File(filepath, 'r') as fast5:
            list_of_names = []

            # get all group names in the Raw/Reads group
            fast5['Raw/Reads'].visit(list_of_names.append)

            # extracts the group name that ends in the pattern 'Read_\d+'
            # where \d+ is any number of digits
            read_group = get_read_group(list_of_names)
            group += read_group

            return fast5[group].attrs['read_id']  # the read_id

    except IOError as err:  # skip file if it cannot be opened
        sys.stderr.write('{} could not be opened. \
            Skipping...\n'.format(filepath))
        raise err


def extract_read_ids(fastq_path):
    """Extracts the all the read ids from the fastq file.
    The read_id is the first part of the header and starts with an @

    Args:
            fastq_path (str): Path to the fastq file to extract read_ids from.

    Returns:
            read_ids (set[str]): A set of all the read ids in the fastq file.

    """
    # get the absolute path
    fastq = os.path.abspath(os.path.realpath(fastq_path))

    read_ids = set()
    with open(fastq_path, 'r') as fastq:
        for line in fastq:
            if line.startswith('@'):  # i.e if line is header
            	# split the line on spaces, take the first element, remove @
                read_id = line.split(' ')[0].replace('@', '')
                read_ids.add(read_id)

    return read_ids


def main():
    parser = argparse.ArgumentParser(
            description = "Outputs paths of all the fast5 files from a \
            given directory that are contained within a fastq file.",
            prog='fast5_in_fastq',
            )

    parser.add_argument(
            "-i", "--fast5_dir",
            help="Directory of fast5 files you want to query. Program will \
            walk recursively through subdirectories.",
            type=str, required=True)

    parser.add_argument(
            "-f", "--fastq",
            help = "Fastq file.",
            required=True)

    parser.add_argument(
            "-o", "--output",
            help = "Filename to write fast5 paths to. If nothing is entered, \
            it will write the paths to STDOUT.", 
            default=None)

    args = parser.parse_args()

    # if no output is given, write to stdout
    out_file = args.output or sys.stdout

    get_fast5_dirs(args.fast5_dir, args.fastq, out_file)


if __name__ == '__main__':
    main()
