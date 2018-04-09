

import os
import re
import sys
import h5py
import time
import progressbar


def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")

def get_read_group(list_of_names):
    """Extracts the correct group name for the group containing the read_id"""
    for name in list_of_names:
        if re.search(r'Read_\d+$', name):
            return name


def get_read_and_run_id(filepath):
    """Extracts the read_id and run_id from a given fast5 file.
    If the file cannot be opened, it will be skipped.

    Args:
        filepath (str): Path to the fast5 file.

    Returns:
        (tuple[str, str]): The read_id and run_id

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
            try:
                run_id = fast5['UniqueGlobalKey/tracking_id/'].attrs['run_id']
            except KeyError as err: # skip malformed files that don't have read id
                print("""{} has a malformed read id. Skipping...""".format(
                    filepath), file=sys.stderr)
                raise err
            return fast5[group].attrs['read_id'].decode(), run_id.decode()  # the read_id

    except IOError as err:  # skip file if it cannot be opened
        print("""{} could not be opened. Skipping...""".format(
            filepath), file=sys.stderr)
        raise err

def _get_file_extension(filepath):
    """returns the file extension of filepath arg"""
    return os.path.splitext(filepath)[1][1:]

def extract_read_ids(ref_path_list):
    """Extracts the all the read ids from the fastq file.
    The read_id is the first part of the header and starts with an @

    Args:
        ref_path (str): Path to the fastq or BAM/SAM file to extract
        read_ids from.

    Returns:
        read_ids (set[str]): A set of all the read ids in the fastq file.

    """
    # get the file extension
    # It is possible to used mixed file types for this
    read_ids = set()
    for ref_path in ref_path_list:
        extension = _get_file_extension(ref_path)
        if extension in {'bam', 'sam'}:
            read_ids |= get_sam_read_ids(ref_path)
        elif extension in {'fq', 'fastq'}:
            read_ids |= get_fastq_read_ids(ref_path)
        else:
            raise Exception('{} is not a supported file format. Supported file '
                            'types are: .fastq, .sam, and .bam'.format(extension))
    return read_ids

def get_fastq_read_ids(ref_path):
    """Extracts the read ids from a fastq file."""
    read_ids = set()
    with open(ref_path, 'r') as ref:
        for line in ref:
            if line.startswith('@'):  # i.e if line is header
                # split the line on spaces, take the first element, remove @
                read_id = line.split(' ')[0].replace('@', '')
                read_ids.add(read_id)

    return read_ids

def _clean_sambam_id(inputname):
    """Sometimes there are additional characters in the fast5 names added
    on by albacore or MinKnow. They have variable length, so this
    attempts to clean the name to match what is stored by the fast5 files.

    There are 5 fields. The first has a variable length.
    [7x or 8x az09]-[4x az09]-[4x az09]-[4x az09]-[12x az09]

    0688dd3-160d-4e2c-8af8-71c66c8db127
    7e33249c-144c-44e2-af45-ed977f6972d8
    67cbf79c-e341-4d5d-97b7-f3d6c91d9a85
    """
    #just grab the first five things when splitting with dashes
    splitname = inputname.split("-")[0:5]
    #The last one might have extra characters, unknown. We're relying
    # on the 5th field to consistently have 12 characters to extract
    # the correct id
    splitname[4] = splitname[4][0:12]
    return "-".join(splitname)

def get_sam_read_ids(ref_path):
    """Extract the read ids from a BAM or SAM file."""
    import pysam

    read_ids = set()
    with pysam.AlignmentFile(ref_path, 'r', ignore_truncation=True) as ref:
        for read in ref:
            # query_name is the query template name
            # - sometimes there are additional characters after the id names
            #    so we should use python string processing to split them up
            cleanname = _clean_sambam_id(read.query_name)
            read_ids.add(cleanname)

    return read_ids

def get_fastq_run_ids(references):
    """This method returns a set of fastq run ids if there are any fastq files.
    If there aren't any fastq files, it just returns and empty set."""
    fastq_run_ids = set()
    for this_ref in references:
        extension = _get_file_extension(this_ref)
        if extension in {'fastq', 'fq'}:
            with open(this_ref, 'r') as ref:
                for line in ref:
                    if line.startswith('@'):  # i.e if line is header
                        line_as_list = line.split(' ')
                        for field in line_as_list:
                            if field.startswith('runid='):
                                fastq_run_ids.add(field.strip().replace('runid=', ''))

    return fastq_run_ids

def get_fast5_paths(fast5_dir):
    """Input is a single directory that could contain fast5 files.
    Output is a set containing all of the paths to the fast5 files in those directories.
    """
    # Programmer's notes (@conchoecia): This method only takes a single
    #  directory as an argument rather than a list of directories. The reason
    #  for this is that a single-element list containing a single filepath as a
    #  string is iterated through like "/path/to/file.txt".split("/"), where
    #  every slash is another element. Only using one path as a string avoids
    #  misinterpreting the path and trying to scan the root directory, '/".
    filepaths = []
    widgets = ["{0}     - Searched at least ".format(len(timestamp())*" "),
               progressbar.Counter(),
               " files. ",
               progressbar.Timer()]
    bar = progressbar.ProgressBar(widgets=widgets, max_value=progressbar.UnknownLength)
    i = 0
    for root, dirs, files in os.walk(fast5_dir):
        for this_file in files:
            i += 1
            bar.update(i)
            if this_file.endswith(".fast5"):
                filepath = os.path.join(root, this_file)
                filepaths.append(filepath)
    # forces the pointer to the next line
    print(file=sys.stderr)
    return set(filepaths)

def main(args):
    """The main method for the program. Runs each command independently.
    Steps:
    1) Collect the arguments.
    2) Determine the outfile
    3) Get a set of read ids that map to the reference.
    4) Get a set of fastq run ids if the references contain any fastq files.
    5) Get a set of fast5 filepaths to look through in step 6.
    6) Iterate through fast5 filepaths and save paths containing a mapped read.
    """

    # Step 2, determine the outfile
    # if no output is given, write to stdout
    out_file = args.output or sys.stdout

    # Step 3, get a set of the read ids in the fastq or BAM/SAM file
    #  use the list comprehension to nicely print out the list of references.
    print("{0} - Looking for reads to extract in the following ref files:".format(
        timestamp()), file=sys.stderr)
    for ref_file in [os.path.split(x)[1] for x in args.reference]:
        print("{0}   - {1}".format(
            len(timestamp())*" ", ref_file), file=sys.stderr)

    read_ids = extract_read_ids(args.reference)
    print("{0} - We found {1} unique read ids".format(timestamp(), len(read_ids)),
          file=sys.stderr)
    print("{0} - To check for formatting normalcy, here are the first five:".format(
        timestamp()), file=sys.stderr)
    for i in range(5):
        print("{0}   - {1}".format(
            len(timestamp())*" ", list(read_ids)[i]), file=sys.stderr)

    # Step 4, figure out if there is a fastq file in the list of references.
    #  If so, make a set of all the run ids.
    #  This is used to avoid conflicts with fast5 files names and run ids.
    fastq_run_ids = get_fastq_run_ids(args.reference)
    # if there are any fastq run ids, print them out nicely
    if fastq_run_ids:
        print("{0} - We found the following fastq run ids:".format(
            timestamp()), file=sys.stderr)
        for run_id in sorted(fastq_run_ids):
            print("{0}   - {1}".format(
                len(timestamp())*" ", run_id), file=sys.stderr)

    # Step 5
    # collect all of the filepaths into a (very large) set object
    # We can run through this iteratively or through parallelization
    print("{0} - Now collecting fast5 filepaths to process.".format(
        timestamp()), file=sys.stderr)
    fast5_filepaths = set()
    for this_path in args.fast5_dir:
        print("{0}   - {1}".format(
            len(timestamp())*" ", this_path), file=sys.stderr)
        new_pathset = get_fast5_paths(this_path)
        print("{0}     - Found {1} .fast5 files.".format(
            len(timestamp())*" ", len(new_pathset)), file=sys.stderr)
        fast5_filepaths |= new_pathset
    print("{0}   - Found {1} .fast5 files total.".format(
        len(timestamp())*" ", len(fast5_filepaths)), file=sys.stderr)

    # Step 6
    # Iterate through fast5 filepaths and save paths containing a mapped read.
    final_filepath_list = []
    bar = progressbar.ProgressBar()
    print("{0} - Now looking through all .fast5 paths for matches.".format(
        timestamp()), file=sys.stderr)
    i = 0 # the number of files we looked at
    for filepath in bar(sorted(fast5_filepaths)):
        i += 1
        bar.update(i)
        try:
            fast5_read_id, fast5_run_id = get_read_and_run_id(filepath)
        except Exception:  # file cannot be opened or is malformed
            continue

        # if the file is from a mapped read
        if fast5_read_id in read_ids:
            #print("found a match, {}".format(fast5_read_id))
            # if fastq, make sure read and fastq are from the same
            # experiment. Small chance read ids could be the same from
            # diff. experiments. If not, skip file.
            if fastq_run_ids:
                if fast5_run_id not in fastq_run_ids:
                    continue
            final_filepath_list.append(filepath)

    # Step 7
    # Convert the final filepath list to a set and print
    final_filepath_set = set(final_filepath_list)
    for thisfile in sorted(final_filepath_set):
        print(thisfile, file=out_file)

    print("\n{0} - {1} files found.".format(
        timestamp(), len(final_filepath_set)),
          file=sys.stderr)

