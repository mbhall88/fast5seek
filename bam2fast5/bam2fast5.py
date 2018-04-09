import os
import re
import sys
import warnings
import logging
import pysam
from typing import List, Set, Generator

# suppress annoying warning coming from this libraries use of h5py
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import h5py


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
            except KeyError as err:  # skip malformed files that don't have read id
                print("""{} has a malformed read id. Skipping...""".format(
                    filepath), file=sys.stderr)
                raise err
            return fast5[group].attrs[
                       'read_id'].decode(), run_id.decode()  # the read_id

    except IOError as err:  # skip file if it cannot be opened
        print("""{} could not be opened. Skipping...""".format(
            filepath), file=sys.stderr)
        raise err


def _get_file_extension(filepath: str) -> str:
    """Returns the file extension of filepath arg without fullstop. If file is
    gzipped, drops the gz."""
    name, ext = os.path.splitext(filepath)
    if ext == '.gz':
        return os.path.splitext(name)[-1][1:]
    else:
        return ext[1:]


def extract_read_ids(ref_path_list: List[str], mapped: bool) -> Set[str]:
    """Extracts the all the read ids from the fastq file.
    The read_id is the first part of the header and starts with an @

    Args:
        ref_path_list: Paths to the fastq or BAM/SAM files to extract
        read_ids from.
        mapped: Only extract read ids from mapped reads in BAM/SAM files.

    Returns:
        read_ids: A set of all the read ids in the fastq file.

    """
    read_ids = set()
    for ref_path in ref_path_list:
        extension = _get_file_extension(ref_path)
        if extension in {'bam', 'sam'}:
            read_ids |= get_sam_read_ids(ref_path, mapped)
        elif extension in {'fq', 'fastq'}:
            read_ids |= get_fastq_read_ids(ref_path)
        else:
            logging.error(
                ' {0} is not a supported file format. Supported file types are:'
                ' fastq, sam, and bam.\n\tSkipping {1}'.format(extension,
                                                               ref_path))
    return read_ids


def get_fastq_read_ids(ref_path: str) -> Set[str]:
    """Extracts the read ids from a fastq file."""
    read_ids = set()
    with open(ref_path, 'r') as ref:
        for line in ref:
            if line.startswith('@'):  # i.e if line is header
                # split the line on spaces, take the first element, remove @
                read_id = line.split(' ')[0].replace('@', '')
                read_ids.add(read_id)

    return read_ids


def _clean_sambam_id(inputname: str) -> str:
    """Sometimes there are additional characters in the fast5 names added
    on by albacore or MinKnow. They have variable length, so this
    attempts to clean the name to match what is stored by the fast5 files.

    There are 5 fields. The first has a variable length.
    [7x or 8x az09]-[4x az09]-[4x az09]-[4x az09]-[12x az09]

    0688dd3-160d-4e2c-8af8-71c66c8db127
    7e33249c-144c-44e2-af45-ed977f6972d8
    67cbf79c-e341-4d5d-97b7-f3d6c91d9a85
    """
    # just grab the first five things when splitting with dashes
    splitname = inputname.split("-")[0:5]
    # The last one might have extra characters, unknown. We're relying
    # on the 5th field to consistently have 12 characters to extract
    # the correct id
    splitname[4] = splitname[4][0:12]
    return "-".join(splitname)


def get_sam_read_ids(ref_path: str, mapped: bool) -> Set[str]:
    """Extract the read ids from a BAM or SAM file.

    :param ref_path: Path to SAM/BAM file.
    :param mapped: Only extract read ids from mapped reads in BAM/SAM files.

    :return A set of read ids.
    """
    read_ids = set()
    with pysam.AlignmentFile(ref_path, 'r', ignore_truncation=True) as ref:
        for read in ref:
            read_is_mapped = not read.is_unmapped
            if mapped and read_is_mapped:
                cleanname = _clean_sambam_id(read.query_name)
            elif not mapped:  # use wants all reads, mapped or unmapped
                cleanname = _clean_sambam_id(read.query_name)
            else:  # only want mapped reads but this read is unmapped
                continue
            read_ids.add(cleanname)

    return read_ids


def get_fastq_run_ids(references: List[str]) -> Set[str]:
    """This method returns a set of fastq run ids if there are any fastq files.
    If there aren't any fastq files, it just returns and empty set."""
    fastq_run_ids = set()
    for this_ref in references:
        extension = _get_file_extension(this_ref)
        if extension not in {'fastq', 'fq'}:
            continue
        with open(this_ref, 'r') as ref:
            for line in ref:
                if not line.startswith('@'):
                    continue
                line_as_list = line.split(' ')
                for field in line_as_list:
                    if field.startswith('runid='):
                        fastq_run_ids.add(field.strip().replace('runid=', ''))
    return fastq_run_ids


def scantree(path: str, ext: str) -> Generator:
    """Recursively scans a directory and returns file paths ending in a given
    extension.
    :param path: Directory to scan.
    :param ext: Yield files with this extension.
    :returns Yields path to each file ending in extension.
    """
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            for nested_entry in scantree(entry.path, ext):
                yield nested_entry
        elif entry.is_file() and entry.name.endswith(ext):
            yield entry.path


def main(args):
    """The main method for the program. Runs each command independently.
    Steps:
    1) Determine the outfile
    2) Get a set of read ids that map to the reference.
    3) Get a set of fastq run ids if the references contain any fastq files.
    4) Get a set of fast5 filepaths to look through in step 6.
    5) Iterate through fast5 filepaths and save paths containing a mapped read.
    """
    # Step 1, determine the outfile
    # if no output is given, write to stdout
    out_file = args.output or sys.stdout

    # Step 2, get a set of the read ids in the fastq or BAM/SAM file
    logging.info(" Looking for reads to extract in the following ref files:")
    for ref_file in [os.path.split(ref_path)[1] for ref_path in args.reference]:
        logging.info(" {}".format(ref_file))

    read_ids = extract_read_ids(args.reference, args.mapped)
    logging.info(" Found {} unique read ids".format(len(read_ids)))
    logging.debug(" To check for formatting normality, here are the first five:")
    for i in range(5):
        logging.debug("   - {0}".format(list(read_ids)[i]))

    # Step 3, figure out if there is a fastq file in the list of references.
    #  If so, make a set of all the run ids.
    #  This is used to avoid conflicts with fast5 files names and run ids.
    fastq_run_ids = get_fastq_run_ids(args.reference)

    if fastq_run_ids:
        logging.debug(" Found the following fastq run ids:")
        for run_id in sorted(fastq_run_ids):
            logging.debug("\t- {}".format(run_id))

    # Step 4
    # collect all of the fast5 filepaths
    fast5_filepaths = set()
    for this_path in args.fast5_dir:
        new_pathset = list(scantree(this_path, ext='.fast5'))
        fast5_filepaths |= new_pathset
    logging.info(" Found {0} fast5 files.".format(len(fast5_filepaths)))

    # Step 6
    # Iterate through fast5 filepaths and save paths containing a mapped read.
    final_filepath_list = []
    bar = progressbar.ProgressBar()
    print("{0} - Now looking through all .fast5 paths for matches.".format(
        timestamp()), file=sys.stderr)
    i = 0  # the number of files we looked at
    for filepath in bar(sorted(fast5_filepaths)):
        i += 1
        bar.update(i)
        try:
            fast5_read_id, fast5_run_id = get_read_and_run_id(filepath)
        except Exception:  # file cannot be opened or is malformed
            continue

        # if the file is from a mapped read
        if fast5_read_id in read_ids:
            # print("found a match, {}".format(fast5_read_id))
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
