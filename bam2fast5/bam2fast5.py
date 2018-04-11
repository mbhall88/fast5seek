import os
import sys
import warnings
import logging
import pysam
from typing import List, Set, Generator, Tuple

# suppress annoying warning coming from h5py
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from ont_fast5_api.fast5_file import Fast5File


def get_fast5_read_id(fast5: Fast5File, filepath: str) -> str:
    """Extracts the read id from a given fast5 file."""
    read_id = ''
    read_id_map = list(fast5.status.read_id_map.keys())
    if len(read_id_map) == 1:
        read_id = read_id_map[0]
    elif len(read_id_map) > 1:
        logging.warning(" More than one read id found for {}\n"
                        "Skipping this file.".format(filepath))
    else:  # no read id
        logging.warning(" No read id found for {}\n"
                        " Skipping this file.".format(filepath))
    return read_id


def get_fast5_run_id(fast5: Fast5File, filepath: str) -> str:
    """Extracts the run id from a given fast5 file."""
    run_id = fast5.get_tracking_id().get('run_id', '')
    if run_id == '':
        logging.warning(" No run id found for {}\nFile can still be "
                        "used if read id is present".format(filepath))
    return run_id


def get_read_and_run_id(filepath: str) -> Tuple[str, str]:
    """Extracts the read_id and run_id from a given fast5 file.
    If the file cannot be opened, it will be skipped.

    :param filepath: Path to the fast5 file.

    :returns The read_id and run_id as a tuple
    """
    read_id = ''
    run_id = ''
    try:
        fast5 = Fast5File(filepath)
        read_id = get_fast5_read_id(fast5, filepath)
        run_id = get_fast5_run_id(fast5, filepath)
    except OSError:  # issue trying to open fast5 file
        logging.error(" Error when trying to open {}\n"
                      "Skipping...".format(filepath))

    return read_id, run_id


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
    logging.info(" Looking for reads to extract in the following ref files:")
    for ref_file in [os.path.split(ref_path)[1] for ref_path in ref_path_list]:
        logging.info(" {}".format(ref_file))

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
    logging.info(" Found {} unique read ids".format(len(read_ids)))

    return read_ids


def get_fastq_read_ids(ref_path: str) -> Set[str]:
    """Extracts the read ids from a fastq file."""
    read_ids = set()
    with pysam.FastxFile(ref_path) as fastq:
        for entry in fastq:
            read_ids.add(entry.name.strip())

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
        with pysam.FastxFile(this_ref) as fastq:
            for entry in fastq:
                if not entry.comment:  # no run id information in fastq header
                    continue
                comments = entry.comment.split(' ')
                for field in comments:
                    if field.startswith('runid='):
                        fastq_run_ids.add(field.strip().replace('runid=', ''))
    if fastq_run_ids:
        logging.debug(" Found the following fastq run ids:")
        for run_id in sorted(fastq_run_ids):
            logging.debug("\t- {}".format(run_id))

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


def collect_fast5_filepaths(fast5_dir: List[str]) -> Set[str]:
    """Get a list of all fast5 filepaths under the given directory/directories.

    :param fast5_dir: Directory or directories to search recursively under.

    :returns Set of all the fast5 filepaths under fast5_dir
    """
    fast5_filepaths = set()
    for this_path in fast5_dir:
        new_pathset = scantree(this_path, ext='.fast5')
        fast5_filepaths |= set(new_pathset)
    logging.info(" Found {0} fast5 files.".format(len(fast5_filepaths)))

    return fast5_filepaths


def main(args):
    """The main method for the program.
    Steps:
    1) Determine the outfile or write to stdout if none given
    2) Get a set of read ids that map to the reference (fastq, sam, bam).
    3) Determine if there are any fastq files with run ids in header.
       This is used to avoid conflicts with fast5 file names and run ids.
    4) Get a set of fast5 filepaths to look through.
    5) Iterate through fast5 filepaths and save paths containing a mapped read.
    """
    # Step 1
    out_file = args.output or sys.stdout

    # Step 2
    read_ids = extract_read_ids(args.reference, args.mapped)

    # Step 3
    fastq_run_ids = get_fastq_run_ids(args.reference)

    # Step 4
    fast5_filepaths = collect_fast5_filepaths(args.fast5_dir)

    # Step 5
    final_filepath_list = []

    i = 0  # the number of files we looked at
    for filepath in fast5_filepaths:
        i += 1
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
