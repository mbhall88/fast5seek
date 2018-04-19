#!/usr/bin/env python
from __future__ import print_function

"""
Outputs paths of all the fast5 files from a
given directory that are contained within a fastq or BAM/SAM file.

Usage:

It's pretty straight-forward to use:

    ./fast5seek -i <fast5_dir> -r <in.fastq|in.bam|in.sam> -o <out.txt>

The script will walk down into subdirectories as well, so you can just give it your directory containing all your files.

What it does is read in `<in.fastq|in.bam|in.sam>` and extract the read id from each header. It then goes through all the fast5 files under `<fast_dir>` and checks whether their read id is in the set of read ids from `<in.fastq|in.bam|in.sam>`. If it is, the path to the file is written to it's own line in `<out.txt>`.

Contributors:
Michael Hall (github@mbhall88)
Darrin Schultz (github@conchoecia)
"""

import argparse
import os
import re
import sys
import h5py
import time
import progressbar

# Now import stuff for parallelization
from multiprocessing import cpu_count
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

class FullPathsList(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                [os.path.abspath(os.path.expanduser(value)) for value in values])

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--fast5_dir",
        action=FullPathsList,
        help="""Directory of fast5 files you want to query. Program will
        walk recursively through subdirectories.""",
        type=str,
        nargs = "+",
        required=True)

    parser.add_argument(
        "-r", "--reference",
        action=FullPathsList,
        help="""Fastq or BAM/SAM file.""",
        nargs = "+",
        required=True)

    parser.add_argument(
        "-o", "--output",
        action=FullPaths,
        help="""Filename to write fast5 paths to. If nothing is entered,
        it will write the paths to STDOUT.""",
        default=None)

    return parser.parse_args()

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

def available_threads():
    threads = cpu_count()
    if threads >= 90:
        return 90
    else:
        return threads

def parallel_process(array, function, n_jobs=2, use_kwargs=False,
                     front_num=3, disable = False, **kwargs):
    """
        A parallel version of the map function with a progress bar.
        A great implementation by http://danshiebler.com/2016-09-14-parallel-progress-bar/

        Args:
            array (array-like): An array to iterate over.
            function (function): A python function to apply to the elements of array
            n_jobs (int, default=16): The number of cores to use
            use_kwargs (boolean, default=False): Whether to consider the elements of array as dictionaries of
                keyword arguments to function
            front_num (int, default=3): The number of iterations to run serially before kicking off the parallel job.
                Useful for catching bugs
        Returns:
            [function(array[0]), function(array[1]), ...]

        For jobs that have a quick compute time and for many threads, it is
         probably more efficient to split this process up so that it runs in
         chunks rather than take a lot of overhead creating new jobs. This also
         needs to be optimized for list unwrapping
    """
    #We run the first few iterations serially to catch bugs
    if front_num > 0:
        front = [function(**a) if use_kwargs else function(a) for a in array[:front_num]]
    #If we set n_jobs to 1, just run a list comprehension. This is useful for benchmarking and debugging.
    if n_jobs==1:
        return front + [function(**a) if use_kwargs else function(a) for a in tqdm(array[front_num:])]
    #Assemble the workers
    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
        #Pass the elements of array into function
        if use_kwargs:
            futures = [pool.submit(function, **a) for a in array[front_num:]]
        else:
            futures = [pool.submit(function, a) for a in array[front_num:]]
        kwargs = {
            'total': len(futures),
            'unit': 'it',
            'unit_scale': True,
            'leave': True
        }
        #Print out the progress as tasks complete
        for f in tqdm(as_completed(futures), disable = disable, **kwargs):
            pass
    out = []
    #Get the results from the futures.
    for i, future in tqdm(enumerate(futures), disable=disable):
        try:
            out.append(future.result())
        except Exception as e:
            out.append(e)
    if front_num > 0:
        return front + out
    else:
        return out

def process_filepath_chunk(args):
    """Takes in a list of files and checks if they are one of the mapped files"""
    fast5_filepaths = args["fast5_filepaths"]
    read_ids = args["read_ids"]
    fastq_run_ids = args["fastq_run_ids"]

    final_filepath_list = []
    #print("in instance {}".format(i), file=sys.stderr)
    for filepath in fast5_filepaths:
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
    return final_filepath_list

def determine_pool_size(job_vector):
    """This function determines how large of a pool to make based on the
    system resources currently available and how many jobs there are to complete
    """
    available_threads = cpu_count()
    total_jobs = len(job_vector)
    threads_to_pass = total_jobs
    if total_jobs >= available_threads:
        threads_to_pass = available_threads
    return available_threads, total_jobs, threads_to_pass

def chunks(l, n):
    """Yield successive n-sized chunks from l.
    https://stackoverflow.com/questions/312443/"""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def process_a_2k_chunk(chunk_size):
    global fast5_filepaths
    global final_filepath_list
    chunk_size = int(chunk_size)
    instances = list(chunks(fast5_filepaths[:2000], chunk_size))
    available_threads, total_jobs, pool_size = determine_pool_size(instances)
    start = time.time()
    results = parallel_process([{'fast5_filepaths': instances[i],
                    'read_ids': read_ids,
                    'fastq_run_ids': fastq_run_ids }
                            for i in range(len(instances))],
                process_filepath_chunk, n_jobs = available_threads,
                use_kwargs = False, front_num=0, disable = True)
    final_filepath_list += [item for sublist in results for item in sublist]
    end = time.time()
    time_elapsed = end - start
    fast5_filepaths = fast5_filepaths[2000:]
    print("""{0} - It took {1:.5f}s per file with {2} files per process.""".format(
        len(timestamp())*" ", time_elapsed/2000, chunk_size), file=sys.stderr)
    return time_elapsed

def main():
    """The main method for the program. Runs each command independently.
    Steps:
    1) Collect the arguments.
    2) Determine the outfile
    3) Get a set of read ids that map to the reference.
    4) Get a set of fastq run ids if the references contain any fastq files.
    5) Get a set of fast5 filepaths to look through in step 6.
    6) Find the optimum chunk size to process all the filepaths given the number
        of available cores and the user's filepath structure.
    7) Iterate through fast5 filepaths and save paths containing a mapped read.
    8) Report how many reads were found
    """
    # Step 1, collect the arguments
    args = parse_arguments()

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

    global read_ids
    read_ids = extract_read_ids(args.reference)
    print("{0} - We found {1:,} unique read ids".format(timestamp(), len(read_ids)),
          file=sys.stderr)
    print("{0} - To check for formatting normalcy, here are the first five:".format(
        timestamp()), file=sys.stderr)
    for i in range(5):
        print("{0}   - {1}".format(
            len(timestamp())*" ", list(read_ids)[i]), file=sys.stderr)

    # Step 4, figure out if there is a fastq file in the list of references.
    #  If so, make a set of all the run ids.
    #  This is used to avoid conflicts with fast5 files names and run ids.
    global fastq_run_ids
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
    print("""{0}   - Please note that this is slow if the fast5 files are """.format(
        len(timestamp())*" "), file=sys.stderr)
    print("""{0}     in many directories.""".format(
        len(timestamp())*" "), file=sys.stderr)

    fast5_filepaths_set = set()
    for this_path in args.fast5_dir:
        print("{0}   - {1}".format(
            len(timestamp())*" ", this_path), file=sys.stderr)
        new_pathset = get_fast5_paths(this_path)
        print("{0}     - Found {1:,} .fast5 files.".format(
            len(timestamp())*" ", len(new_pathset)), file=sys.stderr)
        fast5_filepaths_set |= new_pathset
    print("{0}   - Sorting and converting unique filepaths to a list.".format(
        len(timestamp())*" "), file=sys.stderr)
    global fast5_filepaths
    global final_filepath_list

    final_filepath_list = []
    fast5_filepaths = list(sorted(fast5_filepaths_set))
    print("{0}   - Found {1:,} .fast5 files total.".format(
        len(timestamp())*" ", len(fast5_filepaths)), file=sys.stderr)

    # Step 6
    # Find the optimum chunk size to process all the filepaths given the number
    #  of available cores and the user's filepath structure.
    # For each chunk, the program will go through 2000 reads and time how long
    #  per read. The program will then select the fastest chunk size to proceed.

    # Use a binary search
    print("""{0} - Now finding the optimum number of fast5s per thread.""".format(
        timestamp()), file=sys.stderr)
    chunk_speed_x = [15, 20, 25, 50, 75, 100, 125, 150,
                     200, 225, 250, 275, 300, 400, 500]
    chunk_speed_y = []
    chunk_results_final = []
    start_all_fast5 = time.time()

    for this_chunk_speed in chunk_speed_x:
        len(fast5_filepaths) > 0
        chunk_speed_y.append(process_a_2k_chunk(this_chunk_speed))


    fastest_chunk_size = chunk_speed_x[chunk_speed_y.index(min(chunk_speed_y))]

    print("""{0} - {1} files per process was the fastest.""".format(
        len(timestamp())*" ", fastest_chunk_size), file=sys.stderr)

    # Step 7
    # Iterate through fast5 filepaths and save paths containing a mapped read.
    #process each file sequentially using max number of threads
    #determine the pool size to work with the unique sample names
    chunk_size = fastest_chunk_size
    instances = list(chunks(fast5_filepaths, chunk_size))
    available_threads, total_jobs, pool_size = determine_pool_size(instances)
    print("""{0} - Now opening all fast5 files to look for matches.""".format(
        timestamp()), file=sys.stderr)
    print("""{0} - There are {1:,} threads available.""".format(
        len(timestamp())*" ", available_threads), file=sys.stderr)
    print("""{0}   - There are {1:,} jobs.""".format(
        len(timestamp())*" ", total_jobs), file=sys.stderr)
    print("""{0}   - Making a pool with {1:,} threads""".format(
        len(timestamp())*" ", pool_size), file = sys.stderr)
    print("""{0}   - Please wait while we compile these jobs.""".format(
        len(timestamp())*" "), file = sys.stderr)

    results = parallel_process([{'fast5_filepaths': instances[i],
                        'read_ids': read_ids,
                        'fastq_run_ids': fastq_run_ids }
                                for i in range(len(instances))],
                    process_filepath_chunk, n_jobs = available_threads,
                    use_kwargs = False, front_num=3)
    print("""{0}   - Converting results to a usable form.""".format(
        len(timestamp())*" "), file = sys.stderr)
    final_filepath_list += [item for sublist in results for item in sublist]
    end_all_fast5 = time.time()
    time_elapsed = end_all_fast5 - start_all_fast5
    print("""{0}   - The .fast5 extraction and search took {1:.2f} seconds.""".format(
        len(timestamp())*" ", time_elapsed), file = sys.stderr)

    # Step 7
    # Convert the final filepath list to a set and print
    print("""{0} - Final Steps""".format(
        timestamp()), file=sys.stderr)
    print("""{0}   - Gathering unique file names.""".format(
        len(timestamp())*" "), file = sys.stderr)
    final_filepath_set = set(final_filepath_list)
    print("""{0}   - Printing unique filenames to stdout.""".format(
        len(timestamp())*" "), file = sys.stderr)
    print("""{0}     This step will take the longest if you have used """.format(
        len(timestamp())*" "), file = sys.stderr)
    print("""{0}     <fast5seek> | xargs cp -t subset_dir/ to copy files.""".format(
        len(timestamp())*" "), file = sys.stderr)
    for thisfile in sorted(final_filepath_set):
        print(thisfile, file=out_file)

    # Step 8
    # Final report on how many files found
    print("{0} - {1} files found.".format(
        timestamp(), len(final_filepath_set)),
          file=sys.stderr)
    print("""{0}   - This was {1:2f}% of the total reads mapped.""".format(
        len(timestamp())*" ", len(final_filepath_set)/len(read_ids)), file = sys.stderr)


if __name__ == '__main__':
    print("{0} - Starting fast5seek.".format(
        timestamp()), file=sys.stderr)
    main()
    print("{0} - Done with fast5seek. Bye.".format(
        timestamp()), file=sys.stderr)
