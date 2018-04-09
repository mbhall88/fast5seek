"""
Outputs paths of all the fast5 files from a
given directory that are contained within a fastq or BAM/SAM file.

Usage:

It's pretty straight-forward to use:

    ./bam2fast5 -i <fast5_dir> -r <in.fastq|in.bam|in.sam> -o <out.txt>

The script will walk down into subdirectories as well, so you can just give it your directory containing all your files.

What it does is read in `<in.fastq|in.bam|in.sam>` and extract the read id from each header. It then goes through all the fast5 files under `<fast_dir>` and checks whether their read id is in the set of read ids from `<in.fastq|in.bam|in.sam>`. If it is, the path to the file is written to it's own line in `<out.txt>`.

Please see the github page for more detailed instructions. https://github.com/mbhall88/bam2fast5/

Contributors:
Michael Hall (github@mbhall88)
Darrin Schultz (github@conchoecia)
"""
import argparse
import os
import sys
import logging
from bam2fast5 import bam2fast5



def setup_logging(level):
    """Sets up the logging based on cli option passed."""
    logging_levels = {
        0: "NOTSET",
        1: "CRITICAL",
        2: "ERROR",
        3: "WARNING",
        4: "INFO",
        5: "DEBUG"
    }
    log_level = logging_levels.get(level)
    logging.basicConfig(level=log_level,
                        format='[%(asctime)s]:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))


class FullPathsList(argparse.Action):
    """Expand user- and relative-paths"""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                [os.path.abspath(os.path.expanduser(value)) for value in
                 values])


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--fast5_dir",
        action=FullPathsList,
        help="""Directory of fast5 files you want to query. Program will
        walk recursively through subdirectories.""",
        type=str,
        nargs="+",
        required=True)

    parser.add_argument(
        "-r", "--reference",
        action=FullPathsList,
        help="""Fastq or BAM/SAM file.""",
        nargs="+",
        required=True)

    parser.add_argument(
        "-o", "--output",
        action=FullPaths,
        help="""Filename to write fast5 paths to. If nothing is entered,
        it will write the paths to STDOUT.""",
        default=None)

    parser.add_argument(
        "--log_level",
        help="Level of logging. 0 is none, 5 is for debugging. Default is 4 "
             "which will report info, warnings, errors, and critical "
             "information.",
        default=4,
        type=int,
        choices=range(6))

    args = parser.parse_args()

    setup_logging(args.log_level)
    logging.info(" Starting bam2fast5.")
    bam2fast5.main(args)
    logging.info(" Done with bam2fast5. Bye.")

if __name__ == '__main__':
    sys.exit(cli())
