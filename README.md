# fast5seek
[![PyPI status](https://img.shields.io/pypi/v/fast5seek.svg)](https://pypi.python.org/pypi/fast5seek)
[![Build Status](https://travis-ci.org/mbhall88/fast5seek.svg?branch=master)](https://travis-ci.org/mbhall88/fast5seek)
[![GitHub license](https://img.shields.io/github/license/mbhall88/fast5seek.svg)](https://github.com/mbhall88/fast5seek/blob/master/LICENSE)
![Twitter Follow](https://img.shields.io/twitter/follow/mbhall88.svg?style=social&logo=twitter&label=Follow)  

This program takes a directory (or multiple) of fast5 files along with any number 
of fastq, SAM, or BAM files. The output is the full paths for all fast5 files 
present in the fastq, BAM, or SAM files that are also in the provided fast5 directory(s).

### Installation

Using python3, run

```bash
pip3 install fast5seek
```

### Usage

It's pretty straight-forward to use:

    fast5seek -i /path/to/fast5s -r in.fastq in.bam in.sam -o out.txt

This will write all fast5 paths to a text file called `out.txt` - with each path 
on a new line.

What it does is read in the `in.fastq in.bam in.sam` files and
extracts the read id from each record. It then goes through all the
fast5 files under `/path/to/fast5s` and checks whether their read id is in
the set of read ids from `<in.fastq|in.bam|in.sam>`. If it is, the
path to the file is written to it's own line in `out.txt`.

If no output (`-o`) is given, it will write the output to `stdout`.

There is also an option to only search for read ids among mapped records in a 
BAM or SAM file - `-m/--mapped`. 

Gzipped fastq files can also be provided.

#### Full Usage

```
usage: fast5seek [-h] -i FAST5_DIR [FAST5_DIR ...] -r REFERENCE
                 [REFERENCE ...] [-o OUTPUT] [-m] [--log_level {0,1,2,3,4,5}]
                 [--no_progress_bar]

Outputs paths of all the fast5 files from a given directory that are contained within a fastq or BAM/SAM file.

Please see the github page for more detailed instructions.
https://github.com/mbhall88/fast5seek/

Contributors:
Michael Hall (github@mbhall88)
Darrin Schultz (github@conchoecia)

optional arguments:
  -h, --help            show this help message and exit
  -i FAST5_DIR [FAST5_DIR ...], --fast5_dir FAST5_DIR [FAST5_DIR ...]
                        Directory of fast5 files you want to query. Program
                        will walk recursively through subdirectories.
  -r REFERENCE [REFERENCE ...], --reference REFERENCE [REFERENCE ...]
                        Fastq or BAM/SAM file(s).
  -o OUTPUT, --output OUTPUT
                        Filename to write fast5 paths to. If nothing is
                        entered, it will write the paths to STDOUT.
  -m, --mapped          Only extract read ids for mapped reads in BAM/SAM
                        files.
  --log_level {0,1,2,3,4,5}
                        Level of logging. 0 is none, 5 is for debugging.
                        Default is 4 which will report info, warnings, errors,
                        and critical information.
  --no_progress_bar     Do not display progress bar.
```

#### Multiple Inputs

It is possible to use multiple directories/files as
arguments. No need to merge bam|fastq|sam files.

    fast5seek -i /myfast5/dir/1/ /other/fast5/dir/2/ -r reads.sorted.bam reads2.bam


For example, if all of your fast5 directories contain the prefix
`myfast5_` and the reference files contain `.sorted.bam`, you can use
wildcards to find them all if they are in the same directory.

    fast5seek -i myfast5_* -r *.sorted.bam


#### Piping Commands

If you wanted to pipe these paths into another program, you could do something like

    mkdir subset_dir/
    fast5seek -i /path/to/fast5s/ -r in.fastq | xargs cp -t subset_dir/

The above example would copy the `fast5` files that are found in your `fastq` to `subset_dir/`.

#### Recommended Usage

However because of the computationally intensive step required to open
`fast5` files, we recommend that you first save the output of
`fast5seek` to a file for safekeeping, then proceed with analysis like so:

    mkdir subset_dir/
    fast5seek -i /path/to/fast5s/ -r in.fastq -o mapped_reads.txt
    cat mapped_reads.txt | xargs cp -t subset_dir/

### Contact

If there are any issues with the program please [open an issue above](https://github.com/mbhall88/fast5seek/issues).

### Contributors

Michael Hall @mbhall88  
Darrin Schultz @conchoecia
