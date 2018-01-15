# bam2fast5

This is a script designed for Nanopore (ONT) data that generates a
list of file paths for `fast5` files that are contained within
`fastq`, `BAM`, or `SAM` files of interest.

The output can be piped to another unix command to copy the fast5
files to a new directory or to save a list of file paths. This is
great for situations when you only need to share some, but not all,
reads of an ONT dataset.

### Which Python do I need?

This script should work for both python 2 and python 3.

### Installation

To install the script with git, you can clone it with github. Then,
change to the directory where you installed it (`cd`) and make the
script executable (`chmod`)

    git clone https://github.com/mbhall88/bam2fast5.git
    cd bam2fast5 && chmod +x bam2fast5

### Usage

It's pretty straight-forward to use:

    ./bam2fast5 -i <fast5_dir> -r <in.fastq|in.bam|in.sam> -o <out.txt>

The script will walk down into subdirectories as well, so you can just
give it your directory containing all your files.

What it does is read in the `<in.fastq|in.bam|in.sam>` files and
extract the read id from each header. It then goes through all the
fast5 files under `<fast_dir>` and checks whether their read id is in
the set of read ids from `<in.fastq|in.bam|in.sam>`. If it is, the
path to the file is written to it's own line in `<out.txt>`.

If no output (`-o`) is given, it will write the output to `stdout`.

#### Multiple Inputs

It is possible to use multiple directories/files as
arguments. No need to merge bam|fastq|sam files.

    ./bam2fast5 -i /myfast5/dir/1/ /other/fast5/dir/2/ -r reads.sorted.bam reads2.bam


For example, if all of your fast5 directories contain the prefix
`myfast5_` and the reference files contain `.sorted.bam`, you can use
wildcards to find them all if they are in the same directory.

    ./bam2fast5 -i myfast5_* -r *.sorted.bam

You can also mix reference file types in the arguments. For example,
if you happen to have a sam file, a fastq file, and a bam file that
contain reads you would like to find fast5 files for, they can all be
processed simultaneously like so.

    ./bam2fast5 -i myfast5_* -r mapped.sorted.bam mapped2.sam filtered_reads.fastq

#### Fastq

Currently this program does not support gzipped fastq files. Fastq
files can end with `.fq` or `.fastq`.

#### Piping Commands

So if you wanted to pipe these paths into another program, you could do something like

    mkdir subset_dir/
    ./bam2fast5 -i </path/to/fast5s/> -r <in.fastq> | xargs cp -t subset_dir/

The above example would copy the `fast5` files that are found in your `fastq` to `subset_dir/`.

#### Recommended Usage

However because of the computationally intensive step required to open
`fast5` files, we recommend that you first save the output of
`bam2fast5` to a file for safekeeping, then proceed with analysis like so:

    mkdir subset_dir/
    ./bam2fast5 -i /path/to/fast5s/ -r in.fastq > mapped_reads.txt
    cat mapped_reads.txt | xargs cp -t subset_dir/

For example, it took 37 minutes to look for mapped reads in 1.87
million fast5 files on a single processor. This same process took 10
minutes using a parallelized version of the program with 90 cores with
spinning disk drives. Faster processing speeds are possible if your
fast5 files are stored on solid state drives.

### Parallelization

We are currently developed a parallelized version of the program to
speed up the analysis. It is called `bam2fast5_parallelized` and
uses the same arguments as `bam2fast5`. So far it is 60-70% faster
than the single-processor version for large datasets.

### Dependencies
You need to have [`h5py`](https://github.com/h5py/h5py). You'll also need [`pysam`](https://github.com/pysam-developers/pysam) if you're going to be using SAM or BAM files.

    pip install h5py
    pip install pysam

### Contact

If there are any issues with the program please [open an issue above](https://github.com/mbhall88/bam2fast5/issues).

### Contributors

Michael Hall @mbhall88
Darrin Schultz @conchoecia
