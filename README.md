# fast5_in_fastq
Yes, another Nanopore script. Because the world doesn't have enough of them.

    git clone https://github.com/mbhall88/fast5_in_fastq.git
    cd fast5_in_fastq

This script can be used to generate a list of file paths for `fast5` files that are contained within a `fastq` file of interest.

It's pretty straight-forward to use:

    ./fast5_in_fastq -i <fast5_dir> -f <in.fastq> -o <out.txt>

The script will walk down into subdirectories as well, so you can just give it your directory containing all your files.

If no output (`-o`) is given, it will write the output to `stdout`. 

So if you wanted to pipe these paths into another program, you could do something like

    mkdir subset_dir/
    ./fast5_in_fastq -i </path/to/fast5s/> -f <in.fastq> | xargs cp -t subset_dir/


The above example would copy the fast5 files that are found in your fastq to `subset_dir/`.

If there are any issues (with the program) let me know. 
