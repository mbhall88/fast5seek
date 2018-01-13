- 20180112 @conchoecia change notes
  - Changed the parallelized version of the program to use tqdm instead
    of progressbar. This shows the progress of parallelized processes.
  - Added a heuristic to find the optimum number of files per process.
    This should speed up the analyses drastically if there are many processes
    on the computer being used.

- 20171223 @conchoecia change notes
  - enable the user to input multiple reference files and multiple
    directories to search in
  - I did some preliminary testing of bottlenecks in the code and
    found that it only took 0.23 seconds to save a list of 68k files,
    but 1166 seconds to do the fast5 io and processing. This program
    would greatly benefit from parallelization.
  - Removed the extra conditional from :
     fast5_read_id in read_ids *and filepath not in filepaths*
     - Going to use a set to save the final filepaths if printing to output.
       Will only print output at the very end to ensure there are no duplicates.
       I will also use a set if printing to stdout to avoid this problem. This
       will also speed up the processing time.
  - Made more feedback for the user as well as a timed progressbar

- 20171221 @conchoecia change notes
  - In `extract_read_ids()`, removed the os.path.abspath() command since FullPaths
    action takes care of this.
  - made a `_get_file_extension(args.reference)` since that code is recycled a few times
  - removed `abs_path = os.path.abspath(os.path.realpath(args.fast5_dir))` since FUllPaths
    figures this out
  - Cleaned up the way that read ids are parsed from the fast5 files and dealt with
    inconsistencies in file naming.
  - Made timestamps and useful progress feedback for the user print to stderr.
  - Made the program work with both py2 and py3
  - refactored the code and made it more readable.
  - made all the messages say bam2fast5 since that makes more sense than fast5_in_ref (at least to me!)
