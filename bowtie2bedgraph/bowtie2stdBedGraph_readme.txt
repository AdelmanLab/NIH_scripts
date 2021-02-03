Usage: bowtie2bedgraph.pl [options] [input file] [output file prefix]
        -p              indicate input file is the result of a paired-end
                        alignment
        -a              indicate input is in SAM format (sorted by query name)
        -o [nbmv]       specify additional files to be generated: n=normalized,
                        b=binned, m=merged, v=binned and merged
        -b [integer>=1] specify bin size, requires -o b or -o v, default: 25
        -s [integer>=0] specify number of bps to shift ChIP-seq hits prior to
                        merging, requires -o m or -o v, default: 75
        -l [string]     provide list of chromosome names and lengths to prevent
                        hits from being shifted beyond the end of the
                        chromosome, requires -o b, -o m, or -o v
        -r [string]     rather than specifying a chromosome length list,
                        specify a UCSC reference genome identifier (e.g. mm9),
                        to have chromosome lengths fetched automatically,
                        requires -o b, -o m, or -o v
        -t [integer>=1] specify number of nt trimmed prior to alignment
        -x              request stranded output files be swapped,
                        forward->reverse and reverse->forward; for paired-end
                        input report counts of End 2 rather than End 1 5' ends
        -D              disable appending "chr" to beginning of chromosome names
                        and renaming of mitochondrial chromosome
        -M [string]     string used to match mitochondrial chromosome name,
                        default: "mit"
        -u [integer>=1] remove duplicate reads exceeding specified maximum
                        to retain per position
        -i              perform de-duplication in a strand-independent manner,
                        requires paired-end input (-p)

See bowtie2bedgraph.doc for most usage details.  The differences between bowtie2bedgraph and
bowtie2stdBedGraph are as follows:

-The output file produced is in standard bedGraph format, with 0-based start coordinate and 1-based end coordinate
-Paired-end alignments are accepted as input with the option -p
	-By default, the output produced will describe counts of End1 read 5' mapping locations, to instead
	 report the 5' mapping location of End2 reads, while swapping strand identifiers, use the -x option
	 (e.g. for startRNA, the 5prRNA file is produced by default, to automatically generate the 3prRNA
	 file, use -x along with paired-end input)
-SAM file input is accepted with the option -a
-De-duplication can be performed with the -u option.  Using -u 1 results in true removal of all duplicates, a higher
 value, -u N, results in the retention of up to N copies
	-For paired end input, two read pairs are considered duplicates if the 5' mapping locations of both End 1 and
	 End 2 are identical, for single-end input, reads with the same 5' mapping position are considered duplicates
	-For paired end input, de-duplication is performed in a strand-specific manner by default, i.e. two pairs aligned
	 such that the End 1 5' position of the first is the same as the End 2 5' position of the second with a similar
	 relationship between the End 2 5' position of the first and the End 1 5' position of the second, will not be
	 considered duplicates
	-To perform de-duplication in a strand-independent manner, use -i.  With this option, in the example described
	 above, the two pairs would be considered duplicates
