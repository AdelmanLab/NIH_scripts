## extract_fragments

Within this folder are included text files for chromosome sizes (*dm3_chr_size.txt*, *hg19_chr_size.txt*, *hg38_chr_size.txt*, *mm9_chr_size.txt*, *mm10_chr_size.txt*). These can be used with extract_fragments.pl to make sure the bedGraph files do not include bins that extend past the end of a chromosome, thus generating an error when uploaded to UCSC.

### Usage: extract_fragments.pl

Given a paired end bowtie output file, the program produces a BED file, a standard bedGraph file, or both, which contain those fragments that fall between a user-specified minimum and maximum length.  The BED file contains an entry for every fragment found, while the bedGraph file contains the count of fragment centers falling in bins of user specified length.

The output file produced is in standard bedGraph format, with 0-based start coordinate and 1-based end coordinate

The program is run with the following UNIX command line statement:
```
perl extract_fragments.pl [options] [input file name] [output file prefix] [instrument: miseq|hiseq]
```

### Options:
```
-a	Input is in SAM format.  The file must be sorted by query name rather than mapping coordinate.
-o	Desired output files.  Acceptable arguments are “b”, specifying BED output only, “g”, specifying BEDGRAPH output only, and “a” specifying both BED and BEDGRAPH output.  If a desired output file type is not specified, the default value of “b” will be used.  Example:  -o g.
-min	Minimum fragment length to extract.  If a minimum is not specified, the default value of 100 will be used.  Example:  -min 50.
-max	Maximum fragment length to extract. If a maximum is not specified, the default value of 200 will be used.  Example:  -max 500.
-b	Desired bin size.  The specified bin size may be any number greater than zero.  If the bin size is not specified, the default value of 25 will be used.  Example:  -b 50.
-l	Specify chromosome length list file.   The lengths specified by the file will be used to prevent the final bin of each chromosome from extending beyond the chromosome’s 3’ end.  The file must consist of two columns separated by white space, the first containing chromosome names, and the second, corresponding chromosome lengths.  Example:  -l hg18_chromosome_list.
-r	Specify reference genome.  This option may be used as an alternative to specifying a chromosome length list file (-l), and will automatically fetch the list of chromosomes and lengths from UCSC.  Example: -r hg18.
-D	Disable fixing of chromosome names.  By default, “chr” is appended to all chromosome names lacking that prefix, and any chromosome matching “mit” is replaced with “chrM”.  This is intended to improve compatibility with the UCSC browser when alignments are performed utilizing an index with bare numbers for chromosome names, such as the default fly genome.  Setting this option disables that behavior.
-M	Specify the search string utilized to identify the mitochondrial genome.  By default, any chromosome matching “mit” is replaced with “chrM”, the user may provide a suitable alternative.
-u	Remove duplicate reads exceeding specified maximum value. Examples: -u 1 (true de-duplication), -u 5.
```

### Input File Name:
The name of the bowtie output file to be processed.

#### Output File Prefix:
The specified prefix will be the basis of each output file name:  *[output file prefix].bed* and *[output file prefix].bedgraph*.  Additionally, this prefix will be used to generate track names: *[output file prefix]_bed* and *[output file prefix]_bedGraph*.

#### Sample Commands:
```perl extract_fragments.pl bowtie_file output_file```
Produces the file *output_file.bed*, containing all fragments of length 100 to 200, for data generated on a Miseq instrument.

```perl extract_fragments.pl -o g -b 15 bowtie_file output_file```
Produces the file *output_file.bedgraph*, containing the counts of all fragments of length 100 to 200 in each bin of 15 nucleotides.

```perl extract_fragments.pl -o g -b 15 -l dm3_chr_list bowtie_file output_file```
Produces the same output as the previous example, but uses the value contained in dm3_chr_list to truncate the final bin of each chromosome.  For dm3 chromosome 2L, this shifts the final bin’s end from position 23,011,545 to 23,011,544.

```perl extract_fragments.pl -o a -b 50 -min 200 -max 500 bowtie_file output_file```
Produces output_file.bed, containing all fragments of length 200 to 500, and output_file.bedgraph, containing the counts of all fragments of length 200 to 500 in each bin of 50 nucleotides.


