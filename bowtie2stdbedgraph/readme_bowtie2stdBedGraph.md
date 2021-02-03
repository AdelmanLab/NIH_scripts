## Bowtie2stdBedGraph

Within this folder are included text files for chromosome sizes (dm3_chr_size.txt, hg19_chr_size.txt, mm9_chr_size.txt). These can be used with bowtie2stdBedgraph to make sure the bedGraph files do not include bins that extend past the end of a chromosome, thus generating an error when uploaded to UCSC.

### Usage: bowtie2stdbedgraph.pl

Given a bowtie output file, the program produces a pair of bedgraph files containing the number of alignments per location, per strand.  The program will optionally produce paired files containing counts normalized to millions of total aligned reads and binned alignment counts.  Additionally, a single output file may be requested that contains the counts of both strands which have been shifted downstream, binned, and merged by summing the counts of corresponding locations  (as for ChIP-seq).

The program is run with the following UNIX command line statement:
```
perl bowtie2stdbedgraph.pl [options] [input file name] [desired track name]
```

#### Options:
```
Options:
-o  Desired output files. The program will always produce two files containing the 
    number of aligned reads on each strand. Additional files are requested through 
    the use of this option with an argument specifying the file type. Acceptable 
    arguments are “n”, specifying normalized files, “b”, specifying binned files, 
    “m” specifying a strand merged file, and “v” specifying a binned, strand merged 
    file. If multiple additional file types are desired, the arguments must consist 
    of a single string (single word) with no whitespace, for example: -o nb or -o nm.
    If desired, commas may be used between the characters requesting different file 
    types, which may make the command more readable:  -o n,b,m.
    
-b  Desired bin size. This option requires that binned output files be requested 
    through the use of the -o b or -o v option. The specified bin size may be any 
    number greater than one. If the bin size is not specified, the default value of
    25 will be used. Example:  -b 50.
    
-s  Desired shift value. This option requires that a merged output file be requested
    through the use of the -o m or -o v option. The specified shift value may be any
    number greater than zero. If the shift value is not specified, the default value
    of 75 will be used. Example:  -s 100.
    
-l  Specify chromosome length list file. This option requires that a merged output 
    file be requested through the use of the -o b,  -o m, or -o v option. The lengths 
    specified by the file will be used to prevent aligned reads from being shifted 
    beyond any chromosome’s 3’ end. These reads will instead be included in the last 
    valid bin. The file must consist of two columns separated by white space, the 
    first containing chromosome names, and the second, corresponding chromosome 
    lengths. Example:  -l hg19_chr_size
    
-r  Specify reference genome. This option requires that a binned or merged output file
    be requested through the use of the -o b, -o m, or -o v option. This option may be 
    used as an alternative to specifying a chromosome length list file (-l), and will 
    automatically fetch the list of chromosomes and lengths from UCSC. Example -r hg19
    
-t  Trim value. If sequences were trimmed for quality prior to alignment, this option 
    may be used to specify the number of nucleotides trimmed, and shift the 5’ positions
    of the output accordingly. If sequences were trimmed to remove non-genomic 
    sequence, this option should not be used. Example:  -t 5.
    
-x  Swap stranded output files. This option causes all alignments to the forward 
    strand to be placed in the _reverse.bedgraph output file, and all alignments to 
    the reverse strand to be placed in the _forward.bedgraph file. Use of this option
    ensures 3’ RNA alignments are associated with the correct strand.
    
-D  Disable fixing of chromosome names. By default, “chr” is appended to all 
    chromosome names lacking that prefix, and any chromosome matching “mit” is replaced
    with “chrM”. This is intended to improve compatibility with the UCSC browser when
    alignments are performed utilizing an index with bare numbers for chromosome names,
    such as the default fly genome. Setting this option disables that behavior.
    
-M  Specify the search string utilized to identify the mitochondrial genome. By 
    default, any chromosome matching “mit” is replaced with “chrM”, the user may 
    provide a suitable alternative.
```

#### Input File Name:
In addition to specifying the file to use as input, this file name is used as the prototype name for all output files.  These files will include the original file name, with additional information appended indicating the files’ contents, preceded by an underscore.  The appended information is as follows:
```
_forward		  Number of aligned reads, forward strand
_reverse		  Number of aligned reads, reverse strand
_forward_normal	  Normalized number of aligned reads, forward strand
_reverse_normal	  Normalized number of aligned reads, reverse strand
_forward_binned	  Number of aligned reads per bin, forward strand
_reverse_binned	  Number of aligned reads per bin, reverse strand
_merged		      Number of aligned reads per shifted, merged, bin
```

All output files will also have appended to them the extension .bedgraph.

#### Desired Track Name:
The track name will be included in the header of each file.  The information appended to the name of each file will be appended to each track name as well.

#### Sample Commands:
```
perl bowtie2stdbedgraph.pl bowtie_file fake_track
```
Produces two files named bowtie_file_forward.bedgraph and bowtie_file_reverse.bedgraph with respective track names fake_track_forward and fake_track_reverse.
  

```
perl bowtie2bedgraph.pl –o n bowtie_file fake_track
```
Produces four files:  the two described above, bowtie_file_forward_normal.bedgraph, and bowtie_file_reverse_normal.bedgraph.  The corresponding track names are fake_track_forward_normal and fake_track_reverse_normal.  
  

```
perl bowtie2stdbedgraph.pl -o nbm -b 50 -s 100 bowtie_file fake_track
```
Produces seven files:  the four described above, bowtie_file_forward_binned.bedgraph, bowtie_file_reverse_binned.bedgraph, and bowtie_file_merged.bedgraph.  The corresponding track names are fake_track_forward_binned, fake_track_reverse_binned, and fake_track_merged.

```
perl bowtie2stdbedgraph.pl -o nbm -b 50 -s 100 -l fake_lengths bowtie_file fake_track
```
Produces the seven files listed above.  If any aligned read is shifted beyond a chromosome’s 3’ end, as specified by the fake_lengths file, it will be included in the last valid bin. 

```
perl bowtie2stdbedgraph.pl -t 5 bowtie_file fake_track
```
Produces the two files from the first example, but the position values will be shifted by 5 in the negative direction for the forward strand, and shifted by 5 in the positive direction for the reverse strand.

