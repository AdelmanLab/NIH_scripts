## trim_and_filter_SE

This script is designed to enable trimming of single-end sequencing reads to a given length prior to mapping, and for selecting reads above a user defined average base quality score. 

Usage:
```
Syntax  = ./trim_and_filter_SE.pl [-h] -i <input_fastq> -a <Keep_Pos1> -b <Keep_PosN> -m <MinimumAvgBQS> -q <BaseQualScale> -o <OutputNameRoot>

Required Input:
  -i <input_fastq>      Full path to input ~.fastq file.
  -a <Keep_Pos1>        Read base to retain as first position.
  -b <Keep_PosN>        Read base to retain as last position.
  -m <MinimumAvgBQS>    Minimum average base quality score.
  -q <BaseQualScale>    Input data encoded in illumina or sanger scale for base qualities?
  -o <OutputRootName>   Root name of output files.

Other parameters:
  -h                    Print this help info.

```
Output:
A new ~.fastq file will be generated.  The file name will include the retained region &
minimum average base quality score used for filtering.  The file name format is:
```
  <OutputNameRoot>.trim_<Keep_Pos1>_<Keep_PosN>.minQS_<MinimumAvgBQS>.fastq
```
A file describing how many reads failed the filter will be written to:
```
  <OutputNameRoot>.trim_<Keep_Pos1>_<Keep_PosN>.minQS_<MinimumAvgBQS>.FilterStats.txt
```
