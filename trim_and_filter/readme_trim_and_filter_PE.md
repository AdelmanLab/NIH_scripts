These scripts can be used for quality filtering and trimming of single-end or paired-end FASTQs:

Usage:

Syntax  = ./filter_SE.pl [-h] -i <input_fastq> -m <MinimumAvgBQS> -q <BaseQualScale> -o <OutputNameRoot>

Required Input:
  -i <input_fastq>      Full path to input ~.fastq file.
  -m <MinimumAvgBQS>    Minimum average base quality score.
  -q <BaseQualScale>    Input data encoded in illumina or sanger scale for base qualities? Accepted values are "sanger" or "illumina".
  -o <OutputRootName>   Root name of output files.

Other parameters:
  -h                    Print this help info.

Output:
A new ~.fastq file will be generated.  The file name will include the retained region &
minimum average base quality score used for filtering.  The file name format is:
  <OutputNameRoot>.minQS_<MinimumAvgBQS>.fastq
A file describing how many reads failed the filter will be written to:
  <OutputNameRoot>.minQS_<MinimumAvgBQS>.FilterStats.txt

Usage:

Syntax  = ./filter_PE.pl [-h] -1 <Mate1_fastq> -2 <Mate2_fastq> -m <MinimumAvgBQS> -q <BaseQualScale> -o <OutputNameRoot>

Required Input:
  -1 <Mate1_fastq>      Full path to ~.fastq files for mate 1.
  -2 <Mate2_fastq>      Full path to ~.fastq files for mate 2.
  -m <MinimumAvgBQS>    Minimum average base quality score.
  -q <BaseQualScale>    Input data encoded in illumina or sanger scale for base qualities?
  -o <OutputRootName>   Root name of output files.

Other parameters:
  -h                    Print this help info.

Output:
A new ~.fastq file will be generated for each mate.  The file name will include the retained
region & minimum average base quality score used for filtering.  The file name format is:
  <OutputNameRoot>.minQS_<MinimumAvgBQS>.1.fastq
  <OutputNameRoot>.minQS_<MinimumAvgBQS>.2.fastq
A file describing how many read pairs failed the filter (& why) will be written to:
  <OutputNameRoot>.minQS_<MinimumAvgBQS>.FilterStats.txt

Usage:

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

Output:
A new ~.fastq file will be generated.  The file name will include the retained region &
minimum average base quality score used for filtering.  The file name format is:
  <OutputNameRoot>.trim_<Keep_Pos1>_<Keep_PosN>.minQS_<MinimumAvgBQS>.fastq
A file describing how many reads failed the filter will be written to:
  <OutputNameRoot>.trim_<Keep_Pos1>_<Keep_PosN>.minQS_<MinimumAvgBQS>.FilterStats.txt

Usage:

Syntax  = ./trim_and_filter_PE.pl [-h] -1 <Mate1_fastq> -2 <Mate2_fastq> -a <Mate1_Pos1> -b <Mate1_PosN> -c <Mate2_Pos1> -d <Mate2_PosN> -m <MinimumAvgBQS> -q <BaseQualScale> -o <OutputNameRoot>

Required Input:
  -1 <Mate1_fastq>      Full path to ~.fastq files for mate 1.
  -2 <Mate2_fastq>      Full path to ~.fastq files for mate 2.
  -a <Mate1_Pos1>       Read base to retain as first position of mate1.
  -b <Mate1_PosN>       Read base to retain as last position of mate1.
  -c <Mate2_Pos1>       Read base to retain as first position of mate2.
  -d <Mate2_PosN>       Read base to retain as last position of mate2.
  -m <MinimumAvgBQS>    Minimum average base quality score.
  -q <BaseQualScale>    Input data encoded in illumina or sanger scale for base qualities?
  -o <OutputRootName>   Root name of output files.

Other parameters:
  -h                    Print this help info.

Output:
A new ~.fastq file will be generated for each mate.  The file name will include the retained
region & minimum average base quality score used for filtering.  The file name format is:
  <OutputNameRoot>.1.trim_<Mate1_Pos1>_<Mate1_PosN>.minQS_<MinimumAvgBQS>.fastq
  <OutputNameRoot>.2.trim_<Mate2_Pos1>_<Mate2_PosN>.minQS_<MinimumAvgBQS>.fastq
A file describing how many read pairs failed the filter (& why) will be written to:
  <OutputNameRoot>.mate1_<Mate1_Pos1>_<Mate1_PosN>.mate2_<Mate2_Pos1>_<Mate2_PosN>.minQS_<MinimumAvgBQS>.FilterStats.txt
