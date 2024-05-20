#!/usr/bin/perl

use Getopt::Std;

# This version is for SE data.

## Step 1:  Get input variables & check for missing data files or parameters.
getopts('i:a:b:m:o:q:h',\%option);
$errors = "";
if (exists($option{i})) {
  if ($option{i} eq "") { $errors .= "Empty required argument:  -i <input_fastq>\n"; }
  else {
    $inputfq = $option{i};
    unless (-s $inputfq) { $errors .= "User-specified input file does not exist:  $inputfq\n"; }
  }
}
else { $errors .= "Missing required argument:  -i <input_fastq>\n"; }
if (exists($option{a})) {
  if ($option{a} eq "") { $errors .= "Empty required argument:  -a <Keep_Pos1>\n"; }
  else { $keep_p1 = $option{a}; }
}
else { $errors .= "Missing required argument:  -a <Keep_Pos1>\n"; }
if (exists($option{b})) {
  if ($option{b} eq "") { $errors .= "Empty required argument:  -b <Keep_PosN>\n"; }
  else { $keep_p2 = $option{b}; }
}
else { $errors .= "Missing required argument:  -b <Keep_PosN>\n"; }
if (exists($option{m})) {
  if ($option{m} eq "") { $errors .= "Empty required argument:  -m <MinimumAvgBQS>\n"; }
  else { $minAQS = $option{m}; }
}
else { $errors .= "Missing required argument:  -m <MinimumAvgBQS>\n"; }
if (exists($option{o})) {
  if ($option{o} eq "") { $errors .= "Empty required argument:  -o <OutputNameRoot>\n"; }
  else { $root = $option{o}; }
}
else { $errors .= "Missing required argument:  -o <OutputNameRoot>\n"; }
if (exists($option{q})) {
  if ($option{q} eq "") { $errors .= "Empty required argument:  -q <BaseQualScale>\n"; }
  else {
    if ($option{q} eq "illumina") { $qscale = "ilmn"; }
    elsif ($option{q} eq "sanger") { $qscale = "sang"; }
    else { $errors .= "BaseQualScale must be defined as \'illumina\' or \'sanger\'.\n"; }
  }
}
else { $errors .= "Missing required argument:  -q <BaseQualScale>\n"; }


if (exists $option{h}) { $errors  = "Print help.\n"; }

if ($errors ne "") {
  if ($errors !~ /^Print help/) { print "\nERRORS...\n$errors\n"; }
  print "\nUsage:\n\n";
  print "Syntax  = ./trim_and_filter_SE.pl [-h] -i <input_fastq> -a <Keep_Pos1> -b <Keep_PosN> -m <MinimumAvgBQS> -q <BaseQualScale> -o <OutputNameRoot>\n\n";
  print "Required Input:\n";
  print "  -i <input_fastq>      Full path to input ~.fastq file.\n";
  print "  -a <Keep_Pos1>        Read base to retain as first position.\n";
  print "  -b <Keep_PosN>        Read base to retain as last position.\n";
  print "  -m <MinimumAvgBQS>    Minimum average base quality score.\n";
  print "  -q <BaseQualScale>    Input data encoded in illumina or sanger scale for base qualities?\n";
  print "  -o <OutputRootName>   Root name of output files.\n\n";
  print "Other parameters:\n";
  print "  -h                    Print this help info.\n\n";
  print "Output:\n";
  print "A new ~.fastq file will be generated.  The file name will include the retained region &\n";
  print "minimum average base quality score used for filtering.  The file name format is:\n";
  print "  <OutputNameRoot>.trim_<Keep_Pos1>_<Keep_PosN>.minQS_<MinimumAvgBQS>.fastq\n";
  print "A file describing how many reads failed the filter will be written to:\n";
  print "  <OutputNameRoot>.trim_<Keep_Pos1>_<Keep_PosN>.minQS_<MinimumAvgBQS>.FilterStats.txt\n\n\n";
  exit;
}



## Step 2:  Open input & output fastq files.  Initialize counters.

if (substr($inputfq, -3) eq ".gz") { open(IN1, "zcat $inputfq |"); } else { open(IN1, "$inputfq"); }
$outfq = "$root.trim_"."$keep_p1"."_$keep_p2.minQS_$minAQS.fastq";
$statfile = "$root.trim_"."$keep_p1"."_$keep_p2.minQS_$minAQS.FilterStats.txt";
open(OUT1, ">$outfq");
$fail = 0; $passtot = 0; $readtot = 0;



## Step 3:  Trim & filter reads.  Write passing reads to output.
$cont = 1;
$readlen = $keep_p2 - $keep_p1 + 1;
while ($cont == 1) {
  $line1a = <IN1>; $line2a = <IN1>; $line3a = <IN1>; $line4a = <IN1>;
  chomp $line2a; chomp $line4a;
  if ($line1a eq "") { $cont = 0; next; }
  $keeptrimS = substr($line2a, ($keep_p1-1), $readlen);
  $keeptrimQ = substr($line4a, ($keep_p1-1), $readlen);
  $sum = 0;
  for ($i=0; $i<length($keeptrimQ); $i++) { 
    $b = substr($keeptrimQ, $i, 1);
    if ($qscale eq "ilmn") { $sum += ord($b)-64; }
    if ($qscale eq "sang") { $sum += ord($b)-33; }
  }
  $avg = $sum/length($keeptrimQ);
  $readtot++;
  if ($avg < $minAQS) { $fail++; }
  else {
    $passtot++;
    print OUT1 "$line1a$keeptrimS\n+\n$keeptrimQ\n";
  }
}
close(IN1); close(OUT1); 



## Step 4:  Write summary stats.
open(STAT, ">$statfile");
print STAT "Total input reads\t$readtot\n";
$p = sprintf("%.2f", 100*$passtot/$readtot);
print STAT "Reads passing filter\t$passtot\t$p%\n";
$p = sprintf("%.2f", 100*$fail/$readtot);
print STAT "Reads failing filter\t$fail\t$p%\n";
close(STAT);

