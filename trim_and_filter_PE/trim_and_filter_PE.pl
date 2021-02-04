#!/usr/bin/perl

use Getopt::Std;

# This version is for PE data.

## Step 1:  Get input variables & check for missing data files or parameters.
getopts('1:2:a:b:c:d:m:o:q:h',\%option);
$errors = "";
if (exists($option{1})) {
  if ($option{1} eq "") { $errors .= "Empty required argument:  -1 <Mate1_fastq>\n"; }
  else {
    $mate1fq = $option{1};
    unless (-s $mate1fq) { $errors .= "User-specified input file does not exist:  $mate1fq\n"; }
  }
}
else { $errors .= "Missing required argument:  -1 <Mate1_fastq>\n"; }
if (exists($option{2})) {
  if ($option{2} eq "") { $errors .= "Empty required argument:  -2 <Mate2_fastq>\n"; }
  else {
    $mate2fq = $option{2};
    unless (-s $mate2fq) { $errors .= "User-specified input file does not exist:  $mate2fq\n"; }
  }
}
else { $errors .= "Missing required argument:  -2 <Mate2_fastq>\n"; }
if (exists($option{a})) {
  if ($option{a} eq "") { $errors .= "Empty required argument:  -a <Mate1_Pos1>\n"; }
  else { $mate1_p1 = $option{a}; }
}
else { $errors .= "Missing required argument:  -a <Mate1_Pos1>\n"; }
if (exists($option{b})) {
  if ($option{b} eq "") { $errors .= "Empty required argument:  -b <Mate1_PosN>\n"; }
  else { $mate1_p2 = $option{b}; }
}
else { $errors .= "Missing required argument:  -b <Mate1_PosN>\n"; }
if (exists($option{c})) {
  if ($option{c} eq "") { $errors .= "Empty required argument:  -c <Mate2_Pos1>\n"; }
  else { $mate2_p1 = $option{c}; }
}
else { $errors .= "Missing required argument:  -c <Mate2_Pos1>\n"; }
if (exists($option{d})) {
  if ($option{d} eq "") { $errors .= "Empty required argument:  -d <Mate2_PosN>\n"; }
  else { $mate2_p2 = $option{d}; }
}
else { $errors .= "Missing required argument:  -d <Mate2_PosN>\n"; }
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
  print "Syntax  = ./trim_and_filter_PE.pl [-h] -1 <Mate1_fastq> -2 <Mate2_fastq> -a <Mate1_Pos1> -b <Mate1_PosN> -c <Mate2_Pos1> -d <Mate2_PosN> -m <MinimumAvgBQS> -q <BaseQualScale> -o <OutputNameRoot>\n\n";
  print "Required Input:\n";
  print "  -1 <Mate1_fastq>      Full path to ~.fastq files for mate 1.\n";
  print "  -2 <Mate2_fastq>      Full path to ~.fastq files for mate 2.\n";
  print "  -a <Mate1_Pos1>       Read base to retain as first position of mate1.\n";
  print "  -b <Mate1_PosN>       Read base to retain as last position of mate1.\n";
  print "  -c <Mate2_Pos1>       Read base to retain as first position of mate2.\n";
  print "  -d <Mate2_PosN>       Read base to retain as last position of mate2.\n";
  print "  -m <MinimumAvgBQS>    Minimum average base quality score.\n";
  print "  -q <BaseQualScale>    Input data encoded in illumina or sanger scale for base qualities?\n";
  print "  -o <OutputRootName>   Root name of output files.\n\n";
  print "Other parameters:\n";
  print "  -h                    Print this help info.\n\n";
  print "Output:\n";
  print "A new ~.fastq file will be generated for each mate.  The file name will include the retained\n";
  print "region & minimum average base quality score used for filtering.  The file name format is:\n";
  print "  <OutputNameRoot>.1.trim_<Mate1_Pos1>_<Mate1_PosN>.minQS_<MinimumAvgBQS>.fastq\n";
  print "  <OutputNameRoot>.2.trim_<Mate2_Pos1>_<Mate2_PosN>.minQS_<MinimumAvgBQS>.fastq\n";
  print "A file describing how many read pairs failed the filter (& why) will be written to:\n";
  print "  <OutputNameRoot>.mate1_<Mate1_Pos1>_<Mate1_PosN>.mate2_<Mate2_Pos1>_<Mate2_PosN>.minQS_<MinimumAvgBQS>.FilterStats.txt\n\n\n";
  exit;
}



## Step 2:  Open input & output fastq files.  Initialize counters.

if (substr($mate1fq, -3) eq ".gz") { open(IN1, "zcat $mate1fq |"); } elsif (substr($mate1fq, -4) eq ".bz2") { open(IN1, "bzcat $mate1fq |"); } else { open(IN1, "$mate1fq"); }
if (substr($mate2fq, -3) eq ".gz") { open(IN2, "zcat $mate2fq |"); } elsif (substr($mate2fq, -4) eq ".bz2") { open(IN2, "bzcat $mate2fq |"); } else { open(IN2, "$mate2fq"); }

$outfq1 = "$root.1.trim_"."$mate1_p1"."_$mate1_p2.minQS_$minAQS.fastq";
$outfq2 = "$root.2.trim_"."$mate2_p1"."_$mate2_p2.minQS_$minAQS.fastq";
open(OUT1, ">$outfq1");
open(OUT2, ">$outfq2");

$mate1fail = 0; $mate2fail = 0; $mateBfail = 0; $passtot = 0; $readtot = 0;



## Step 3:  Trim & filter reads.  Write passing reads to output.
$cont = 1;
$mate1_len = $mate1_p2 - $mate1_p1 + 1;
$mate2_len = $mate2_p2 - $mate2_p1 + 1;
while ($cont == 1) {
  $line1a = <IN1>; $line2a = <IN1>; $line3a = <IN1>; $line4a = <IN1>;
  $line1b = <IN2>; $line2b = <IN2>; $line3b = <IN2>; $line4b = <IN2>;
  chomp $line2a; chomp $line4a;
  if ($line1a eq "") { $cont = 0; next; }
  $mate1trimS = substr($line2a, ($mate1_p1-1), $mate1_len);
  $mate1trimQ = substr($line4a, ($mate1_p1-1), $mate1_len);
  $sum = 0;
  for ($i=0; $i<length($mate1trimQ); $i++) { 
    $b = substr($mate1trimQ, $i, 1);
    if ($qscale eq "ilmn") { $sum += ord($b)-64; }
    if ($qscale eq "sang") { $sum += ord($b)-33; }
  }
  $avg1 = $sum/length($mate1trimQ);
  chomp $line2b; chomp $line4b;
  $mate2trimS = substr($line2b, ($mate2_p1-1), $mate2_len);
  $mate2trimQ = substr($line4b, ($mate2_p1-1), $mate2_len);
  $sum = 0;
  for ($i=0; $i<length($mate2trimQ); $i++) { 
    $b = substr($mate2trimQ, $i, 1);
    if ($qscale eq "ilmn") { $sum += ord($b)-64; }
    if ($qscale eq "sang") { $sum += ord($b)-33; }
  }
  $avg2 = $sum/length($mate2trimQ);
  $readtot++;
  if ($avg1 < $minAQS && $avg2 < $minAQS) { $mateBfail++; }
  if ($avg1 >= $minAQS && $avg2 < $minAQS) { $mate2fail++; }
  if ($avg1 < $minAQS && $avg2 >= $minAQS) { $mate1fail++; }
  if ($avg1 >= $minAQS && $avg2 >= $minAQS) {
    $passtot++;
    print OUT1 "$line1a$mate1trimS\n+\n$mate1trimQ\n";
    print OUT2 "$line1b$mate2trimS\n+\n$mate2trimQ\n";
  }
}
close(IN1); close(IN2); close(OUT1); close(OUT2);



## Step 4:  Write summary stats.
$statfile = "$root.mate1_"."$mate1_p1"."_$mate1_p2.mate2_"."$mate2_p1"."_$mate2_p2.minQS_$minAQS.FilterStats.txt";
open(STAT, ">$statfile");
print STAT "Total input read pairs\t$readtot\n";
$p = sprintf("%.2f", 100*$passtot/$readtot);
print STAT "Read pairs passing filter\t$passtot\t$p%\n";
$p = sprintf("%.2f", 100*$mate1fail/$readtot);
print STAT "Failed filter, mate 1 < avgQS threshold\t$mate1fail\t$p%\n";
$p = sprintf("%.2f", 100*$mate2fail/$readtot);
print STAT "Failed filter, mate 2 < avgQS threshold\t$mate2fail\t$p%\n";
$p = sprintf("%.2f", 100*$mateBfail/$readtot);
print STAT "Failed filter, both mates < avgQS threshold\t$mateBfail\t$p%\n";
close(STAT);

