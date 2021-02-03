## bedgraphs2stdBedGraph

Usage: 
```
bedgraphs2stdBedGraph [output prefix]
```
This script finds all standard bedGraphs and bedgraphs with 1-based start coordinates in
the working directory (inferring the format from the file name), and generates a merged
file in standard bedGraph format.  With a single file, this can be used to convert from
bedgraph with 1-based start to standard bedGraph.  A log file is also produced that includes
a list of the input files.
