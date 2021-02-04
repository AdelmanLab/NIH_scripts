# Adelman/NIH_scripts
A collection of perl scripts for NGS analysis, written by Adam Burholder, Sara Grimm and David Fargo of the National Institute of Environmental Health Sciences (NIEHS). These scripts have been used in a number of collaborative studies betwen the Adelman lab and the NIEHS Center for Integrative Bioinformatics.

Examples:
Nechaev et al. Science 2010, Gilchrist et al. Genes Dev 2012, Core et al. Cell Reports 2012, Henriques et al. Mol Cell 2013, Scruggs et al. Mol Cell 2015, Williams et al. Mol Cell 2015, Henriques et al. Genes Dev 2018.


## things to put in this readme
These modular perl scripts allow a user to perform common tasks involved in visualization and analysis of NGS data. They are provided here for general use. 


### notes (BJM)
- most of the scripts had a readme text and/or Word file in its O2 folder. I've just migrated these over into README markdown files. Some of these are pretty simple (and for simple scripts) but seemed helpul to include this documentation.
- Q: regrading bowtie2stdbedGraph script documentation. It seems like this is an update to a previous bowtie2bedGraph script? The Word document ("bowtie2stdbedgraph.doc") seems to detail the old script, but the readme txt file ("bowtie2stdBedGraph_readme.txt") detailed the updates made in the new one. I've stitched these together into the README.md file, so I think this now covers all the options in the sciprt... but please check this readme carefully. Also, the documentation noted the chromosome length files, so I've left those in the folder.
- for the new Mac make_heatmap compilier, I don't know which versions of MacOS this is good for, so I just said "MacOS" without specifying the version.
- I made gitignore file because my MAC kept creating .DS_Store files and then pushes these to the repository (I did most of this locally). Once I've finished what I need to do on my MAC I can remove the gitignore file.
