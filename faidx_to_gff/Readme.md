This function takes as input a fasta index (fai) file, and using tha creates a very simple gtf under some *heavy* assumptions.
It should be used only in cases of "emergency", i.e. when you have the sequence of a transcriptome, but not the genome, nor any idea on how the genes are divided into exons. If knowing exonic structure is important for you, then you shouldn't be using this function.

The assumptions are the following:

1) Each entry of the fasta index file is a gene. All genes will be treated as monoexonic.


