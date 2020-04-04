This function takes as input a fasta index (fai) file, and using tha creates a very simple gtf with some *heavy* limitations.
It should be used only in cases of "emergency", i.e. when you have the sequence of a transcriptome, but not the genome, nor any idea on how the genes are divided into exons. If knowing exonic structure is important for you, then you shouldn't be using this function.

The **limitaitons** are the following:

1) We expect that each entry of the fasta index file is a gene. All genes will be treated as monoexonic. Each entry of fai file will be explained by two lines. One for the feature "gene", and one for the feautre "exon". Assuming monoexonic entries, the content for the two features will be identical.

2) Column 6 (Score) will always be set as 1000. This can be changed inside the function. 

3) Column 7 (strand) will always be coded as positive (+)

4) Column 8 (frame) will always be set to unknown (".") 

5) Column 9 (attribute) will always be of the form: gene_id=name_of_my_gene.1; Parent=name_of_my_gene 


