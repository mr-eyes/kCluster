# kCluster

## Trial 1
Using pairwise matrix, it succeeded at very small datasets, failed on protein coding and whole transcriptome due to the large number of targets (transcripts), the size of the matrix at that point will reach about 203K x 203K cell, and ordinary RAM would not handle this large matrix.

## Trial 2
Using NetworkX (Graph library) it will work, but consumes a lot of memory and time,
The whole transcriptome took about 10 GB RAM to be constructed in a form of Graph.

## Trial 3
Using hard coded function to get all connected components in a very short time. 
Clustering of the whole transcriptome took about 3 minutes to complete, from the kProcessor output to the clusters file generation.
