# kCluster

## Trial 1
Using pairwise matrix, it succeeded at very small datasets, failed on protein coding and whole transcriptome due to the large number of targets (transcripts), the size of the matrix at that point will reach about 203K x 203K cell, and ordinary RAM would not handle this large matrix.
