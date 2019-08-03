# Example 1

## Download test data

```bash

# Create directory
mkdir example1 && cd example1

# Download
wget test_data.zip

# Extract
unzip test_data.zip

```

Now you will have `seq.fa` and `seq.fa.names

## Description

seq.fa file contains two identical 100bp sequences with except in position 25,50,75 and 100.

running kCluster with kmer-size=25 demonstrates the effect of the processing with virtualQs 5,10,15,20 because in the normal scenario, there will be zero matches.

---

## Indexing

Indexing the Fasta file

`kCluster index_kmers -f seq.fa -n seq.fa.names -k 25`

That will generate multiple files with the prefix **idx_seq** that stores the index.

---

## Pairwise matrix generation

Here we performs pairwise distance matrix generation with minimum virtualQ = 5, maximum virtualQ = kSize = 25, and with a step of 5.
That means we will calculate the containment of kmers (virtualQs) of sizes [5,10,15,20,25] between each two sequences.
As we already have just only two sequences, with mismatches at every 25 nucleotides, that means there are *zero* shared kmers between the two sequences but there are shared virtualQs.

`kCluster pairwise --min-q 5 --max-q 25 --step-q 5 -i idx_seq`

The command will generate sqlite3 database file that stores all the information needed for clustering, including the pairwise similarity matrix.

If we dumped the file using the following command

`kCluster dump --db idx_seq_kCluster.sqlite`

will have this TSV formatted output with the following column names.
- `ID` : Pairwise distance record serial ID
- `seq1` : the ID of seq1, check the `namesMap` file for getting the original sequence header corresponding to that ID. 
- `seq2` : the ID of seq2.
- `min_kmers` : the shortest sequence contains 76 kmers (k = 25)
- `Q_X` : the virtualQ with its size.


**Output:**

| ID | seq1 | seq2 | min_kmers | Q_5 | Q_10 | Q_15 | Q_20 | Q_25 | 
|----|------|------|-----------|-----|------|------|------|------| 
| 1  | 1    | 2    | 76        | 59  | 40   | 25   | 15   | 0    | 

**Explanation**
- There are `59` virtualQs of size `5 bp` shared between `seq1` & `seq2`.
- There are `40` virtualQs of size `10 bp` shared between `seq1` & `seq2`.
- There are `25` virtualQs of size `15 bp` shared between `seq1` & `seq2`.
- There are `15` virtualQs of size `20 bp` shared between `seq1` & `seq2`.
- There are `zero` virtualQs of size `25 bp` shared between `seq1` & `seq2`.

---


## Clustering

You can perform clustering with any number of virtualQs that have been already generated in the pairwise matrix by determining the minimum, maximum and step of virtualQs.

In the following command, we will perform clustering on virtualQs [5,10,15,20,25] with similarity threshold `0`.

That means if there's at least 1% of the sequences is shared between both of the sequences they will be grouped together in one cluster.

`kCluster cluster -m 5 -M 25 -s 5 --db idx_seq_kCluster.sqlite --cutoff 0.0`

The output of the clustering will be just one group contains the sequences `1` and `2`.

**Increasing** the similarity threshold to 90% will split the two sequences into two separate clusters,

`kCluster cluster -m 5 -M 25 -s 5 --db idx_seq_kCluster.sqlite --cutoff 0.9`
