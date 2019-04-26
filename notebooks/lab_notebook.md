# Lab Book :D

- added a modified index function with one more extra parameter (to change Q more freely)
  - new: colored_kDataFrame* index(string seqFileName,string namesFiles, uint64_t kSize, uint64_t Q);

## RUNs

### indexing: min_protein_coding_gencode.v28.transcripts.fa with K=31 & Q=31

time:  4min (indexing + writing on disk)
ram : 9489 MB
index-size: 8.8G

### indexing: min_protein_coding_gencode.v28.transcripts.fa with K=(25) & Q=32

time:  4.27 min (indexing + writing on disk)
ram : 12.12 GB
index-size: 11.9G

### indexing: min_protein_coding_gencode.v28.transcripts.fa with K=31 & Q=28

**FAILED SEGFAULT**

### indexing: min_protein_coding_gencode.v28.transcripts.fa with K=31 & Q=29

**FAILED SEGFAULT**

### indexing: min_protein_coding_gencode.v28.transcripts.fa with K=(31,25) & Q=30

**FAILED SEGFAULT**

### indexing: min_protein_coding_gencode.v28.transcripts.fa with K=(21,25) & Q=31

FAILED SEGFAULT

---

### virtualQs {'minQ': 1, 'maxQ': 31, 'stepQ': 5} - min_protein_coding_gencode.v28.transcripts.fa

Qs = 1,6,11,16,21,26,31
WALL-TIME: 27:41.60
RAM :      14.7 GB


### virtualQs {'minQ': 5, 'maxQ': 29, 'stepQ': 2} - min_protein_coding_gencode.v28.transcripts.fa

TIME:               50 min (ENTERED SWAP)
Total Memory :      15.1 GB
VirtualQs memory:   ~3 GB (will be enhanced)
