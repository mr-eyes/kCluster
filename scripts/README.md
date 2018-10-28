# Scripts

## CD-HIT Clstr file to FASTA file

- **Description:** Append clstr_id to FASTA header and write into new file.

- **Run:**  `python clstr_to_fasta.py <Fasta_file> <clstr_file>`

## Generate relations stats between each pair of transcripts

- **Description:** Convert the kProcessor output to a simple tsv representing shared kmers number and normalized value between each two transcripts.

- **Run:**  `python generate_relations.py <map_index_file> <names_map_file> <fasta_file> <kmer_size>`

## Cluster transcripts based on the their pre-generated relations

- **Description:** Cluster the relations into multiple clusters based on shared kmers number among transcripts with optional % threshold.

- **Run:**  `python generate_relations.py <relations_tsv_file> <names_map_file> <op: cuttof_threshold %>`