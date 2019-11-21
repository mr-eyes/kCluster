mkdir exp1; cd exp1/

#0 Create directories structure
mkdir -p ./{80,85,90,95}%/{k21,k25,k31}/{summaries,details}
mkdir -p kClusters/k{21,25,31}
mkdir -p names_maps
mkdir -p relations
mkdir -p visualizations/{80%,85%,90%,95%}/{k21,k25,k31}

#1 Download Human Transcriptome
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.transcripts.fa.gz

#2 Download GTF: Comprehensive gene annotation | PRI
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz

#3 Extract files
gunzip *.gz

#4 Generate protein_coding sequences
python ../scripts/filter_multifasta.py gencode.v28.transcripts.fa protein_coding

#5 Adding gene locus to protein coding transcripts fasta headers 
python ../scripts/add_locus.py -i protein_coding_gencode.v28.transcripts.fa -g gencode.v28.annotation.gtf -o loci_pc_gencode.v28.transcripts.fa
rm gencode.v28.annotation.gtf

#6 Download and extract clustering results of protein_coding
wget -c https://mr-eyes.github.io/kSpider/cdhit_benchmarking/data/cdhit_pc_clstrs.tar.gz
tar -xvzf cdhit_pc_clstrs.tar.gz

#7 Add clusterID to fasta files
for i in 80 85 90 95; 
do
python ../scripts/clstr_to_fasta.py loci_pc_gencode.v28.transcripts.fa cdhit_pc_transcripts_${i}.fa.clstr ${i};
mv clstr${i}*fa ${i}%/ ;
done;

#8 Delete unnecessary files
rm gencode.v28.transcripts.fa
rm protein_coding_gencode.v28.transcripts.fa
rm *.clstr *tar.gz

#9 Prepare data for indexing
python ../scripts/kProcessor_prepare.py loci_protein_coding_gencode.v28.transcripts.fa
rm loci_pc_gencode.v28.transcripts.fa

#10 kProcessor indexing of protein_coding
for  ksize  in 21 25 31;
do ./../kprocessor/Kprocessor index -i min_loci_protein_coding_gencode.v28.transcripts.fa -o ${ksize}_pc -k ${ksize} --names min_loci_protein_coding_gencode.v28.transcripts.fa.names --method MAP;
done

#11 Generate relations
for  ksize  in 21 25 31;
do python generate_relations.py ${ksize}_pc.map ${ksize}_pc.namesMap min_loci_protein_coding_gencode.v28.transcripts.fa ${ksize}_pc
done

#12 organize
mv *namesMap names_maps/
mv *tsv relations/
rm *map *fa *names

#13 Automate kClustering (Parallel Processing [All Threads])
bash ../scripts/parallelize_kClustering.sh