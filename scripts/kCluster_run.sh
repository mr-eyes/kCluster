START=$1
END=$2
RELATIONS_TSV=$3
NAMESMAP=$4
KMER=$5

for (( THRESHOLD=${START}; THRESHOLD<=${END}; THRESHOLD++ ));
do python ../scripts/relations_clustering.py ${RELATIONS_TSV} ${NAMESMAP} $THRESHOLD;
mv clusters_c${THRESHOLD}*k${KMER}*tsv kClusters/k${KMER}/
done