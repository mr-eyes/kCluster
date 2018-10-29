CD_HIT_CLSTR_FASTA=$1
KCLUSTERS_DIR=$2
KSIZE=$3
START=$4
END=$5

for (( THRESHOLD=${START}; THRESHOLD<=${END}; THRESHOLD++ ));
do 
    echo "Assessing K${KSIZE} | cdhit_${CD_HIT_CLSTR_FASTA:0:2}% | kCluster_${THRESHOLD}%"
    python ../scripts/kCluster_cdhit_compare.py ${CD_HIT_CLSTR_FASTA} ${KCLUSTERS_DIR}/clusters_c${THRESHOLD}.0_k${KSIZE}*tsv cdhit${CD_HIT_CLSTR_FASTA:0:2}_kCluster${THRESHOLD}%_k${KSIZE}_comparison.tsv; 
    mv cdhit${CD_HIT_CLSTR_FASTA:0:2}_kCluster${THRESHOLD}%_k${KSIZE}_comparison.tsv ${CD_HIT_CLSTR_FASTA:0:2}%/k${KSIZE}/details
    mv cdhit${CD_HIT_CLSTR_FASTA:0:2}_kCluster${THRESHOLD}%_k${KSIZE}_comparison_summary.txt ${CD_HIT_CLSTR_FASTA:0:2}%/k${KSIZE}/summaries
done


