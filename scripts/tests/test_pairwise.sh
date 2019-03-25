set -euo pipefail

RUNLOG=runlog.log
echo "Run by `whoami` on `date`" > $RUNLOG # write log while running.

REF_FA=$1  # Reference Fasta
REF_NAMES=$2 # Reference Names

kCluster=/home/mabuelanin/Desktop/kprocessor/refseq_orthodb/kCluster
OUT=tmp_idx/idx_min_test

mkdir -p tmp_idx

export PATH=$PATH:/home/mabuelanin/Desktop/kprocessor/refseq_orthodb/skipmers_effect/kProcessor/build/apps/

# Indexing
echo "*** Indexing ***"
kProcessorApp index -i ${REF_FA} --names ${REF_NAMES} -k 21 -d 1 -m MAP -o ${OUT} 2>> ${RUNLOG}
echo "--------------------" >> ${RUNLOG}

# Building Pairwise
echo "*** Building pairwise matrix ***"
pypy ${kCluster}/scripts/generate_relations.py ${OUT}.map ${OUT}.namesMap 1>> ${RUNLOG}
mv idx*tsv tmp_idx
echo "--------------------" >> ${RUNLOG}

echo "Verifying..."
# Iterate over pairwise
sed 1d ${OUT}.tsv | while IFS=$'\t' read -r gene1 gene2 kmers norm
do 
    # Get genes header
    gene1header=$(grep -E "\s${gene1}$" ${OUT}.namesMap | awk '{print $1}')
    gene2header=$(grep -E "\s${gene2}$" ${OUT}.namesMap | awk '{print $1}')
    
    # Count shared kmers
    seqkit grep -r -p ${gene1header} ${REF_FA} | jellyfish count -t 8 -C -m 21 -s 1000 -o tmp_jelly_seq1 /dev/fd/0
    seqkit grep -r -p ${gene2header} ${REF_FA} | jellyfish count -t 8 -C -m 21 -s 1000 -o tmp_jelly_seq2 /dev/fd/0
    jellyfish dump -t tmp_jelly_seq1 | grep -E "^[^>]" | sort > tmp_kmers_1
    jellyfish dump -t tmp_jelly_seq2 | grep -E "^[^>]" | sort > tmp_kmers_2
    verified_kmers=$(join tmp_kmers_1 tmp_kmers_2 | wc -l)
    gene1_kmers=$(cat tmp_kmers_1 | wc -l)
    gene2_kmers=$(cat tmp_kmers_2 | wc -l)

    # Assert
    if [ $kmers -eq $verified_kmers ];
    then
        echo "CORRECT shared_kmers $gene1header & $gene2header = $verified_kmers" >> ${RUNLOG}
    else
        echo "WRONG shared_kmers $gene1header & $gene2header = $verified_kmers not $kmers" >> ${RUNLOG}
    fi


    # Verify similarity percentage kmers
    : '
    if [ $gene1_kmers -gt $gene2_kmers ];
    then
        verified_norm=$(echo "scale=2; ($verified_kmers/$gene2_kmers)*100" | bc)
    else
        verified_norm=$(echo "scale=2; ($verified_kmers/$gene1_kmers)*100" | bc)
    fi
    #verified_norm=$(printf '%.0f' "$verified_norm")
    #norm=$(printf '%.0f' "$norm")

    # Check if normalized values are equal
    if [ `echo "$verified_norm>$norm"|bc` -eq 1 ];
    then
        echo "CORRECT norm_value $gene1 & $gene2 = $verified_norm" >> ${RUNLOG}
    else
        echo "WRONG norm_value $gene1_kmers & $gene2_kmers = $verified_norm not $norm" >> ${RUNLOG}
    fi
    '

    echo "~~~~~" >> ${RUNLOG}
done

rm tmp_jelly_seq1 tmp_jelly_seq2 tmp_kmers_1 tmp_kmers_2