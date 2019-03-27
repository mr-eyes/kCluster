#set -euo pipefail

RUNLOG=runlog.log
echo "Run by `whoami` on `date`" > $RUNLOG # write log while running.

REF_FA=$1  # Reference Fasta
REF_NAMES=$2 # Reference Names
PREFIX=$3

# set to one in case want jellyFish kmer counting
RUN_JELLY=false

export PATH=$PATH:/home/mabuelanin/Desktop/kprocessor/refseq_orthodb/skipmers_effect/kProcessor/build/apps/

echo "Verifying..."

# Iterate over pairwise
sed 1d ${PREFIX}.tsv | while IFS=$'\t' read -r gene1 gene2 kmers norm
do 
    # Get genes header
    gene1header=$(grep -E "\s${gene1}$" ${PREFIX}.namesMap | awk '{print $1}')
    gene2header=$(grep -E "\s${gene2}$" ${PREFIX}.namesMap | awk '{print $1}')
    echo $gene1header
    echo $gene2header
    echo "-------------------------------"
    
    READS1=$(grep "${gene1header}$" ${REF_NAMES} | awk '{print $1}')
    READS2=$(grep "${gene2header}$" ${REF_NAMES} | awk '{print $1}')

    echo "" > tmp1.fa
    echo "" > tmp2.fa

    if ${RUN_JELLY};
    then
    echo "***Running kmer countin using JellyFish***"
    
    
    for value in $READS1
    do
        seqkit grep -r -p ${value} ${REF_FA} >> tmp1.fa
    done

    for value in $READS2
    do
        seqkit grep -r -p ${value} ${REF_FA} >> tmp2.fa
    done
    
    
    # Count shared kmers
    seqkit seq tmp1.fa | jellyfish count -t 8 -C -m 21 -s 1000 -o tmp_jelly_seq1 /dev/fd/0
    seqkit seq tmp2.fa  | jellyfish count -t 8 -C -m 21 -s 1000 -o tmp_jelly_seq2 /dev/fd/0
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

    fi

    echo "***Check by traversing the index itself***"
    # Check by traversing the index itself
    grep "[-,\,]$gene1," ${PREFIX}.map | awk -F'-' '{print $1}' | sort | uniq > seq1_colors
    grep "[-,\,]$gene2," ${PREFIX}.map | awk -F'-' '{print $1}' | sort | uniq > seq2_colors

    common_colors=$(join seq1_colors seq2_colors)
    echo $common_colors
    rm seq1_colors seq2_colors

    shared_kmers=0

    for value in $common_colors
    do
        cnt=$(grep ":${value}$" ${PREFIX}*map | wc -l)
        shared_kmers=$((shared_kmers+cnt))
    done

    if [ $shared_kmers -eq $kmers ];
    then
        echo "CORRECT [2] shared_kmers $gene1header & $gene2header = $kmers" >> ${RUNLOG}
    else
        echo "WRONG [2] shared_kmers $gene1header & $gene2header = $kmers not $shared_kmers" >> ${RUNLOG}
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

rm tmp_jelly_seq1 tmp_jelly_seq2 tmp_kmers_1 tmp_kmers_2 tmp1.fa tmp2.fa