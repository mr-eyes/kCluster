for cdhit in 80 85 90 95
do
    for kmer in 21 25 31
    do
        for kCluster in 1 26 51 76 # Create four background processes
        do
            bash ../scripts/compare.sh ${cdhit}%/clstr${cdhit}_*.fa kClusters/k${kmer} ${kmer} ${kCluster} $((${kCluster} + 24)) &
        done
    done
done