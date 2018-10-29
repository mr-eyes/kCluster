for cdhit in 80 85 90 95
do
    for kmer in 21 25 31
    do
        echo "Generating cdhit${cdhit}_k${kmer}.html"
        cd ${cdhit}%/k${kmer}/summaries/;
        find *txt | sort -V | xargs grep "" | python ../../../../scripts/visualize_cdhit_benchmarking.py cdhit${cdhit}_k${kmer}.html && mv cdhit${cdhit}_k${kmer}.html ../../../visualizations/${cdhit}%/k${kmer}/ &
        cd ../../../
    done
done