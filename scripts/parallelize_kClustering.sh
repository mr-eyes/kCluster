KCLUSTER=../scripts/kCluster_run.sh

# K21 Clustering

bash ${KCLUSTER} 1 25 relations/k21_pc.tsv names_maps/k21_pc.namesMap 21 &
bash ${KCLUSTER} 26 50 relations/k21_pc.tsv names_maps/k21_pc.namesMap 21 &
bash ${KCLUSTER} 51 75 relations/k21_pc.tsv names_maps/k21_pc.namesMap 21 &
bash ${KCLUSTER} 76 100 relations/k21_pc.tsv names_maps/k21_pc.namesMap 21 &

# K25 Clustering

bash ${KCLUSTER} 1 25 relations/k25_pc.tsv names_maps/k25_pc.namesMap 25 &
bash ${KCLUSTER} 26 50 relations/k25_pc.tsv names_maps/k25_pc.namesMap 25 &
bash ${KCLUSTER} 51 75 relations/k25_pc.tsv names_maps/k25_pc.namesMap 25 &
bash ${KCLUSTER} 76 100 relations/k25_pc.tsv names_maps/k25_pc.namesMap 25 &

# K31 Clustering

bash ${KCLUSTER} 1 25 relations/k31_pc.tsv names_maps/k31_pc.namesMap 31 &
bash ${KCLUSTER} 26 50 relations/k31_pc.tsv names_maps/k31_pc.namesMap 31 &
bash ${KCLUSTER} 51 75 relations/k31_pc.tsv names_maps/k31_pc.namesMap 31 &
bash ${KCLUSTER} 76 100 relations/k31_pc.tsv names_maps/k31_pc.namesMap 31 &
