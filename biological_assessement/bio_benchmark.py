transcript_to_locus = {}
transcript_to_gene = {}
locus_to_transcripts = {}
gene_to_transcripts = {}
locus_to_genes = {}
gene_to_locus = {}

cdhit_to_loci = {}
cdhit_to_genes = {}
transcript_to_cdhit = {}
cdhit_to_transcripts = {}
gene_to_cdhits = {}
locus_to_cdhits = {}

kCluster_to_loci = {}
kCluster_to_genes = {}
transcript_to_kCluster = {}
kCluster_to_transcripts = {}
gene_to_kCluster = {}
locus_to_kCluster = {}

kCluster_to_cdhits = {}
cdhit_to_kClusters = {}


with open("clstr95_loci_protein_coding_gencode.v28.transcripts.fa", "r") as fa:
    for line in fa:
        if line[0] != ">":
            continue

        line = line[1:-2].split("|")

        transcript = line[0]
        gene = line[1]
        locus = line[8]
        cdhit = line[9]

        transcript_to_gene[transcript] = gene
        transcript_to_locus[transcript] = locus
        gene_to_locus[gene] = locus
        transcript_to_cdhit[transcript] = cdhit

        # locus_to_genes
        if locus not in locus_to_genes:
            locus_to_genes[locus] = [gene]
        else:
            locus_to_genes[locus].append(gene)

        # locus_to_cdhits
        if locus not in locus_to_cdhits:
            locus_to_cdhits[locus] = [cdhit]
        else:
            locus_to_cdhits[locus].append(cdhit)

        # locus_to_transcripts
        if locus not in locus_to_transcripts:
            locus_to_transcripts[locus] = [transcript]
        else:
            locus_to_transcripts[locus].append(transcript)

        # gene_to_transcripts
        if gene not in gene_to_transcripts:
            gene_to_transcripts[gene] = [transcript]
        else:
            gene_to_transcripts[gene].append(transcript)

        # cdhit_to_genes
        if cdhit not in cdhit_to_genes:
            cdhit_to_genes[cdhit] = [gene]
        else:
            cdhit_to_genes[cdhit].append(gene)

        # cdhit_to_transcripts
        if cdhit not in cdhit_to_transcripts:
            cdhit_to_transcripts[cdhit] = [transcript]
        else:
            cdhit_to_transcripts[cdhit].append(transcript)

        # cdhit_to_loci
        if cdhit not in cdhit_to_loci:
            cdhit_to_loci[cdhit] = [locus]
        else:
            cdhit_to_loci[cdhit].append(locus)

        # gene_to_cdhits
        if gene not in gene_to_cdhits:
            gene_to_cdhits[gene] = [cdhit]
        else:
            gene_to_cdhits[gene].append(cdhit)


with open("cdhit95_kCluster100%_k25_comparison.tsv", "r") as kC:
    print(next(kC).split())
    for line in kC:
        line = line.split()
        q1 = line[1]
        q2 = line[2]
        kCluster = line[5]
        cdhits = line[6].split()
        kCluster_to_cdhits[kCluster] = cdhits


with open("clusters_c100.0_UNORD_PC.tsv", "r") as clusters:
    print next(clusters)
    for cluster in clusters:
        cluster = cluster.split()
        #print cluster
        kCluster = cluster[0]
        transcripts = cluster[1].split()
        kCluster_to_loci[kCluster] = [transcript_to_locus[tr] for tr in transcripts]
        kCluster_to_genes[kCluster] = [transcript_to_gene[tr] for tr in transcripts]


        for tr in transcripts:
            gene = transcript_to_gene[tr]
            locus = transcript_to_locus[tr]
            gene_to_kCluster[gene] = kCluster
            locus_to_kCluster[locus] = kCluster
            transcript_to_kCluster[tr] = kCluster

        kCluster_to_transcripts[kCluster] = transcripts


for kCluster, cdhits in kCluster_to_cdhits.iteritems():
    for cdhit in cdhits:
        if cdhit in cdhit_to_kClusters:
            cdhit_to_kClusters[cdhit] = [kCluster]
        else:
            cdhit_to_kClusters[cdhit].append(kCluster)
