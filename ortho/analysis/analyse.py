class Ensemble:
    """
    Ensemble Class : Analyise and process Ensemble Sequences
    """

    def __init__(self, file_name):
        self.file_name = file_name

    def transcriptome_gene_symbols(self):
        """
        return dict {transcript_id: gene_symbol} that are found in the headers
        """
        genes_symbols = {}
        total = 0
        found = 0
        with open(self.file_name, "r") as f:
            for line in f:
                total += 1
                if "gene_symbol:" not in line:
                    continue
                found += 1
                transcript_id = line.split(" ")[0].replace(">", "")
                gene_symbol = line.split("gene_symbol:")[1].split(" ")[0]
                genes_symbols[transcript_id] = gene_symbol

        print ("%d Genes symbols found (%d uniq) in %d reads" %
               (found, len(set(genes_symbols)), total))
        return genes_symbols

    def transcriptome_gene_description(self):
        """
        return dict {transcript_id: gene_description} that are found in the headers
        """
        genes_descs = {}
        total = 0
        found = 0
        with open(self.file_name, "r") as f:
            for line in f:
                total += 1
                if "description:" not in line:
                    continue
                found += 1
                transcript_id = line.split(" ")[0].replace(">", "")
                gene_desc = line.split("description:")[1]
                if "[" in gene_desc:
                    gene_desc = gene_desc.split("[")[0]

                genes_descs[transcript_id] = gene_desc

        print ("%d Genes description found (%d uniq) in %d reads" %
               (found, len(set(genes_descs)), total))
        return genes_descs


class OrthoDB:
    def __init__(self):
        pass

    def odb_genes_info(self, path, tax_id, info):
        """
        File: odb*genes.tab
        info = [tax_id,prot_seq_id,uniprot_id,ensemble_gene_name,ncbi_gid,description]
        return : dict {uq_ortho_gene_id : info}
        """

        info_location = {
            "tax_id": 1,
            "prot_seq_id": 2,
            "uniprot_id": 3,
            "ensemble_gene_name": 4,
            "ncbi_gid": 5,
            "description": 6
        }

        if info not in info_location.keys():
            print (
                "Please select from [uq_ortho_gene_id,tax_id,prot_seq_id,uniprot_id,ensemble_gene_name,ncbi_gid,description]")
            return 0

        result = {}

        total = 0
        found = 0

        with open(path, "r") as f:
            next(f)  # skip header
            for line in f:
                line = line.replace("\n", "").split("\t")
                if str(tax_id) not in line[1]:
                    continue

                uq_ortho_gene_id = line[0]
                feature_value = line[info_location[info]]
                total += 1
                if len(feature_value) > 2:
                    found += 1
                    result[uq_ortho_gene_id] = feature_value

        print ("(Tax_ID %d) %d %s found (%d uniq) in %d record" %
               (tax_id, found, info, len(set(result.values())), total))
        return result
