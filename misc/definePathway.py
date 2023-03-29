import json
import pandas as pd

def pathway_define():
    ConsensusPDB_dict = dict()
    for index, row in ConsensusPDB.iterrows():
        ConsensusPDB_dict[row["pathway"]] = {"external_id":row["external_id"],
                                             "source":row["source"],
                                             "hgnc_symbol_ids":row["hgnc_symbol_ids"]}
def index_re_gene2pathway():
    Overall_gene_set = list(set([gene for genes in ConsensusPDB["hgnc_symbol_ids"].str.split(",") for gene in genes]))

    belonged_path_dict = {gene: [row["pathway"]
                                 for _, row in ConsensusPDB.iterrows() if gene in row["hgnc_symbol_ids"].split(",")]
                          for gene in Overall_gene_set}

    with open("../dbs/indexed_ConsensusPDB.json", "w") as indexed_CPDB_fl:
        json.dump(belonged_path_dict,indexed_CPDB_fl)

    return None

ConsensusPDB = pd.read_csv('../dbs/CPDB_pathways_genes_KEGG_WikiPath.tsv',delimiter="\t")

if __name__ == '__main__':
    index_re_gene2pathway()