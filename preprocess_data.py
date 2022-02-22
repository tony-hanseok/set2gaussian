import argparse
import os
from collections import defaultdict
import pandas as pd


def get_protein_gene_mapping(file_path, ensembl2symbols):
    with open(file_path, "r") as f:
        lines = f.readlines()
    protein2gene = {}
    for line in lines[1:]:
        protein_id, preferred, *_ = line.split("\t")

        if preferred.startswith("ENSG"):
            try:
                preferred = ensembl2symbols[preferred]
            except KeyError:
                continue

        protein2gene[protein_id] = preferred.lower()
    return protein2gene


def preprocess_string(string_db_path, ensembl2symbols, threshold):
    protein2gene = get_protein_gene_mapping(
        file_path=os.path.join(string_db_path, "9606.protein.info.v11.5.txt"),
        ensembl2symbols=ensembl2symbols,
    )

    links = []
    genes = []
    with open(os.path.join(string_db_path, "9606.protein.links.v11.5.txt"), "r") as f:
        lines = f.readlines()
    for line in lines[1:]:
        try:
            prot0, prot1, score = line.strip().split()
        except ValueError:
            continue

        score = float(score) / 1000
        if score < threshold:
            continue

        try:
            links.append(
                (protein2gene[prot0], protein2gene[prot1], float(score) / 1000)
            )
            genes.append(protein2gene[prot0])
            genes.append(protein2gene[prot1])
        except KeyError:
            continue
    return links, set(genes)


def preprocess_cpdb(cpdb_file):
    node_set = []
    pathway_desc = []
    genes = []
    with open(cpdb_file, "r") as f:
        lines = f.readlines()

    source_id_dict = defaultdict(lambda: [])
    for line in lines[1:]:
        pathway, external_id, source, gene_symbols = line.strip().split("\t")
        for gene_symbol in gene_symbols.split(","):
            if external_id is None or external_id == "None":
                source_id_dict[source].append(pathway)
                external_id = f"{source}_{len(source_id_dict[source])}"

            gene_symbol = gene_symbol.lower()
            node_set.append((external_id, gene_symbol))
            genes.append(gene_symbol)
        pathway_desc.append((external_id, pathway))
    return node_set, pathway_desc, set(genes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--threshold", type=float, default=0.7)
    parser.add_argument(
        "--prefix", type=str, default="/home/tony/gits/set2gaussian/data"
    )
    parser.add_argument(
        "--gene_table_path",
        type=str,
        default="/mnt/aitrics_ext/ext01/shared/BI/db_integration_v1/tables/gene_table.tsv",
    )
    args = parser.parse_args()

    # 1. Get gene table
    gene_table = pd.read_csv(args.gene_table_path, sep="\t", index_col="hgnc_id")
    ensembl2symbols = gene_table[["ensembl_gene_id", "symbol"]].dropna()
    ensembl2symbols = {
        row["ensembl_gene_id"]: row["symbol"]
        for _, row in ensembl2symbols.dropna().iterrows()
    }

    # 2. Get PPI network and related genes (converting string id to gene symbol)
    string_db_path = "/mnt/aitrics_ext/ext01/shared/BI/public/STRING"
    PPI_net, string_genes = preprocess_string(
        string_db_path=string_db_path,
        ensembl2symbols=ensembl2symbols,
        threshold=args.threshold,
    )

    # 3. Get pathways from CPDB
    cpdb_file = os.path.join(args.prefix, "CPDB_pathways_genes.tab")
    node_set, pathway_desc, cpdb_genes = preprocess_cpdb(cpdb_file=cpdb_file)
    common_genes = string_genes.intersection(cpdb_genes)

    print(f"Threshold: {args.threshold}")
    print(f"STRING only: {len(string_genes -cpdb_genes)}")
    print(f"CPDB only: {len(cpdb_genes - string_genes)}")
    print(f"Common: {len(common_genes)}")

    # 4. Write
    with open(os.path.join(args.prefix, str(args.threshold), f"network.txt"), "w") as f:
        invalid_gene_links = 0
        valid_gene_links = 0
        for gene0, gene1, score in PPI_net:
            if gene0 not in common_genes or gene1 not in common_genes:
                invalid_gene_links += 1
                continue

            f.write(f"{gene0}\t{gene1}\t{score}\n")
            valid_gene_links += 1
        print(f"Valid\t{valid_gene_links}\tInvalid\t{invalid_gene_links}")

    with open(
        os.path.join(args.prefix, str(args.threshold), f"cpdb_node_set.txt"), "w"
    ) as f:
        invalid_gene_pathways = 0
        valid_gene_pathways = 0
        for pathway_id, gene_symbol in node_set:
            if gene_symbol not in common_genes:
                invalid_gene_pathways += 1
                continue

            f.write(f"{pathway_id}\t{gene_symbol}\n")
            valid_gene_pathways += 1
        print(f"Valid\t{valid_gene_pathways}\tInvalid\t{invalid_gene_pathways}")

    with open(os.path.join(args.prefix, str(args.threshold), f"cpdb_desc.txt"), "w") as f:
        for pathway_id, desc in pathway_desc:
            f.write(f"{pathway_id}\t{desc}\n")
