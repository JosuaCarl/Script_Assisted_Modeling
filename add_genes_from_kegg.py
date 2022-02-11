import re
import sys
import os
import cobra
from tqdm import tqdm
import pandas as pd
import memote
import helper_functions as hf
from bioservices.kegg import KEGG
import gffpandas.gffpandas as gffpd

'''
Usage: add_genes_from_kegg.py <path_input-file> <path_output-file> <path_tsv-file> <path_tsv-file-mismatches>  <path-GFF File> 
'''


def main(args):
    # console access
    if len(args) != 8:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    tsv_file_current = args[3]
    tsv_file_not_added = args[4]
    tsv_file_missing = args[5]
    gff_file = args[6]
    memote_report = args[7]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # Read in files
    model = cobra.io.sbml.read_sbml_model(infile)
    df = gffpd.read_gff3(gff_file)
    df_attr = df.attributes_to_columns()
    kegg = KEGG()

    # find organism
    req = kegg.lookfor_organism('finegoldia magna')[0].split(' ')
    entry = req[0]  # 'T00661'
    org_id = req[1]  # 'fma'
    sbo_nr = "SBO:0000243"

    # Kegg genes extraction
    genes = kegg.list(org_id).split("\n")[:-1]
    genome_dict = {g.split("\t")[0].replace("fma:", ""): g.split("\t")[1] for g in genes}

    # Data collection
    genes_current = pd.DataFrame(columns=["id", "locus_tag", "new_locus_tag", "keyword_search", "EC", "new/old"])
    genes_not_added = pd.DataFrame(columns=["id", "locus_tag", "new_locus_tag", "keyword_search", "EC", "annotations"])
    genes_missing_refseq = pd.DataFrame(
        columns=["id", "locus_tag", "new_locus_tag", "keyword_search", "EC", "annotations"])

    # -- Cobra model annotation
    for locus_tag, name in tqdm(list(genome_dict.items())):

        # Get info on gene from KEGG
        gene_dict = kegg.parse(kegg.get(org_id + ":" + locus_tag))

        # checking for keywords and EC numbers
        ec_matches = []
        keywords_matches = []
        if "ORTHOLOGY" in gene_dict:
            for go in gene_dict["ORTHOLOGY"].values():
                ec_matches = re.findall(r"(EC:{1}(?:\d+\.){3}(?:\d+){1})", go)
                keywords_matches = re.findall("DNA|RNA|histidinkinase|hypothetical|signal|putative", go, re.IGNORECASE)
        elif "NAME" in gene_dict:
            for go in gene_dict["NAME"]:
                ec_matches = re.findall(r"(EC:{1}(?:\d+\.){3}(?:\d+){1})", go)
                keywords_matches = re.findall("DNA|RNA|histidinkinase|hypothetical|signal|putative", go, re.IGNORECASE)

        # looks for DBLINKS in KEGG
        if "DBLINKS" in gene_dict:
            db_links = gene_dict["DBLINKS"]
            annotations = {"kegg.genes": f"{org_id}:{locus_tag}", "sbo": sbo_nr,
                           "ncbiprotein": db_links["NCBI-ProteinID"], "uniprot": db_links["UniProt"]}
        else:
            annotations = {"kegg.genes": f"{org_id}:{locus_tag}", "sbo": sbo_nr}

        # Get info on gene from GFF File
        locus_match = df_attr.loc[df_attr["old_locus_tag"] == locus_tag]
        if locus_match.empty:
            genes_missing_refseq.loc[len(genes_current.index)] = ["", locus_tag, "",
                                                                  keywords_matches, ec_matches, annotations]
            continue

        index = locus_match.index[0]
        # barrier for non-protein coding genes
        if df_attr.loc[index, "type"] != "gene" or df_attr.loc[index + 1, "type"] != "CDS":
            continue

        name = df_attr.loc[index + 1, "Name"]
        note = df_attr.loc[index + 1, "Note"]
        id_sbml = name[:-2] + "_" + name[-1]
        new_locus_tag = df_attr.loc[index, "locus_tag"]

        annotations = hf.dict_add_overlap_to_list(annotations, {"refseq": id_sbml[:-2] + "." + id_sbml[-1]})

        # check for critical keywords and occurring ec_codes
        if not ec_matches or keywords_matches:
            genes_not_added.loc[len(genes_current.index)] = [id_sbml, locus_tag, new_locus_tag,
                                                             keywords_matches, ec_matches, annotations]
            continue

        # Annotating old genes and adding new genes from KEGG
        if model.genes.has_id(id_sbml):
            genes_current.loc[len(genes_current.index)] = [id_sbml, locus_tag, new_locus_tag, [], [], "old"]
            model.genes.get_by_id(id_sbml).annotation = \
                hf.dict_add_overlap_to_list(model.genes.get_by_id(id_sbml).annotation, annotations)
            if note is None:
                model.genes.get_by_id(id_sbml).notes.update({"locus tag:": new_locus_tag})
            else:
                model.genes.get_by_id(id_sbml).notes.update({"locus tag:": new_locus_tag, "NCBI note:": note})
            model.genes.get_by_id(id_sbml).name = name

        # Add new gene
        else:
            genes_current.loc[len(genes_current.index)] = [id_sbml, locus_tag, new_locus_tag,
                                                           keywords_matches, ec_matches, "new"]
            gene = cobra.Gene(id=id_sbml, name=f"G_{id_sbml}", functional=True)
            model.genes.add(gene)
            model.genes.get_by_id(id_sbml).annotation = \
                hf.dict_add_overlap_to_list(model.genes.get_by_id(id_sbml).annotation, annotations)
            if note is None:
                model.genes.get_by_id(id_sbml).notes.update({"locus tag:": new_locus_tag})
            else:
                model.genes.get_by_id(id_sbml).notes.update({"locus tag:": new_locus_tag, "NCBI note:": note})
            model.genes.get_by_id(id_sbml).name = name

    # Export model
    cobra.io.sbml.write_sbml_model(model, outfile)
    genes_current.to_csv(tsv_file_current, sep="\t")
    genes_not_added.to_csv(tsv_file_not_added, sep="\t")
    genes_missing_refseq.to_csv(tsv_file_missing, sep="\t")

    # Make memote report
    result = memote.test_model(model, results=True, skip=["test_find_metabolites_not_produced_with_open_bounds"])
    report = memote.snapshot_report(result[1], config=None, html=True)
    with open(memote_report, "w") as handle:
        handle.write(report)


if __name__ == '__main__':
    main(sys.argv)
