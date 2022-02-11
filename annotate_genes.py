import sys
import os
import cobra
from tqdm import tqdm
import memote
import helper_functions as hf
from bioservices.kegg import KEGG
import gffpandas.gffpandas as gffpd

'''
Usage: annotate_genes.py <path_input-file> <path_output-file> <path_tsv-file> <path_tsv-file-mismatches>  <path-GFF File> 
'''


def main(args):
    # console access
    if len(args) != 5:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    gff_file = args[3]
    memote_report = args[4]

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

    # -- Cobra model annotation
    for i in tqdm(range(len(model.genes))):
        annotations = {"sbo": sbo_nr}
        id_sbml = model.genes[i].id
        refseq = id_sbml[:-2] + "." + id_sbml[-1]
        gff_query = df_attr["Name"] == refseq

        if gff_query.any():
            matches = df_attr[gff_query]
            for j, row in matches.iterrows():
                if df_attr.loc[j - 1, "type"] == "gene" and df_attr.loc[j, "type"] == "CDS":
                    locus_tag = df_attr.loc[j - 1, "old_locus_tag"]
                    new_locus_tag = df_attr.loc[j - 1, "locus_tag"]
                    note = df_attr.loc[j, "Note"]
                    name = df_attr.loc[j, "Name"]

                    annotations = hf.dict_add_overlap_to_list(annotations, {"kegg.genes": f"{org_id}:{locus_tag}",
                                                                            "refseq": refseq})

                    model.genes[i].annotation = \
                        hf.dict_add_overlap_to_list(model.genes[i].annotation, annotations)
                    if note is None:
                        model.genes[i].notes.update({"locus tag:": new_locus_tag})
                    else:
                        model.genes[i].notes.update({"locus tag:": new_locus_tag, "NCBI note:": note})
                    model.genes[i].name = name

    # Export model
    cobra.io.sbml.write_sbml_model(model, outfile)

    # Make memote report
    result = memote.test_model(model, results=True, skip=["test_find_metabolites_not_produced_with_open_bounds"])
    report = memote.snapshot_report(result[1], config=None, html=True)
    with open(memote_report, "w") as handle:
        handle.write(report)


if __name__ == '__main__':
    main(sys.argv)
