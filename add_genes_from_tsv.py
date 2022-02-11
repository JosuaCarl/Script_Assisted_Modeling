import sys
import os
import cobra
import pandas as pd
import helper_functions as hf

'''
Usage: add_genes_from_kegg.py <path_input-file> <path_output-file> <path-tsv-file> <path-GFF File> <from-file>
'''


def main(args):
    # console access
    if len(args) < 6:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    tsv_file = args[3]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # Read in files
    model = cobra.io.sbml.read_sbml_model(infile)
    genes_mismatches = pd.read_csv(tsv_file, sep="\t")

    for i, row in genes_mismatches.iterrows():
        gene = cobra.Gene(id=row["id"], name=f"G_{row['id']}", functional=True)  # What is functional ? Coding ?
        model.genes.add(gene)
        model.genes.get_by_id(row["locus_tag"]).annotation = \
            hf.dict_add_overlap_to_list(model.genes.get_by_id(row["locus_tag"]).annotation, row["annotations"])

    # Export model
    cobra.io.sbml.write_sbml_model(model, outfile)


if __name__ == '__main__':
    main(sys.argv)
