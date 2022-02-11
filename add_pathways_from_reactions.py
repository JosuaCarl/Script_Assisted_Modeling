import sys
import os
import cobra
import pandas as pd
from tqdm import tqdm
import memote
import libsbml
import helper_functions as hf
from bioservices.kegg import KEGG

'''
Usage: add_genes_from_kegg.py <path_input-file> <path_output-file>
'''


def main(args):
    # console access
    if len(args) != 5:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    outfile_tsv = args[3]
    memote_report = args[4]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # create Readers and Writers
    reader = libsbml.SBMLReader()
    writer = libsbml.SBMLWriter()
    kegg = KEGG()

    # Read SBML File
    doc = reader.readSBML(infile)
    model = doc.getModel()

    # find organism
    req = kegg.lookfor_organism('finegoldia magna')[0].split(' ')
    org_id = req[1]  # 'fma'

    model_cobra = cobra.io.read_sbml_model(infile)

    changes_pathways = pd.DataFrame(columns=["pos", "gene_id", "pathway"])
    start = 0

    # accessing previous progress
    if os.path.exists(outfile_tsv):
        changes_pathways_old = pd.read_csv(outfile_tsv, sep="\t", index_col=0)
        start = changes_pathways_old.loc[len(changes_pathways_old.index)-1, "pos"]

    # -- Pathway annotation via KEGG
    reac_num = model.getNumReactions()
    for i in tqdm(range(start, reac_num)):
        pathways = dict()
        if "kegg.reaction" in model_cobra.reactions[i].annotation:
            reac_ids = model_cobra.reactions[i].annotation["kegg.reaction"]
            if not isinstance(reac_ids, list):
                reac_ids = [reac_ids]
            for reac_id in reac_ids:
                reaction = kegg.parse(kegg.get(reac_id))
                if "PATHWAY" in reaction:
                    pathway = reaction["PATHWAY"]
                    pathways.update(pathway)

        genes = model_cobra.reactions[i].genes
        for gene in genes:
            if "kegg.genes" in gene.annotation:
                locus_tags = gene.annotation["kegg.genes"]
                if not isinstance(locus_tags, list):
                    locus_tags = [locus_tags]
                for locus_tag in locus_tags:
                    locus_tag = locus_tag.split(":")[1]
                    pathway = kegg.get_pathway_by_gene(locus_tag, org_id)
                    if pathway is not None:
                        pathways.update(pathway)

        if pathways is not None:
            changes_pathways.loc[len(changes_pathways.index)] = [i, model_cobra.reactions[i].id, pathways]
            for pathway_id_kegg in pathways.keys():
                lnk = f"https://identifiers.org/kegg.pathway/{pathway_id_kegg}"
                model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_OCCURS_IN, i)

        if i % 100 == 99:
            # Export model
            doc.setModel(model)
            writer.writeSBML(doc, outfile)

            # Export progress
            if os.path.exists(outfile_tsv):
                changes_pathways_old = pd.read_csv(outfile_tsv, sep="\t", index_col=0)
                changes_pathways = pd.concat([changes_pathways_old, changes_pathways])
                changes_pathways.reset_index(drop=True, inplace=True)
                changes_pathways_old = None
            changes_pathways.to_csv(outfile_tsv, sep="\t")
            changes_pathways = pd.DataFrame(columns=["pos", "gene_id", "pathway"])

    # Export model
    doc.setModel(model)
    writer.writeSBML(doc, outfile)

    # Export changes
    if os.path.exists(outfile_tsv):
        changes_pathways_old = pd.read_csv(outfile_tsv, sep="\t", index_col=0)
        changes_pathways = pd.concat([changes_pathways_old, changes_pathways])
        changes_pathways.reset_index(drop=True, inplace=True)
        changes_pathways_old = None
    changes_pathways.to_csv(outfile_tsv, sep="\t")

    # Import model as cobra model for memote
    model = cobra.io.read_sbml_model(outfile)

    # Make memote report
    result = memote.test_model(model, results=True)
    report = memote.snapshot_report(result[1], config=None, html=True)
    with open(memote_report, "w") as handle:
        handle.write(report)


if __name__ == '__main__':
    main(sys.argv)
