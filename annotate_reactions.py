import sys
import os
import cobra.io
import libsbml
from tqdm import tqdm
import pandas as pd
import re
import memote
from bioservices.kegg import KEGG
import helper_functions as hf

'''
Usage: annotate_reactions.py <path_input-file> <outfile-csv> <program_name> <program_version> 
'''


def main(args):
    # console access
    if len(args) != 5:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    outfile_missing_bigg = args[3]
    memote_report = args[4]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # create Readers and Writers
    reader = libsbml.SBMLReader()
    writer = libsbml.SBMLWriter()

    # Read SBML File
    doc = reader.readSBML(infile)
    model = doc.getModel()

    # Knowledge base preparation
    bigg_db = pd.read_csv("Databases/BiGG/bigg_models_reactions.tsv", sep='\t').fillna("")

    seed_db = pd.read_csv("Databases/SEED/reactions.tsv", header=0, sep="\t")
    seed_db.fillna("", inplace=True)

    kegg = KEGG()
    req = kegg.lookfor_organism('finegoldia magna')[0].split(' ')
    entry = req[0]  # "T00661"
    org_code = req[1]  # 'fma'

    num_reac = model.getNumReactions()

    missing_bigg = pd.DataFrame(columns=["id", "name"])

    # BiGG
    for i in tqdm(range(num_reac)):
        reac_id = str(model.getReaction(i).getId())

        bigg_id = re.sub("^R_", "", reac_id)
        try:
            bigg_entry_idx = bigg_db.loc[bigg_db["bigg_id"] == bigg_id].index[0]
        except IndexError:
            missing_bigg.loc[len(missing_bigg.index)] = [bigg_id, model.getReaction(i).getName()]
            continue

        bigg_dblnks = bigg_db.loc[bigg_entry_idx, "database_links"].split(";")
        for db_lnk in bigg_dblnks:
            if db_lnk != "":
                lnk = db_lnk.split(": ")[1]
                model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, reac_id)

    # Export model
    doc.setModel(model)
    writer.writeSBML(doc, outfile)

    # Export tsv
    missing_bigg.to_csv(outfile_missing_bigg, sep="\t")

    # Cobra import for convenience
    model_cobra = cobra.io.read_sbml_model(outfile)

    # SEED
    for i in tqdm(range(num_reac)):
        reac_id = str(model.getReaction(i).getId())

        if "seed.reaction" in model_cobra.reactions[i].annotation:
            seed_ids = model_cobra.reactions[i].annotation["seed.reaction"]
            if not isinstance(seed_ids, list):
                seed_ids = [seed_ids]
            for seed_id in seed_ids:
                match = seed_db.loc[seed_db['id'] == seed_id]
                lnk = f"https://identifiers.org/seed.reaction/{seed_id}"
                model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, reac_id)
                for idx, row in match.iterrows():
                    for ec in row["ec_numbers"].split("|"):
                        if ec not in [None, ""]:
                            lnk = f"https://identifiers.org/ec-code/{ec}"
                            model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, reac_id)
                    for alias in row['aliases'].split("|"):
                        if alias.startswith("KEGG: "):
                            lnk = f"https://identifiers.org/kegg.reaction/{alias.split(': ')[1]}"
                            model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, reac_id)

        else:
            for db, db_ids in model_cobra.reactions[i].annotation.items():
                if not isinstance(db_ids, list):
                    db_ids = [db_ids]
                for db_id in db_ids:
                    matches = seed_db['aliases'].str.split("|") == db + ": " + db_id
                    if matches.any():
                        for idx, row in seed_db.loc[matches].iterrows():
                            lnk = f"https://identifiers.org/seed.reaction/{row['id']}"
                            model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, reac_id)
                            for ec in row["ec_numbers"].split("|"):
                                if ec not in [None, ""]:
                                    lnk = f"https://identifiers.org/ec-code/{ec}"
                                    model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, reac_id)
                            for alias in row['aliases'].split("|"):
                                if alias.startswith("KEGG: "):
                                    lnk = f"https://identifiers.org/kegg.reaction/{alias.split(': ')[1]}"
                                    model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, reac_id)

    # Export model
    doc.setModel(model)
    writer.writeSBML(doc, outfile)

    # Export tsv
    missing_bigg.to_csv(outfile_missing_bigg, sep="\t")

    # Read in model for memote
    model = cobra.io.read_sbml_model(outfile)

    # Make memote report
    result = memote.test_model(model, results=True)  #, skip=["test_find_metabolites_not_produced_with_open_bounds"])
    report = memote.snapshot_report(result[1], config=None, html=True)
    with open(memote_report, "w") as handle:
        handle.write(report)


if __name__ == '__main__':
    main(sys.argv)
