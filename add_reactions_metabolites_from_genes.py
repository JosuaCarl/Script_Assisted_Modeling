import sys
import os
import cobra
from tqdm import tqdm
import pandas as pd
import re
import helper_functions as hf
import memote
from bioservices.kegg import KEGG

'''
Usage: add_reactions_metabolites_from_genes.py <path_input_sbml-file> <path_output_sbml-file>
<path_output_tsv-file_mismatches_bigg> <path_output_tsv-file_mismatches_locus-tags>
<path_memote-report>
'''


def main(args):
    # console access
    if len(args) != 6:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    outfile_tsv_bigg = args[3]
    outfile_tsv_lt = args[4]
    memote_report = args[5]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # read model
    model = cobra.io.read_sbml_model(infile)
    kegg = KEGG()

    # read in Bigg Reactions DB
    bigg_db = pd.read_csv("Databases/BiGG/bigg_models_reactions.tsv", sep='\t').fillna("")

    # Saves all reactions with no correspondence in BiGG
    mismatches_bigg_reacs = pd.DataFrame(
        columns=["bigg_id", "kegg_id", "coresponding_enzyme", "corresponding_locus_tag"])
    mismatches_locus_tags = pd.DataFrame(columns=["gene_id", "gene_name"])

    # -- Iterates through all genes
    for i in tqdm(range(len(model.genes))):

        # Checks for locus tag
        if "kegg.genes" in model.genes[i].annotation:
            locus_tag = model.genes[i].annotation["kegg.genes"]
        else:
            mismatches_locus_tags.loc[len(mismatches_locus_tags.index)] = [model.genes[i].id, model.genes[i].name]
            continue

        # Extracts all locus tags
        if not isinstance(locus_tag, list):
            locus_tag = [locus_tag]
        for lt in locus_tag:
            gene_dict = kegg.parse(kegg.get(lt))

            # Extracts all EC numbers for enzymes out of orthology
            ec_matches = []
            if "ORTHOLOGY" not in gene_dict:
                continue
            for go in gene_dict["ORTHOLOGY"].values():
                ec_matches = re.findall(r"(EC:{1}(?:\d+\.){3}(?:\d+){1})", go)

            # Extracts all reactions of the enzyme and compares them to BiGG DB
            enzyme_dict = dict()
            bigg_queries = dict()
            for ec_nr in ec_matches:
                enzyme_dict = kegg.parse(kegg.get(ec_nr.split(":")[1]))     # nr (exclusively) lead to enzyme
                bigg_queries[ec_nr] = bigg_db["database_links"].str\
                    .contains(f"EC Number: http://identifiers.org/ec-code/{ec_nr};", regex=False)

                if "ALL_REAC" in enzyme_dict:
                    for reac_nr in enzyme_dict["ALL_REAC"]:
                        bigg_queries[reac_nr] = bigg_db["database_links"].str\
                            .contains(f"KEGG Reaction: http://identifiers.org/kegg.reaction/{reac_nr};",
                                      regex=False)
                else:
                    continue

            # Iterate through all possible match queries
            has_bigg_entry = False
            for kegg_id, bigg_query in bigg_queries.items():
                if bigg_query.any():
                    has_bigg_entry = True
                    bigg_entries = bigg_db.loc[bigg_query]

                    # Examine reactions from BiGG
                    for j, row in bigg_entries.iterrows():

                        if not model.reactions.has_id(row["bigg_id"]):
                            reaction = cobra.Reaction(row["bigg_id"],
                                                      name=row["name"],
                                                      subsystem="")
                            reaction.gene_reaction_rule = model.genes[i].id

                            compounds = re.split(" <-> | -> | <- ", row["reaction_string"])
                            reactants = re.split(" [+] ", compounds[0])
                            products = re.split(" [+] ", compounds[1])

                            reac_metab = dict()
                            for reactant in reactants:
                                reactant = re.split(r"(\d+\.\d+) +", reactant)
                                if len(reactant) == 3:                                  # re.split keeps groups and ""
                                    reac_metab[reactant[2]] = -float(reactant[1])
                                else:
                                    reac_metab[reactant[0]] = -1.0
                            for product in products:
                                product = re.split(r"(\d+\.\d+) +", product)
                                if len(product) == 3:
                                    reac_metab[product[2]] = float(product[1])
                                else:
                                    reac_metab[product[0]] = 1.0

                            # add missing metabolites
                            wrong_compartment = False
                            for k in reac_metab.keys():
                                if k in model.metabolites:
                                    continue
                                m = hf.bigg_request(k[:-2])
                                if not m['charges']:
                                    m['charges'] = [0]

                                # check viable compartments (_c, _e, _p)
                                if re.search(r"_[cep]", k[-2:]) is None:
                                    wrong_compartment = True
                                    break

                                metabolite = cobra.Metabolite(k,
                                                              formula=m["formulae"][0],
                                                              charge=m["charges"][0],
                                                              name=m["name"],
                                                              compartment="C" + k[-2:])
                                model.add_metabolites(metabolite)

                            if wrong_compartment:
                                continue
                            model.add_reactions([reaction])
                            model.reactions.get_by_id(row["bigg_id"]).add_metabolites(reac_metab)

                        if re.search(r"((?:\d+\.){3}(?:\d+){1})", kegg_id) is None:
                            model.reactions.get_by_id(row["bigg_id"]).annotation = hf\
                                .dict_add_overlap_to_list(model.reactions.get_by_id(row["bigg_id"]).annotation,
                                                          {"kegg.reaction": kegg_id})
                        else:
                            model.reactions.get_by_id(row["bigg_id"]).annotation = hf \
                                .dict_add_overlap_to_list(model.reactions.get_by_id(row["bigg_id"]).annotation,
                                                          {"ec-code": kegg_id})
                else:
                    continue

            if not has_bigg_entry:
                reacs = enzyme_dict["ALL_REAC"] if "ALL_REAC" in enzyme_dict else []
                mismatches_bigg_reacs.loc[len(mismatches_bigg_reacs.index)] = ["", reacs,
                                                                               ec_matches, locus_tag]

    # Export model
    cobra.io.write_sbml_model(model, outfile)

    # Export mismatches to tsv
    mismatches_bigg_reacs.to_csv(outfile_tsv_bigg, sep="\t")
    mismatches_locus_tags.to_csv(outfile_tsv_lt, sep="\t")

    # Make memote report
    result = memote.test_model(model, results=True)  #, skip=["test_find_metabolites_not_produced_with_open_bounds"])
    report = memote.snapshot_report(result[1], config=None, html=True)
    with open(memote_report, "w") as handle:
        handle.write(report)


if __name__ == '__main__':
    main(sys.argv)
