import sys
import os
import libsbml
from tqdm import tqdm
import pandas as pd
from itertools import product
from bioservices.kegg import KEGG
from requests.exceptions import HTTPError, RequestException
import helper_functions as hf

'''
Usage: check+annotate_metabolites.py <path_input-file> <outfile-csv> <program_name> <program_version> 
<tolerate_charge_hydrogen_balancing> : -chBal, if +1 charge should correspond to +1 H-atom
Takes formulas from the notes field and fbc-plugin, if none are found, BiGG-DB is searched for a formula. 
If multiple or no possibilities are given in BiGG, a csv-formatted table with these metabolites is returned.
'''


def main(args):
    # console access
    if len(args) < 3:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile_mismatches = args[2]
    outfile_formula_search = args[3]

    tolerate_ch_h_bal = "-chBal" in args

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # create Readers and Writers
    reader = libsbml.SBMLReader()

    # Read SBML File
    doc = reader.readSBML(infile)
    model = doc.getModel()

    # Knowledge base preparation
    # bigg_db = pd.read_csv("Databases/BiGG/bigg_models_metabolites.tsv", sep='\t')

    mnx_db = pd.read_csv("Databases/MetaNetX/chem_prop.tsv", header=351, sep='\t')
    mnx_db.rename(columns={'#ID': 'id'}, inplace=True)
    mnx_db.fillna("", inplace=True)

    seed_db = pd.read_csv("Databases/SEED/compounds.tsv", header=0, sep="\t")
    seed_db.fillna("", inplace=True)

    kegg = KEGG()
    req = kegg.lookfor_organism('finegoldia magna')[0].split(' ')
    entry = req[0]  # "T00661"
    org_code = req[1]  # 'fma'

    # -------- formula check against knowledge bases ---------
    start = 0
    if os.path.exists(outfile_mismatches):
        mismatches_old = pd.read_csv(outfile_mismatches, sep="\t", index_col=0)
        start = max(start, int(mismatches_old.tail(1)["model_index"][1]) + 1)
        mismatches_old = None
    if os.path.exists(outfile_formula_search):
        formula_searches_old = pd.read_csv(outfile_formula_search, sep="\t", index_col=0)
        start = max(start, int(formula_searches_old.tail(1)["model_index"][1]) + 1)
        mismatches_old = None

    mismatches = pd.DataFrame(columns=["model_index", "name", "spec_id", "ids_biocyc", "ids_metanetx", "ids_seed",
                                       "formula_bigg", "formula_biocyc", "formula_metanetx", "formula_seed",
                                       "formula_model",
                                       "charge_bigg", "charge_biocyc", "charge_metanetx", "charge_seed",
                                       "charge_model", "matching_db"])
    formula_search = pd.DataFrame(columns=["model_index", "name", "spec_id", "formula_model",
                                           "ids_biocyc", "ids_mnx", "ids_seed", "ids_kegg"])
    form_comp = [False, False, False, False, False]
    num_spec = model.getNumSpecies()
    for i in tqdm(range(start, num_spec)):
        # --------- Knowledge collection ---------
        spec_id = str(model.getSpecies(i).getId())

        try:
            # Check for formula in model
            formula_model = []
            charge_model = []
            name_model = ""
            if model.getSpecies(i).getPlugin('fbc').isSetChemicalFormula():
                formula_model = str(model.getSpecies(i).getPlugin('fbc').getChemicalFormula())
            if model.getSpecies(i).getPlugin('fbc').isSetCharge():
                charge_model = model.getSpecies(i).getPlugin('fbc').getCharge()
            name_model = model.getSpecies(i).getName()

            # BiGG              formulas  - commented out: extraction from tsv (contains no formula and charge)
            # bigg_query = bigg_db.loc[bigg_db['bigg_id'] == pruned_id]
            bigg_query = hf.bigg_request(spec_id[2:-2], "metabolites")
            formulas_bigg = bigg_query["formulae"]
            charges_bigg = bigg_query["charges"]
            try:
                inchikey = bigg_query["database_links"]["InChi Key"]["id"]
            except KeyError:
                inchikey = False

            # Biocyc
            formulas_biocyc = []
            charges_biocyc = []
            ids_biocyc = []
            biocyc_req = hf.biocyc_request("GCF_000010185", "BIGG", spec_id[2:-2])
            if not biocyc_req[0]["STATUS"] == 1:
                biocyc_req = hf.biocyc_get_from_formula("GCF_000010185", formula_model)
                form_comp[1] = True
            if biocyc_req[0]["STATUS"] == 1:
                for res in biocyc_req[0]["RESULTS"]:
                    ids_biocyc.append(res["ID"])
                    try:
                        biocyc_tree = hf.biocyc_get(id_org="meta", id_biocyc=res["ID"], detail="low")
                        charges_biocyc.append(
                            int(biocyc_tree["ptools-xml"]["Compound"]["cml"]["molecule"]["@formalCharge"]))
                        formulas_biocyc_str = biocyc_tree["ptools-xml"]["Compound"]["cml"]["molecule"]["formula"][
                            "@concise"]
                        formulas_biocyc.append(formulas_biocyc_str.replace(" ", ""))
                    except KeyError:
                        print(spec_id + ": no simple compound.")
                    except HTTPError or RequestException:
                        print(spec_id + " failed in biocyc request.")

            # MetaNetX
            charges_mnx = []
            formulas_mnx = []
            ids_mnx = []
            if inchikey:
                mnx_query = mnx_db.loc[mnx_db['InChiKey'] == inchikey]
            elif name_model != "":
                mnx_query = mnx_db.loc[mnx_db['name'] == name_model]
            elif formula_model != "":
                mnx_query = mnx_db.loc[mnx_db['formula'] == formula_model]
                form_comp[2] = True
            else:
                mnx_query = pd.DataFrame({'formula': [], 'charge': []})
            for idx, row in mnx_query.iterrows():
                ids_mnx.append(row["id"])
                if row["formula"] != "" and row["charge"] != "":
                    formulas_mnx.append(row['formula'])
                    charges_mnx.append(row['charge'])

            # SEED
            formulas_seed = []
            charges_seed = []
            ids_seed = []
            search = seed_db['aliases'].str.contains("BiGG: " + spec_id[2:-2])
            if search.any():
                seed_query = seed_db.loc[search]
            elif inchikey:
                seed_query = seed_db.loc[seed_db['inchikey'] == inchikey]
            elif name_model != "":
                seed_query = seed_db.loc[seed_db['name'] == name_model]
            elif formula_model != "":
                seed_query = seed_db.loc[seed_db['formula'] == formula_model]
                form_comp[3] = True
            else:
                seed_query = pd.DataFrame({'formula': [], 'charge': []})
            for idx, row in seed_query.iterrows():
                ids_seed.append(row['id'])
                formulas_seed.append(row['formula'])
                charges_seed.append(int(row['charge']))

            # KEGG
            ids_kegg = []
            if formula_model != "":
                kegg_query = kegg.find("compound", formula_model, "formula").split("\n")
                form_comp[4] = True
            else:
                kegg_query = pd.DataFrame({'formula': [], 'charge': []})
            for kq in kegg_query:
                ids_kegg.append(kq.split("\t")[0])

        except Exception as e:
            print(spec_id)
            raise e

        # --------- Knowledge vs. current entry - comparison ---------
        matching_dbs = []
        formula_matching_ids = []
        if not model.getSpecies(i).getPlugin('fbc').isSetChemicalFormula():
            mismatches.loc[len(mismatches.index)] = [i, model.getSpecies(i).getName(),
                                                     spec_id, ids_biocyc, ids_mnx, ids_seed,
                                                     formulas_bigg, formulas_biocyc, formulas_mnx, formulas_seed,
                                                     formula_model,
                                                     charges_bigg, charges_biocyc, charges_mnx, charges_seed,
                                                     charge_model,
                                                     matching_dbs]
            continue

        formulas_all = [formulas_bigg, formulas_biocyc, formulas_mnx, formulas_seed]
        charges_all = [charges_bigg, charges_biocyc, charges_mnx, charges_seed]
        ids_all = [spec_id[2:-2], ids_biocyc, ids_mnx, ids_seed]
        comparisons_bool = []
        for j in range(len(formulas_all)):
            if form_comp[j]:
                formula_matching_ids.append(ids_all[j])
            else:
                if tolerate_ch_h_bal:
                    for_cha_product = list(product(formulas_all[j], charges_all[j]))
                else:
                    for_cha_product = formulas_all[j]

                for fcp in for_cha_product:
                    if tolerate_ch_h_bal and charges_all[j] and charge_model:
                        comparison = hf.compare_formulas([fcp[0], formula_model], [int(fcp[1]), charge_model])
                        comparisons_bool.append(comparison)
                    else:
                        comparison = hf.compare_formulas([fcp[0], formula_model])
                        comparisons_bool.append(comparison)

                    if comparison:
                        matching_dbs.append(j)

        # --------- Collection in table ---------
        if True not in comparisons_bool:
            mismatches.loc[len(mismatches.index)] = [i, model.getSpecies(i).getName(),
                                                     spec_id, ids_biocyc, ids_mnx, ids_seed,
                                                     formulas_all[0], formulas_all[1], formulas_all[2], formulas_all[3],
                                                     formula_model,
                                                     charges_all[0], charges_all[1], charges_all[2], charges_all[3],
                                                     charge_model,
                                                     matching_dbs]
        formula_search.loc[len(formula_search.index)] = [i, model.getSpecies(i).getName(), spec_id, formula_model,
                                                         ids_biocyc, ids_mnx, ids_seed, ids_kegg]

        # in between saves
        if i % 50 == 25:
            if os.path.exists(outfile_mismatches):
                mismatches_old = pd.read_csv(outfile_mismatches, sep="\t", index_col=0)
                mismatches = pd.concat([mismatches_old, mismatches])
                mismatches.reset_index(drop=True, inplace=True)
                mismatches_old = None
            mismatches.to_csv(outfile_mismatches, sep="\t")
            mismatches = pd.DataFrame(
                columns=["model_index", "name", "spec_id", "ids_biocyc", "ids_metanetx", "ids_seed",
                         "formula_bigg", "formula_biocyc", "formula_metanetx", "formula_seed",
                         "formula_model",
                         "charge_bigg", "charge_biocyc", "charge_metanetx", "charge_seed",
                         "charge_model", "matching_db"])

            if os.path.exists(outfile_formula_search):
                formula_search_old = pd.read_csv(outfile_formula_search, sep="\t", index_col=0)
                formula_search = pd.concat([formula_search_old, formula_search])
                formula_search.reset_index(drop=True, inplace=True)
                formula_search_old = None
            formula_search.to_csv(outfile_formula_search, sep="\t")
            formula_search = pd.DataFrame(columns=["model_index", "name", "spec_id", "formula_model",
                                                   "ids_biocyc", "ids_mnx", "ids_seed", "ids_kegg"])

    # Exporting mismatches and formula search results
    mismatches_old = pd.read_csv(outfile_mismatches, sep="\t", index_col=0)
    mismatches = pd.concat([mismatches_old, mismatches])
    mismatches.reset_index(drop=True, inplace=True)
    mismatches_old = None
    mismatches.to_csv(outfile_mismatches, sep="\t")
    formula_search_old = pd.read_csv(outfile_formula_search, sep="\t", index_col=0)
    formula_search = pd.concat([formula_search_old, formula_search])
    formula_search.reset_index(drop=True, inplace=True)
    formula_search_old = None
    formula_search.to_csv(outfile_formula_search, sep="\t")


if __name__ == '__main__':
    main(sys.argv)
