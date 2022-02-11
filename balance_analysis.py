import sys
import os
import cobra
import re
import pandas as pd
from tqdm import tqdm
import helper_functions as hf

'''
Usage: balance_analysis.py <path_input-file> <path_output>
The output is a .csv-format table with the columns: "model_index", "imbalances", "reaction_string", "frequent_compound"
The last row of the table contains a list of compounds that consistently produce the same problem
'''


def main(args):
    # console access
    if len(args) != 3:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # create Readers and Writers
    model = cobra.io.sbml.read_sbml_model(infile)

    # check mass balance:
    unbalanced_list = pd.DataFrame(columns=["model_index", "reaction_name", "imbalances", "reaction_string",
                                            "formulas", "frequent_compound"])
    compound_counter = dict()
    for i in tqdm(range(len(model.reactions))):
        reaction = model.reactions[i]
        imbalances = reaction.check_mass_balance()
        formulas = ""
        charges = ""

        # Makes human readable and interpretable lists, when imbalances exist
        if len(imbalances) > 0 and re.search("EX_|sink_|Growth", reaction.id) is None:
            reactants, operator, products = re.split("(-->|<=>|<--)", reaction.build_reaction_string())

            metabolites = list(reaction.metabolites.keys())
            for m in metabolites:
                if m.id in compound_counter.keys():
                    compound_counter[m.id].append(i)
                else:
                    compound_counter[m.id] = [i]

                formulas = formulas + m.formula + f"({m.charge}),\n"

            formulas = formulas[:-2]
            reactants, operator, products = re.split("(-->|<=>|<--)", reaction.build_reaction_string())
            reac_str = reactants + operator + "\n" + products
            unbalanced_list.loc[len(unbalanced_list.index)] = [i, reaction.id, imbalances, reac_str, formulas, ""]

    # Extract frequent compounds, that are present with the same imbalances in multiple cases
    frequent_compounds = []
    for cckey, ccvalues in compound_counter.items():

        list_imbalances = []
        for pos in ccvalues:

            imbalances = list(unbalanced_list.loc[unbalanced_list.model_index == pos, "imbalances"])[0]
            negative_imbalances = dict()
            for k, v in imbalances.items():
                negative_imbalances[k] = -v

            if imbalances in list_imbalances or negative_imbalances in list_imbalances:
                frequent_compounds.append(cckey)
                unbalanced_list.loc[unbalanced_list.model_index == pos, "frequent_compound"] = \
                    unbalanced_list.loc[unbalanced_list.model_index == pos, "frequent_compound"] + cckey + ", "

            list_imbalances.append(imbalances)

    # Make table human-readable
    frequent_compounds = hf.delete_doubles(frequent_compounds)
    fc_str = ""
    for i in range(len(frequent_compounds)):
        if i % 5 == 4:
            fc_str = fc_str + str(frequent_compounds[i]) + ",\n"
        else:
            fc_str = fc_str + str(frequent_compounds[i]) + ","

    for i in range(len(unbalanced_list.index)):
        imbalances_str = ""
        for key in unbalanced_list.loc[i, "imbalances"].keys():
            imbalances_str = imbalances_str + str(key) + ": " + str(unbalanced_list.loc[i, "imbalances"][key]) + "\n"
        imbalances_str = imbalances_str[:-1]
        unbalanced_list.loc[i, "imbalances"] = imbalances_str

    unbalanced_list.loc[len(unbalanced_list.index)] = [-1, "Recurring imbalances:",
                                                       fc_str, len(frequent_compounds), "",  ""]

    # export list
    unbalanced_list.to_csv(outfile, sep="\t")


if __name__ == '__main__':
    main(sys.argv)
