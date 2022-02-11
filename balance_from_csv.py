import sys
import os
import libsbml
from tqdm import tqdm
import pandas as pd
import helper_functions as hf

'''
Usage: amend_formulas.py <path_input-file> <path_output-file> <infile-csv> <program_name> <program_version> 
'''


def main(args):
    # console access
    if len(args) < 6:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    infile_csv = args[3]
    program_name = args[4]
    program_version = args[5]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # create Readers and Writers
    reader = libsbml.SBMLReader()
    writer = libsbml.SBMLWriter()

    # Read SBML File
    doc = reader.readSBML(infile)
    model = doc.getModel()

    # parse through table
    table = pd.read_csv(infile_csv)
    for i in tqdm(range(1, len(table["id"]))):
        meta_id = table["id"][i]
        try:
            if table["id"][i].startswith("R_"):
                if table["foundation"][i] == "SEED":
                    lnk = f"https://identifiers.org/seed.reaction:{table['db_id'][i]}"
                    model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, meta_id)
                elif table["foundation"][i] == "BiGG":
                    lnk = f"https://identifiers.org/bigg.reaction:{table['db_id'][i]}"
                    model = hf.add_link_annotation_reaction(model, lnk, libsbml.BQB_IS, meta_id)

                if table["change_type"][i] == "product":
                    new_comp = table["new"][i].split(" ")
                    comp_nr = float(new_comp[0])
                    comp = new_comp[1]
                    if comp_nr == 0:
                        model.getReaction(meta_id).removeProduct("M_" + comp)
                    else:
                        species = model.getSpecies("M_" + comp)
                        model.getReaction(meta_id).removeProduct("M_" + comp)
                        model.getReaction(meta_id).addProduct(species, comp_nr)
                    note_str = f"Changed product from {table['old'][i]} to {table['new'][i]}. Source: {table['foundation'][i]}"
                    model = hf.add_note_reaction(model, note_str, meta_id)

                if table["change_type"][i] == "reactant":
                    new_comp = table["new"][i].split(" ")
                    comp_nr = float(new_comp[0])
                    comp = new_comp[1]
                    if comp_nr == 0:
                        model.getReaction(meta_id).removeReactant("M_" + comp)
                    else:
                        species = model.getSpecies("M_" + comp)
                        model.getReaction(meta_id).removeReactant("M_" + comp)
                        model.getReaction(meta_id).addReactant(species, comp_nr)
                    note_str = f"Changed reactant from {table['old'][i]} to {table['new'][i]}. Source: {table['foundation'][i]}"
                    model = hf.add_note_reaction(model, note_str, meta_id)

                if type(table["notes"][i]) == str and table["notes"][i] != "":
                    model = hf.add_note_reaction(model, table["notes"][i], meta_id)

                if type(table["eco"][i]) == str and table["eco"][i] != "":
                    link = f"https://identifiers.org/eco/{table['eco'][i]}"
                    model = hf.add_link_annotation_species(model, link, libsbml.BQB_IS, meta_id)

            else:
                meta_id = "M_" + table["id"][i]
                if table["change_type"][i] == "charge":
                    note_str = f"Changed charge from {table['old'][i]} to {table['new'][i]}. Source: {table['foundation'][i]}"
                    model = hf.add_note_species(model, note_str, meta_id)
                    model.getSpecies(meta_id).getPlugin('fbc').setCharge(int(table["new"][i]))
                elif table["change_type"][i] == "formula":
                    note_str = f"Changed formula from {table['old'][i]} to {table['new'][i]}. Source: {table['foundation'][i]}"
                    model = hf.add_note_species(model, note_str, meta_id)
                    model.getSpecies(meta_id).getPlugin('fbc').setChemicalFormula(table["new"][i])

                if table["foundation"][i] == "SEED":
                    link = f"https://identifiers.org/seed.compound:{table['db_id'][i]}"
                    model = hf.add_link_annotation_species(model, link, libsbml.BQB_IS, meta_id)
                elif table["foundation"][i] == "BiGG":
                    link = f"https://identifiers.org/bigg.metabolite:{table['db_id'][i]}"
                    model = hf.add_link_annotation_species(model, link, libsbml.BQB_IS, meta_id)
                elif table["foundation"][i] == "MetaCyc":
                    link = f"https://identifiers.org/metacyc.compound:{table['db_id'][i]}"
                    model = hf.add_link_annotation_species(model, link, libsbml.BQB_IS, meta_id)
                elif table["foundation"][i] == "MetaNetX":
                    link = f"https://identifiers.org/metanetx.chemical:{table['db_id'][i]}"
                    model = hf.add_link_annotation_species(model, link, libsbml.BQB_IS, meta_id)
                elif table["foundation"][i] == "KEGG":
                    link = f"https://identifiers.org/kegg.compound:{table['db_id'][i]}"
                    model = hf.add_link_annotation_species(model, link, libsbml.BQB_IS, meta_id)

                if type(table["notes"][i]) == str and table["notes"][i] != "":
                    model = hf.add_note_species(model, table["notes"][i], meta_id)

                if type(table["eco"][i]) == str and table["eco"][i] != "":
                    link = f"https://identifiers.org/eco/{table['eco'][i]}"
                    model = hf.add_link_annotation_species(model, link, libsbml.BQB_IS_DESCRIBED_BY, meta_id)
        except AttributeError as ae:
            print(meta_id)
            raise ae

    # Saving new model
    doc.setModel(model)
    writer.setProgramName(program_name)
    writer.setProgramVersion(program_version)
    writer.writeSBML(doc, outfile)


if __name__ == '__main__':
    main(sys.argv)
