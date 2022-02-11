import sys
import os
import libsbml
from tqdm import tqdm
import pandas as pd
import helper_functions as hf

'''
Usage: amend_formulas.py <path_input-file> <path_output-file> <outfile-csv> <program_name> <program_version> 
<tolerate_charge_hydrogen_balancing> : -chBal, if +1 charge should correspond to +1 H-atom
Takes formulas from the notes field and fbc-plugin, if none are found, BiGG-DB is searched for a formula. 
If multiple or no possibilities are given in BiGG, a csv-formatted table with these metabolites is returned.
'''


def main(args):
    # console access
    if len(args) < 6:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    outfile_csv = args[3]
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

    # Check for errors
    if doc.getNumErrors() > 0:
        # ++ Add constant and Boundary Condition (=False) to the model ++
        for i in range(0, model.getNumSpecies()):
            model.getSpecies(i).setConstant(False)
            model.getSpecies(i).setBoundaryCondition(False)

        # Set newly created model and write to file
        doc.setModel(model)
        writer.setProgramName(program_name)
        writer.setProgramVersion(program_version)
        writer.writeSBML(doc, outfile)

    doc = reader.readSBML(outfile)
    print("Document errors: " + str(doc.getNumErrors()))
    doc.printErrors()

    # Use BiGG Database for formulae check, if none is given
    mismatches = pd.DataFrame(columns=["model_index", "id", "name", "formula_bigg", "formula_model"])
    num_spec = model.getNumSpecies()
    for i in tqdm(range(num_spec)):

        # check fbc plugin for formula
        if model.getSpecies(i).getPlugin('fbc').isSetChemicalFormula():
            continue

        meta_id = str(model.getSpecies(i).getMetaId())
        pruned_id = meta_id[2:-2]
        metabolite_info = hf.bigg_request(pruned_id, "metabolites")
        formulas_bigg = metabolite_info["formulae"]
        if len(formulas_bigg) == 1:
            model.getSpecies(i).getPlugin('fbc').setChemicalFormula(formulas_bigg)
            note_str = f"Changed formula from '' to {formulas_bigg[0]}. Source: BiGG"
            model = hf.add_note_species(model, note_str, meta_id)
            lnk = f"https://identifiers.org/bigg.reaction:{pruned_id}"
            model = hf.add_link_annotation_species(model, lnk, meta_id)
        else:
            mismatches.loc[len(mismatches.index)] = [i, meta_id, model.getSpecies(i).getName(),
                                                     formulas_bigg, ""]

    # Saving new model
    doc.setModel(model)
    writer.setProgramName(program_name)
    writer.setProgramVersion(program_version)
    writer.writeSBML(doc, outfile)

    # Exporting mismatches
    mismatches.to_csv(outfile_csv)


if __name__ == '__main__':
    main(sys.argv)
