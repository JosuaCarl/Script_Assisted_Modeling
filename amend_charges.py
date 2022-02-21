import sys
import os
import libsbml
from tqdm import tqdm
import pandas as pd
import helper_functions as hf

'''
Usage: amend_charges.py <path_input_sbml-file> <path_output_sbml-file>
<path_out_tsv-file_mismatches>
Transfers charges from Notes to the fbc-Plugin Annotation. If none is found, BiGG-DB is used for a search. If BiGG 
contains multiple or no charges, the metabolite is returned as a csv-formatted file.
'''


def main(args):
    # console access
    if len(args) < 6:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    outfile_tsv = args[3]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    reader = libsbml.SBMLReader()
    writer = libsbml.SBMLWriter()

    # Conversion to Fbc here
    command = "python convertCobraToFbc.py " + infile + " " + outfile
    os.system(command)

    # Read SBML File
    doc = reader.readSBML(outfile)
    model = doc.getModel()

    # Use BiGG Database for charge annotation, if none is given
    mismatches = pd.DataFrame(columns=["model_index", "id", "name", "charge_bigg", "charge_model", "formula_model"])
    num_spec = model.getNumSpecies()
    for i in tqdm(range(num_spec)):

        # check fbc plugin for formula
        if model.getSpecies(i).getPlugin('fbc').isSetCharge():
            continue

        meta_id = str(model.getSpecies(i).getMetaId())
        pruned_id = meta_id[2:-2]
        metabolite_info = hf.bigg_request(pruned_id, "metabolites")
        charges_bigg = metabolite_info["charges"]

        if len(charges_bigg) == 1:
            model.getSpecies(i).getPlugin('fbc').setCharge(charges_bigg[0])
            note_str = f"Changed charge from '' to {charges_bigg[0]}. Source: BiGG"
            model = hf.add_note_species(model, note_str, meta_id)
            lnk = f"https://identifiers.org/bigg.reaction:{pruned_id}"
            model = hf.add_link_annotation_species(model, lnk, meta_id)

        else:
            mismatches.loc[len(mismatches.index)] = [i, meta_id, model.getSpecies(i).getName(),
                                                     charges_bigg, "",
                                                     model.getSpecies(meta_id).getPlugin('fbc').getChemicalFormula()]

    # Saving new model
    doc.setModel(model)
    writer.writeSBML(doc, outfile)

    # Exporting mismatches
    mismatches.to_csv(outfile_tsv, sep="\t")


if __name__ == '__main__':
    main(sys.argv)
