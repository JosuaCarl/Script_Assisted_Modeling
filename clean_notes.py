import sys
import os
import libsbml
from tqdm import tqdm
import pandas as pd
import helper_functions as hf

'''
Usage: clean_notes.py <path_input-file> <path_output-file> <infile-csv> <program_name> <program_version> 
'''


def main(args):
    # console access
    if len(args) < 5:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    program_name = args[3]
    program_version = args[4]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # create Readers and Writers
    reader = libsbml.SBMLReader()
    writer = libsbml.SBMLWriter()

    # Read SBML File
    doc = reader.readSBML(infile)
    model = doc.getModel()

    num_spec = model.getNumSpecies()
    for i in tqdm(range(num_spec)):
        if not model.getSpecies(i).isSetNotes():
            continue
        for j in reversed(range(model.getSpecies(i).getNotes().getChild(0).getNumChildren())):
            note_str = model.getSpecies(i).getNotes().getChild(0).getChild(j).getChild(0).toXMLString().split(" ")

            if note_str[0] == "FORMULA:":
                if not model.getSpecies(i).getPlugin('fbc').isSetChemicalFormula():
                    model.getSpecies(i).getPlugin('fbc').setChemicalFormula(note_str[1])

            elif note_str[0] == "CHARGE:":
                if not model.getSpecies(i).getPlugin('fbc').isSetCharge():
                    model.getSpecies(i).getPlugin('fbc').setCharge(int(note_str[1]))

            elif note_str[0] == "SBOTerm:":
                if not model.getSpecies(i).isSetSBOTerm():
                    model.getSpecies(i).setSBOTerm(note_str[1])

            else:
                continue

            model.getSpecies(i).getNotes().getChild(0).removeChild(j)

    num_reac = model.getNumReactions()
    for i in tqdm(range(num_reac)):
        if not model.getReaction(i).isSetNotes():
            continue
        for j in reversed(range(model.getReaction(i).getNotes().getChild(0).getNumChildren())):
            note_str = model.getReaction(i).getNotes().getChild(0).getChild(j).getChild(0).toXMLString().split(" ")

            if note_str[0] == "SBOTerm:":
                if not model.getReaction(i).isSetSBOTerm():
                    model.getReaction(i).setSBOTerm(note_str[1])

            else:
                continue

            model.getReaction(i).getNotes().getChild(0).removeChild(j)

    num_comp = model.getNumCompartments()
    for i in tqdm(range(num_comp)):
        if not model.getCompartment(i).isSetNotes():
            continue
        for j in reversed(range(model.getCompartment(i).getNotes().getChild(0).getNumChildren())):
            note_str = model.getCompartment(i).getNotes().getChild(0).getChild(j).getChild(0).toXMLString().split(" ")

            if note_str[0] == "SBOTerm:":
                if not model.getCompartment(i).isSetSBOTerm():
                    model.getCompartment(i).setSBOTerm(note_str[1])

            else:
                continue

            model.getCompartment(i).getNotes().getChild(0).removeChild(j)

    num_gene_prod = model.getPlugin('fbc').getNumGeneProducts()
    for i in tqdm(range(num_gene_prod)):
        if not model.getPlugin('fbc').getGeneProduct(i).isSetNotes():
            continue
        for j in reversed(range(model.getPlugin('fbc').getGeneProduct(i).getNotes().getChild(0).getNumChildren())):
            note_str = model.getPlugin('fbc').getGeneProduct(i).getNotes().getChild(0).getChild(j).getChild(0).toXMLString().split(" ")

            if note_str[0] == "SBOTerm:":
                if not model.getPlugin('fbc').getGeneProduct(0).isSetSBOTerm():
                    model.getPlugin('fbc').getGeneProduct(i).setSBOTerm(note_str[1])

            else:
                continue

            model.getPlugin('fbc').getGeneProduct(i).getNotes().getChild(0).removeChild(j)

    # Saving new model
    doc.setModel(model)
    writer.setProgramName(program_name)
    writer.setProgramVersion(program_version)
    writer.writeSBML(doc, outfile)


if __name__ == '__main__':
    main(sys.argv)
