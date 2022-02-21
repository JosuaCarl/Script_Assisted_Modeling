import sys
import os
import libsbml
from tqdm import tqdm

'''
Usage: annotate_links_from_mp.py <path_model-polisher_sbml-file> <path_input_sbml-file> <path_output_sbml-file>
Extracts annotations form ModelPolisher file and adds it to previous file.
'''


def main(args):
    # console access
    if len(args) != 6:
        print(f"Arguments: {len(args)}")
        print(main.__doc__)
        sys.exit(1)

    infile_mp = args[1]
    infile_of = args[2]
    outfile = args[3]

    if not os.path.exists(infile_mp):
        print("[Error] %s : No such file." % infile_mp)
        sys.exit(1)

    # create Readers and Writers
    reader = libsbml.SBMLReader()
    writer = libsbml.SBMLWriter()

    # Read SBML File
    doc_mp = reader.readSBML(infile_mp)
    model_mp = doc_mp.getModel()

    doc_nf = reader.readSBML(infile_of)
    model_nf = doc_nf.getModel()

    # Compare Annotations
    elements_mp = list(model_mp.getListOfAllElements())
    for k in tqdm(range(len(elements_mp))):
        e_id = elements_mp[k].getMetaId()

        # checks whether the element-id is valid and CV Terms exist for this id
        if model_mp.getElementByMetaId(e_id) is None or model_nf.getElementByMetaId(e_id) is None or\
                model_mp.getElementByMetaId(e_id).getCVTerm(0) is None:
            continue

        # Make list of all annotated links in ModelPolisher Model
        list_mp = []
        for i in range(model_mp.getElementByMetaId(e_id).getCVTerm(0).getNumResources()):
            list_mp.append(model_mp.getElementByMetaId(e_id).getCVTerm(0).getResources().getValue(i))

        # Makes List of all previously annotated links
        list_nf = []
        if not model_nf.getElementByMetaId(e_id).getCVTerm(0) is None:
            for j in range(model_nf.getElementByMetaId(e_id).getCVTerm(0).getNumResources()):
                list_nf.append(model_nf.getElementByMetaId(e_id).getCVTerm(0).getResources().getValue(j))

        for lnk in list_mp:
            if not list_nf.__contains__(lnk) and lnk.startswith("http"):
                element = model_nf.getElementByMetaId(e_id)

                if element.getAnnotation() is None:
                    c = libsbml.CVTerm()
                    c.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER)

                    # setting right biological qualifier
                    if e_id.startswith("C_"):
                        c.setBiologicalQualifierType(libsbml.BQB_IS_DESCRIBED_BY)
                    else:
                        c.setBiologicalQualifierType(libsbml.BQB_IS)
                    c.addResource(lnk)

                    model_nf.getElementByMetaId(e_id).addCVTerm(c)

                elif element.getCVTerm(0) is None:
                    triple = libsbml.XMLTriple("li", "rdf", "rdf")
                    attributes = libsbml.XMLAttributes()
                    attributes.add("resource", lnk, "li", "rdf")
                    li_node = libsbml.XMLNode(triple, attributes)

                    model_nf.getElementByMetaId(e_id).getAnnotation().getChild(0).getChild(0).getChild(0).getChild(0).\
                        addChild(li_node)

                else:
                    model_nf.getElementByMetaId(e_id).getCVTerm(0).addResource(lnk)

    # Saving new model
    doc_nf.setModel(model_nf)
    writer.writeSBML(doc_nf, outfile)


if __name__ == '__main__':
    main(sys.argv)
