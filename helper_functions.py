"""
Useful functions for processing SMBL Data
"""
import requests
import sys
import re
import libsbml
import xmltodict
import json


def delete_doubles(arr):
    """
    :param arr: list()
    :return: Given list, without duplicated entries
    """
    arr2 = []
    for element in arr:
        if not arr2.__contains__(element):
            arr2.append(element)
    return arr2


def compare_formulas(formulas, charges=[]):
    """
    compares formulas
    :param charges: charges of molecules as list of integers
    :param charge_hydrogen_balance_accepted: boolean, that indecates, wether differences in the number of H atoms are
    accepted, if the difference is accounted for in the charge
    :param formulas: list of formulas
    :return: True, if all formulas have the same components with the same amount
    <formula> must be string: Upper case character marks new Element  (Mg12Ag2  = [Mg12, Ag2] & MG12AG2 = [M,G12,A,G2])
    """
    # Separate Components of Formula
    formulas_split = []
    for formula in formulas:
        formula_split = []

        # separates the atoms
        for char in formula:
            if char.isupper():
                formula_split.append(char)
            else:
                formula_split[len(formula_split) - 1] = formula_split[len(formula_split) - 1] + char

        # adds "1" to formula, if no number is given at the end
        for i in range(len(formula_split)):
            if re.search("[0-9]", formula_split[i]) is None:
                formula_split[i] = formula_split[i] + "1"

        # adds separated formulas to a list
        formulas_split.append(formula_split)

    # Iterates through all formulas
    for j in range(len(formulas_split) - 1):
        for component in formulas_split[j]:
            # accounts for hydrogen - charge relationship
            if charges and not re.search("^H(?![a-z])+([0-9])*", component) is None:
                component = int(component.split("H")[1])
                component = component + (charges[j + 1] - charges[j])
                component = "H" + str(component)

            # Check next element for current element
            if component not in formulas_split[j + 1]:
                return False

        # Check whether all components were in formula
        if len(formulas_split[j]) != len(formulas_split[j + 1]):
            return False

    return True


# ++ Get Metabolite Data from BiGG Database ++
def bigg_request(_id: str, search_type: str = "metabolites"):
    """
    Requests an entry from the BIGG Database
    :param _id: str e.g. "nh3"
    :param search_type: str e.g. "metabolites"
    :return: decoded .json into dictionary
    """
    custom_request = "http://bigg.ucsd.edu/api/v2/universal/" + search_type + "/" + _id
    req = requests.get(custom_request, headers={"Content-Type": "application/json"})

    if not req.ok:
        req.raise_for_status()
        sys.exit()

    decoded_req = req.json()
    return decoded_req


# ++ Get Metabolite Data from Biocyc Database ++
def biocyc_request(id_org: str, db: str, id_db: str):
    """
    Requests an entry from the BioCyc DB
    :param db: Database e.g. BIGG, SEED,..
    :param id_db: ID from Database e.g. atp, cpd0001
    :param id_org: ID of organism e.g. GCF_000010185
    :return: decoded .json into dictionary
    """
    custom_request = f"https://websvc.biocyc.org/{id_org}/foreignid?ids={db}:{id_db}&fmt=json"
    req = requests.get(custom_request, headers={"Content-Type": "application/json"})

    if not req.ok:
        req.raise_for_status()
        sys.exit()

    try:
        decoded_req = req.json()
    except json.decoder.JSONDecodeError:
        assert id_org != "meta"
        decoded_req = biocyc_request("meta", db, id_db)
    return decoded_req


# KEGG Request Function
def kegg_get(org_id: str, kegg_id: str):
    request_url = f"http://rest.kegg.jp/get/{org_id}:{kegg_id}"
    req = requests.get(request_url).text.split("\n")
    return req


# ++ Get Metabolite Data from BioCyc Database ++
def biocyc_get(id_org: str, id_biocyc: str, detail: str = "full"):
    """
    Requests an entry from the BioCyc DB
    :param detail: either none, low or full, defaults to full
    :param id_biocyc: ID of object e.g. ATP
    :param id_org: ID of organism e.g. GCF_000010185
    :return: decoded .xml into dictionary
    """
    custom_request = f"https://websvc.biocyc.org/getxml?id={id_org}:{id_biocyc}&detail={detail}"
    req = requests.get(custom_request)

    if not req.ok:
        req.raise_for_status()
        sys.exit()

    decoded_req = xmltodict.parse(req.content)
    return decoded_req


def biocyc_get_from_formula(id_org: str, formula: str):
    """
    Requests an entry from the BioCyc DB
    :param formula:
    :param detail: either none, low or full, defaults to full
    :param id_biocyc: ID of object e.g. ATP
    :param id_org: ID of organism e.g. GCF_000010185
    :return: decoded .xml into dictionary
    """
    custom_request = f"https://websvc.biocyc.org/{id_org}/CF?cfs={formula}&fmt=json"
    req = requests.get(custom_request)

    if not req.ok:
        req.raise_for_status()
        sys.exit()

    try:
        decoded_req = req.json()
    except json.decoder.JSONDecodeError:
        assert id_org != "meta"
        decoded_req = biocyc_get_from_formula("meta", formula)
    return decoded_req


def make_cv_term(link: str, qual_type=libsbml.BQB_IS):
    """
    :param qual_type:
    :param link: string that is added to CV-Term
    :return: libsbml.CVTerm
    This method is not generic, but only creates species and reaction standard CV Terms.
    """
    c = libsbml.CVTerm()
    c.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER)
    c.setBiologicalQualifierType(qual_type)
    c.addResource(link)
    return c


def add_link_annotation_species(model, lnk, qual_type, s_id):
    """
    :param qual_type: libsbml.QUALIFIER
    :param model: libsbml.model
    :param lnk: string
    :param s_id: string
    :return: libsbml.model
    """
    cv_term = make_cv_term(lnk, qual_type)

    # eliminate duplicates
    list_cv = []
    for i in range(model.getSpecies(s_id).getNumCVTerms()):
        list_cv.append(model.getSpecies(s_id).getCVTerm(i))

    if cv_term not in list_cv:
        model.getSpecies(s_id).addCVTerm(cv_term)
    return model


def add_link_annotation_reaction(model, lnk, qual_type, s_id):
    """
    :param qual_type: libsbml.QUALIFIER
    :param model: libsbml.model
    :param lnk: string
    :param s_id: string
    :return: libsbml.model
    """
    cv_term = make_cv_term(lnk, qual_type)

    # eliminate duplicates
    list_cv = []
    for k in range(model.getReaction(s_id).getNumCVTerms()):
        list_cv.append(model.getReaction(s_id).getCVTerm(k))

    if cv_term not in list_cv:
        model.getReaction(s_id).addCVTerm(cv_term)
    return model


def add_note_species(model, note: str, fbc_id):
    """
    :param fbc_id: str
    :param model: libsbml.model
    :param note: str
    :return: libsbml.model
    """
    str_note = f"<body  xmlns=\"http://www.w3.org/1999/xhtml\">\n  <p>{note}</p>\n  </body>"
    if not model.getSpecies(fbc_id).isSetNotes():
        model.getSpecies(fbc_id).setNotes(str_note)
    else:
        notes_curent = model.getSpecies(fbc_id).getNotes().toXMLString()
        if note not in notes_curent:
            model.getSpecies(fbc_id).appendNotes(str_note)
    return model


def add_note_gene_product(model, note: str, fbc_id):
    """
    :param fbc_id: str
    :param model: libsbml.model
    :param note: str
    :return: libsbml.model
    """
    str_note = f"<body  xmlns=\"http://www.w3.org/1999/xhtml\">\n  <p>{note}</p>\n  </body>"
    if not model.getPlugin('fbc').getGeneProduct(fbc_id).isSetNotes():
        model.getPlugin('fbc').getGeneProduct(fbc_id).setNotes(str_note)
    else:
        notes_curent = model.getPlugin('fbc').getGeneProduct(fbc_id).getNotes().toXMLString()
        if note not in notes_curent:
            model.getPlugin('fbc').getGeneProduct(fbc_id).appendNotes(str_note)
    return model


def add_note_reaction(model, note: str, fbc_id):
    """
    :param model: libsbml.model
    :param note: str
    :param fbc_id: str
    :return:
    """
    str_note = f"<body  xmlns=\"http://www.w3.org/1999/xhtml\">\n  <p>{note}</p>\n  </body>"
    if not model.getReaction(fbc_id).isSetNotes():
        model.getReaction(fbc_id).setNotes(str_note)
    else:
        notes_curent = model.getReaction(fbc_id).getNotes().toXMLString()
        if note not in notes_curent:
            model.getReaction(fbc_id).appendNotes(str_note)
    return model


def dict_add_overlap_to_list(orig_dict, extend_dict):
    for k, v in extend_dict.items():
        if k not in orig_dict:
            orig_dict[k] = v
        else:
            if hasattr(orig_dict[k], '__iter__') and not isinstance(orig_dict[k], str):
                orig_dict[k] = set(orig_dict[k])
            else:
                orig_dict[k] = {orig_dict[k]}

            if hasattr(v, '__iter__') and not isinstance(v, str):
                orig_dict[k] |= set(v)
            else:
                orig_dict[k] |= {v}
            orig_dict[k] = list(orig_dict[k])
    return orig_dict
