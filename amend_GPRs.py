import sys
import os
import cobra
from tqdm import tqdm
import gffpandas.gffpandas as gffpd
import re
import memote
from bioservices.kegg import KEGG

'''
Usage: amend_GPRs.py <path_input_sbml-file> <path_output_sbml-file> <path GFF file> <path_output-memote>
Adds Gene Protein Reaction rules to model. Based on heavy heuristics (naming).
'''


def main(args):
    # console access
    if len(args) != 5:
        print(main.__doc__)
        sys.exit(1)

    infile = args[1]
    outfile = args[2]
    gff_file = args[3]
    memote_report = args[4]

    if not os.path.exists(infile):
        print("[Error] %s : No such file." % infile)
        sys.exit(1)

    # read model
    model = cobra.io.read_sbml_model(infile)
    kegg = KEGG()

    # read GFF file
    df = gffpd.read_gff3(gff_file)
    df_attr = df.attributes_to_columns()
    org_id = "FMA"

    for i in tqdm(range(len(model.reactions))):

        # Extract all genes(locus-tag) into list
        if "ec-code" in model.reactions[i].annotation:
            ec_codes = model.reactions[i].annotation["ec-code"]
            if not isinstance(ec_codes, list):
                ec_codes = [ec_codes]
            locus_tags = []
            for ec in ec_codes:
                kegg_enzyme = kegg.parse(kegg.get(ec))
                if "GENES" in kegg_enzyme and org_id in kegg_enzyme["GENES"]:
                    locus_tags.append(kegg_enzyme["GENES"][org_id])
            if not locus_tags:
                continue

            # extract all subunits into dict with key being the enzyme
            enzyme_subunit = dict()
            for locus_tag in locus_tags:
                kegg_gene = kegg.parse(kegg.get(org_id.lower() + ":" + locus_tag))

                # Get sbml_name of gene from GFF File
                locus_match = df_attr.loc[df_attr["old_locus_tag"] == locus_tag]
                if locus_match.empty:
                    continue
                index = locus_match.index[0]
                name = df_attr.loc[index + 1, "Name"]
                id_sbml = name[:-2] + "_" + name[-1]

                if "NAME" in kegg_gene:
                    if not isinstance(kegg_gene["NAME"], list):
                        kegg_gene["NAME"] = [kegg_gene["NAME"]]
                    for name in kegg_gene["NAME"]:
                        if "subunit" in name:
                            enzyme_name = re.sub(r"( [A-Z])? subunit", "", name)
                            if enzyme_name in enzyme_subunit:
                                enzyme_subunit[enzyme_name] += f" and {id_sbml}"
                            else:
                                enzyme_subunit[enzyme_name] = id_sbml
                        else:
                            enzyme_subunit[name] = id_sbml

            # make GPR String
            gpr_str = ""
            for subunit_str in enzyme_subunit.values():
                if "and" in subunit_str:
                    gpr_str += f" or ( {subunit_str} )"
                else:
                    gpr_str += f" or {subunit_str}"
            gpr_str = gpr_str[4:]    # crops off first or incl. whitespaces

            # Sanity check
            if gpr_str == "":
                continue

            # assigns GPR
            if model.reactions[i].gene_reaction_rule == "" or model.reactions[i].gene_reaction_rule is None:
                model.reactions[i].gene_reaction_rule = gpr_str
            else:
                model.reactions[i].gene_reaction_rule += f" or {gpr_str}"

    # Export model
    cobra.io.write_sbml_model(model, outfile)

    # Make memote report
    result = memote.test_model(model, results=True)  #, skip=["test_find_metabolites_not_produced_with_open_bounds"])
    report = memote.snapshot_report(result[1], config=None, html=True)
    with open(memote_report, "w") as handle:
        handle.write(report)


if __name__ == '__main__':
    main(sys.argv)
