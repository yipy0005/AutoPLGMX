from typing import Dict, Optional, Union, List

amino_acids_3_letter_codes = [
    "ARG",
    "HIS",
    "HIP",
    "HID",
    "HIE",
    "LYS",
    "ASP",
    "GLU",
    "SER",
    "THR",
    "ASN",
    "GLN",
    "CYS",
    "GLY",
    "PRO",
    "ALA",
    "VAL",
    "ILE",
    "LEU",
    "MET",
    "PHE",
    "TYR",
    "TRP",
    "TPO",
]


def get_pdb_line_components(pdb_line: str) -> Dict[str, Union[int, float, str]]:
    # Length of an ATOM line should have 78 characters
    pdb_line_terms: Dict[str, Union[int, float, str]] = {}

    if pdb_line:
        if (pdb_line.split())[0] in ["ATOM", "HETATM"]:
            pdb_line_terms["RECORD TYPE"] = pdb_line[:6].replace(" ", "")
            pdb_line_terms["ATOM SERIAL NUMBER"] = (
                int(pdb_line[6:11].replace(" ", ""))
                if pdb_line[6:11].replace(" ", "") != ""
                else pdb_line[6:11].replace(" ", "")
            )
            pdb_line_terms["ATOM NAME"] = pdb_line[12:16].replace(" ", "")
            pdb_line_terms["ALTERNATE LOCATION INDICATOR"] = pdb_line[16].replace(
                " ", ""
            )
            pdb_line_terms["RESIDUE NAME"] = pdb_line[17:20].replace(" ", "")
            pdb_line_terms["CHAIN ID"] = pdb_line[21].replace(" ", "")
            pdb_line_terms["RESIDUE SEQ NUMBER"] = (
                int(pdb_line[22:26].replace(" ", ""))
                if pdb_line[22:26].replace(" ", "") != ""
                else pdb_line[22:26].replace(" ", "")
            )
            pdb_line_terms["CODE FOR INSERTIONS OF RESIDUES"] = pdb_line[26].replace(
                " ", ""
            )
            pdb_line_terms["X COORDINATE"] = (
                float(pdb_line[30:38].replace(" ", ""))
                if pdb_line[30:38].replace(" ", "") != ""
                else pdb_line[30:38].replace(" ", "")
            )
            pdb_line_terms["Y COORDINATE"] = (
                float(pdb_line[38:46].replace(" ", ""))
                if pdb_line[38:46].replace(" ", "") != ""
                else pdb_line[38:46].replace(" ", "")
            )
            pdb_line_terms["Z COORDINATE"] = (
                float(pdb_line[46:54].replace(" ", ""))
                if pdb_line[46:54].replace(" ", "") != ""
                else pdb_line[46:54].replace(" ", "")
            )
            pdb_line_terms["OCCUPANCY"] = (
                float(pdb_line[54:60].replace(" ", ""))
                if pdb_line[54:60].replace(" ", "") != ""
                else pdb_line[54:60].replace(" ", "")
            )
            pdb_line_terms["TEMPERATURE FACTOR"] = (
                float(pdb_line[60:66].replace(" ", ""))
                if pdb_line[60:66].replace(" ", "") != ""
                else pdb_line[60:66].replace(" ", "")
            )
            pdb_line_terms["SEGMENT ID"] = pdb_line[72:76].replace(" ", "")
            pdb_line_terms["ELEMENT SYMBOL"] = pdb_line[76:78].replace(" ", "")
        elif (pdb_line.split())[0] == "TER":
            pdb_line_terms["RECORD TYPE"] = pdb_line[:6].replace(" ", "")
            try:
                pdb_line_terms["ATOM SERIAL NUMBER"] = int(
                    pdb_line[6:11].replace(" ", "")
                )
                pdb_line_terms["RESIDUE NAME"] = pdb_line[17:20].replace(" ", "")
                pdb_line_terms["CHAIN ID"] = pdb_line[21].replace(" ", "")
                pdb_line_terms["RESIDUE SEQ NO"] = (
                    int(pdb_line[22:26].replace(" ", ""))
                    if pdb_line[22:26].replace(" ", "") != ""
                    else pdb_line[22:26].replace(" ", "")
                )
                # pdb_line_terms["CODE FOR INSERTIONS OF RESIDUES"] = pdb_line[
                #     26
                # ].replace(" ", "")
            except ValueError:
                pass
        elif (pdb_line.split())[0] == "CONECT":
            pdb_line_terms["RECORD TYPE"] = pdb_line[:6].replace(" ", "")
            pdb_line_terms["ATOM SERIAL NUMBER"] = pdb_line[6:11].replace(" ", "")
            pdb_line_terms["BONDED ATOM SERIAL NUMBER"] = pdb_line[12:].strip("\n")

    return pdb_line_terms


def assemble_pdb_line_components(
    pdb_line_components: Dict[str, Union[int, float, str]]
) -> Optional[str]:
    # Length of an ATOM line should have 78 characters
    if pdb_line_components["RECORD TYPE"] in ["ATOM", "HETATM"]:
        for key in [
            "X COORDINATE",
            "Y COORDINATE",
            "Z COORDINATE",
            "OCCUPANCY",
            "TEMPERATURE FACTOR",
        ]:
            pdb_line_components[key] = f"{float(pdb_line_components[key]):.3f}"

        for key in [
            "OCCUPANCY",
            "TEMPERATURE FACTOR",
        ]:
            pdb_line_components[key] = f"{float(pdb_line_components[key]):.2f}"

        pdb_line = (
            f'{pdb_line_components["RECORD TYPE"]: <6}'
            f'{pdb_line_components["ATOM SERIAL NUMBER"]: >5} '
            f'{pdb_line_components["ATOM NAME"]: <4}'
            f'{pdb_line_components["ALTERNATE LOCATION INDICATOR"]: ^1}'
            f'{pdb_line_components["RESIDUE NAME"]: >3} '
            f'{pdb_line_components["CHAIN ID"]: ^1}'
            f'{pdb_line_components["RESIDUE SEQ NUMBER"]: >4}'
            f'{pdb_line_components["CODE FOR INSERTIONS OF RESIDUES"]: ^1}   '
            f'{pdb_line_components["X COORDINATE"]: >8}'
            f'{pdb_line_components["Y COORDINATE"]: >8}'
            f'{pdb_line_components["Z COORDINATE"]: >8}'
            f'{pdb_line_components["OCCUPANCY"]: >6}'
            f'{pdb_line_components["TEMPERATURE FACTOR"]: >6}      '
            f'{pdb_line_components["SEGMENT ID"]: <4}'
            f'{pdb_line_components["ELEMENT SYMBOL"]: >2}'
        )

        assert len(pdb_line) >= 78, f"{pdb_line} {len(pdb_line)}"
        return pdb_line


def get_missing_terms(pdb_line_terms: Dict[str, Union[int, float, str]]) -> List[str]:

    important_columns = [
        "RECORD TYPE",
        "ATOM SERIAL NUMBER",
        "ATOM NAME",
        "RESIDUE NAME",
        "RESIDUE SEQ NUMBER",
        "X COORDINATE",
        "Y COORDINATE",
        "Z COORDINATE",
        "OCCUPANCY",
        "TEMPERATURE FACTOR",
        "ELEMENT SYMBOL",
    ]

    return [
        key
        for key in pdb_line_terms
        if pdb_line_terms[key] == "" and key in important_columns
    ]


def check_missing_pdb_info(pdb_contents: List[str]) -> Dict[str, str]:
    problematic_lines: Dict[str, str] = {}
    for line in pdb_contents:
        pdb_line_terms = get_pdb_line_components(line)

        all_missing_terms = get_missing_terms(pdb_line_terms)
        if len(all_missing_terms) > 0:
            for missing_term in all_missing_terms:
                problematic_lines[line] = f"{missing_term}"

    return problematic_lines
