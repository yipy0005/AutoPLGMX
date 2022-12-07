import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List

from gql import Client, gql
from gql.transport.requests import RequestsHTTPTransport

from check_pdb import get_pdb_line_components

amino_acids_3_letter_codes = [
    "ARG",
    "HIS",
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
]


transport = RequestsHTTPTransport(
    url="https://data.rcsb.org/graphql",
)

client = Client(transport=transport)

protein_pdb_query = gql(
    """
    query get_pdb_info($pdb_id: String!){
        entry(entry_id: $pdb_id) {
            rcsb_binding_affinity {
            comp_id
            }
        }
    }
    """
)


@dataclass
class Commandline_Operation:
    params: List[str]
    input_files: Dict[str, Path]
    output_files: Dict[str, Path]
    input_for_subprocess: str = ""
    # additional_gmx_args: str = ""

    def run_gmx_cmdline(self) -> subprocess.CompletedProcess[str]:
        self.params.insert(0, "gmx")

        for file_group in [self.input_files, self.output_files]:
            for key, value in file_group.items():
                self.params.append(key)
                self.params.append(str(value))

        return subprocess.run(
            self.params, input=self.input_for_subprocess, encoding="utf-8"
        )


def download_pdb(pdb_id: str = "", pdb_file: Path = Path("")) -> List[str]:
    """This function downloads the PDB file from RCSB using its API and GraphQL.

    Args:
        pdb_code (str): the 4-letter code of the PDB file to be downloaded

    Returns:
        list[str]: the contents of the downloaded PDB file as a list of strings
    """

    record_types: list[str] = [
        "ATOM",
        "HETATM",
        "TER",
        "CONECT",
        "MASTER",
        "END",
    ]

    if pdb_id != "" and pdb_file == Path(""):
        pdb_contents: list[str] = (
            subprocess.run(
                f"pdb_fetch {pdb_id}",
                shell=True,
                capture_output=True,
                encoding="UTF-8,",
            )
        ).stdout.split("\n")

    elif not pdb_id and pdb_file != Path(""):
        with open(pdb_file, "r") as pdb_f:
            pdb_contents: list[str] = pdb_f.readlines()

    else:
        sys.exit("Please specify ONLY the PDB ID or the PDB FILE PATH!")

    return [
        line
        for line in pdb_contents
        if len(line) != 0 and (line.split())[0] in record_types
    ]


def checkConsecutive(lst: List[int]) -> List[int]:
    """
    This function checks for consecutive numbers in a list.

    :param lst: a list of numbers
    :return: a list of missing numbers
    """

    missing_numbers = [
        number
        for number in sorted(lst) + list(range(min(lst), max(lst) + 1))
        if number not in sorted(lst)
        or number not in list(range(min(lst), max(lst) + 1))
    ]

    key_missing_numbers: list[int] = []
    for i in range(len(missing_numbers)):
        try:
            if missing_numbers[i + 1] - missing_numbers[i] > 1:
                key_missing_numbers.append(missing_numbers[i + 1])
        except IndexError:
            pass

    try:
        key_missing_numbers.insert(0, missing_numbers[0])
        key_missing_numbers.insert(0, int("0"))
        key_missing_numbers.append(int("9999"))
    except IndexError:
        key_missing_numbers.extend((int("0"), int("9999")))
    return key_missing_numbers


def get_lines(lines_lst: List[str], first_line: str = "") -> List[str]:
    output_lst: list[str] = []
    for line_idx in range(len(lines_lst)):
        if first_line != "" and first_line in lines_lst[line_idx]:
            output_lst.append(lines_lst[line_idx])
            try:
                while len(lines_lst[line_idx]) != 0:
                    line_idx += 1
                    if len(lines_lst[line_idx]) != 0:
                        output_lst.append(lines_lst[line_idx])
                    else:
                        break
            except IndexError:
                pass

    return output_lst


def read_topol(pdb_type: str, topol_file: Path, ligid: str) -> Dict[str, List[str]]:
    topol_contents: Dict[str, List[str]] = {}

    with open(topol_file, "r") as topol_f:
        lines = [line.strip("\n") for line in topol_f.readlines()]

        topol_contents["MOLECULE TYPE"] = get_lines(lines, "[ moleculetype ]")
        topol_contents["ATOMS"] = get_lines(lines, "[ atoms ]")
        topol_contents["BONDS"] = get_lines(lines, "[ bonds ]")
        topol_contents["PAIRS"] = get_lines(lines, "[ pairs ]")
        topol_contents["ANGLES"] = get_lines(lines, "[ angles ]")
        topol_contents["DIHEDRALS1"] = get_lines(
            lines,
            (
                ";  ai    aj    ak    al funct            c0            c1            c2"
                "            c3            c4            c5"
            ),
        )
        if len(topol_contents["DIHEDRALS1"]) != 0:
            topol_contents["DIHEDRALS1"].insert(0, "[ dihedrals ]")
        topol_contents["DIHEDRALS2"] = get_lines(
            lines,
            ";  ai    aj    ak    al funct            c0            c1            c2            c3",
        )
        if len(topol_contents["DIHEDRALS2"]) != 0:
            topol_contents["DIHEDRALS2"].insert(0, "[ dihedrals ]")
        topol_contents["CMAP"] = get_lines(lines, "[ cmap ]")
        topol_contents["POSRE"] = get_lines(lines, "; Include Position restraint file")

        if Path.exists(Path("posre.itp")) is False:
            topol_contents["CHAIN TOPOLOGIES"] = get_lines(
                lines, "; Include chain topologies"
            )
        else:
            topol_contents["CHAIN TOPOLOGIES"] = []
        topol_contents["SYSTEM"] = get_lines(lines, "[ system ]")
        topol_contents["MOLECULES"] = get_lines(lines, "[ molecules ]")

        if pdb_type in {"COMPLEX", "PROTEIN"}:
            topol_contents["FORCEFIELD PARAMS"] = get_lines(
                lines, "; Include forcefield parameters"
            )
            topol_contents["WATER TOPOLOGY"] = get_lines(
                lines, "; Include water topology"
            )
            topol_contents["WATER POSRE"] = get_lines(lines, "#ifdef POSRES_WATER")
            topol_contents["IONS TOPOLOGY"] = get_lines(
                lines, "; Include topology for ions"
            )

        else:
            topol_contents["FORCEFIELD PARAMS"] = [
                "; Include forcefield parameters",
                '#include "amber99sb.ff/forcefield.itp"',
            ]
            topol_contents["WATER TOPOLOGY"] = [
                "; Include water topology",
                '#include "amber99sb.ff/tip3p.itp"',
            ]
            topol_contents["WATER POSRE"] = [
                "#ifdef POSRES_WATER",
                "; Position restraint for each water oxygen",
                "[ position_restraints ]",
                ";  i funct       fcx        fcy        fcz",
                "1    1       1000       1000       1000",
                "#endif",
            ]
            topol_contents["IONS TOPOLOGY"] = [
                "; Include topology for ions",
                '#include "amber99sb.ff/ions.itp"',
            ]

        if pdb_type in {"COMPLEX", "LIGAND"}:
            topol_contents["LIGAND TOPOLOGY"] = [
                "; Include ligand topology",
                f'#include "{ligid}_fix.acpype/{ligid}_fix_GMX.itp"',
            ]
            topol_contents["LIGAND POSRE"] = [
                "; Ligand position restraints",
                "#ifdef POSRES",
                f'#include "{ligid}_fix.acpype/posre_{ligid}_fix.itp"',
                "#endif",
            ]

            flat_lst: list[Any] = []
            for line in topol_contents["MOLECULES"]:
                terms = line.split()
                flat_lst.extend(iter(terms))
            if f"{ligid}_fix" not in flat_lst:
                topol_contents["MOLECULES"].append(f"{ligid}_fix        1")

    return topol_contents


def write_prm_file(lig_itp_file: Path) -> List[str]:
    """Extracts ligand topology information from the ligand.itp file

    Args:
        lig_itp_file (Path): Path to ligand.itp file

    Returns:
        list[str]: returns contents of the ligand topology information as a string
    """

    with open(lig_itp_file, "r") as lig_itp_f:
        lines = lig_itp_f.readlines()
        return lines[:1] + lines[10:]


def check_pdb_type(pdb_contents: List[str]) -> str:

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

    resname_lst: list[str] = []
    for line in pdb_contents:
        pdb_terms_dict = get_pdb_line_components(line)
        try:
            resname_lst.append(str(pdb_terms_dict["RESIDUE NAME"]))
        except KeyError:
            pass
    resname_lst = list(set(resname_lst))

    check = all(item in amino_acids_3_letter_codes for item in resname_lst)

    if check:
        return "PROTEIN"
    count = sum(item in amino_acids_3_letter_codes for item in resname_lst)
    return "LIGAND" if count == 0 else "COMPLEX"


def calc_net_charge(pdb_contents: List[str]) -> int:
    net_charge: int = 0
    for line in pdb_contents:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_chrg = ((line.split()[-1])[1:])[::-1]
            try:
                if atom_chrg[0] == "+":
                    net_charge += int(atom_chrg[1])
                elif atom_chrg[0] == "-":
                    net_charge -= int(atom_chrg[1])
            except IndexError:
                pass

    return net_charge
