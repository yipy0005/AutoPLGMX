import string
import subprocess
from typing import Dict, Optional, Tuple, Union

from gql import Client, gql
from gql.transport.requests import RequestsHTTPTransport
from pipe import Pipe  # type: ignore

from check_pdb import assemble_pdb_line_components, get_pdb_line_components
from utils import checkConsecutive

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


def check_missing_residues(pdb_contents: list[str]) -> bool:
    """
    Checks if there are any missing residues in the pdb file.

    :param pdb_contents: A list of the lines in the pdb file.
    :return: True if there are missing residues, False otherwise.
    """

    residue_no_lst: list[int] = []
    for line in pdb_contents[:-1]:
        terms = line.split()
        if "ATOM" in terms[0]:
            try:
                residue_no_lst.append(int(terms[5]))
            except ValueError:
                pass

    return checkConsecutive(list(set(residue_no_lst))) is None


@Pipe
def remove_alternate_positions(pdb_contents: list[str]) -> list[str]:
    """
    Remove alternate positions from a PDB file.

    Parameters
    ----------
    pdb_contents : list
        A list of strings representing the contents of a PDB file.

    Returns
    -------
    list
        A list of strings representing the contents of a PDB file with alternate positions removed.
    """

    return (
        subprocess.run(
            ["pdb_selaltloc"],
            input="\n".join(pdb_contents),
            capture_output=True,
            text=True,
        )
    ).stdout.split("\n")


@Pipe
def delete_hydrogens(pdb_contents: list[str]) -> list[str]:
    """
    Takes a list of strings representing a PDB file and returns the same list with all hydrogen atoms
    removed.
    :param pdb_contents: list of strings representing a PDB file
    :return: list of strings representing a PDB file with all hydrogen atoms removed
    """

    edited_pdb_contents: list[str] = []
    for line in pdb_contents:
        if (
            line.split()[0] in ["ATOM", "HETATM"]
            and line.split()[-1] != "H"
            or line.split()[0] == "CONECT"
        ):
            edited_pdb_contents.append(line)

    return edited_pdb_contents


@Pipe
def delete_hetatms(pdb_contents: list[str]) -> list[str]:
    """
    Delete all HETATM and CONECT lines from a PDB file.

    Args:
        pdb_contents (list[str]): The contents of a PDB file.

    Returns:
        list[str]: The contents of a PDB file with all HETATM and CONECT lines removed.
    """

    edited_pdb_contents: list[str] = []
    for line in pdb_contents:
        if (get_pdb_line_components(line))["RECORD TYPE"] == "CONECT":
            edited_pdb_contents.append(line)
        elif (get_pdb_line_components(line))[
            "RESIDUE NAME"
        ] in amino_acids_3_letter_codes:
            edited_pdb_contents.append(line)

    return edited_pdb_contents


@Pipe
def assign_chain_ids(pdb_contents: list[str]) -> list[str]:
    edited_pdb_contents: Optional[list[str]] = []

    all_residue_numbers: list[int] = []

    for line in pdb_contents:
        if line.startswith("CONECT") is False:
            pdb_line_components = get_pdb_line_components(line)
            all_residue_numbers.append(int(pdb_line_components["RESIDUE SEQ NUMBER"]))

    all_residue_numbers = list(set(all_residue_numbers))
    missing_residue_numbers: list[int] = checkConsecutive(all_residue_numbers)

    for i in range(len(missing_residue_numbers)):
        try:
            for line in pdb_contents:
                pdb_line_components: Optional[
                    Dict[str, Union[int, float, str]]
                ] = get_pdb_line_components(line)
                if line.startswith("CONECT") is False:
                    if (
                        int(pdb_line_components["RESIDUE SEQ NUMBER"])
                        > missing_residue_numbers[i]
                        and int(pdb_line_components["RESIDUE SEQ NUMBER"])
                        < missing_residue_numbers[i + 1]
                    ):
                        pdb_line_components["CHAIN ID"] = string.ascii_uppercase[i]
                        edited_pdb_contents.append(
                            str(assemble_pdb_line_components(pdb_line_components))
                        )
        except IndexError:
            pass

    return edited_pdb_contents


@Pipe
def adding_ter(pdb_contents: list[str]) -> list[str]:
    """
    Add TER to the end of each chain in the PDB file.

    :param pdb_contents: list of strings
    :return: list of strings
    """
    edited_pdb_contents: list[str] = []

    for i in range(len(pdb_contents)):
        current_pdb_line_components = get_pdb_line_components(pdb_contents[i])
        edited_pdb_contents.append(pdb_contents[i])

        if i == len(pdb_contents) - 1:
            break

        next_pdb_line_components = get_pdb_line_components(pdb_contents[i + 1])
        if (
            next_pdb_line_components["RECORD TYPE"] != "TER"
            and int(next_pdb_line_components["RESIDUE SEQ NUMBER"])
            - int(current_pdb_line_components["RESIDUE SEQ NUMBER"])
            > 1
        ):
            edited_pdb_contents.append("TER")

    if pdb_contents[-1] != "TER" or (pdb_contents[-1].split())[0] != "TER":
        edited_pdb_contents.append("TER")

    return edited_pdb_contents


def get_bound_ligand(pdb_id: str, pdb_contents: list[str]) -> Tuple[str, list[str]]:
    """This function extracts only the lines in the PDB file pertaining to the ligand molecule.

    Args:
        pdb_id (str): the 4-letter code of the PDB file

    Returns:
        Tuple[str, list[str]]: returns the residue name/id of the ligand molecule and a list containing
        the lines of the ligand molecule in the PDB file
    """

    if pdb_id != "":
        response = client.execute(  # type: ignore
            protein_pdb_query, variable_values={"pdb_id": pdb_id}
        )

        try:
            bound_ligand_pdbid: str = response["entry"]["rcsb_binding_affinity"][0][
                "comp_id"
            ]

            bound_ligand_pdb: list[str] = []

            for line in pdb_contents:
                if "HETATM" in line and bound_ligand_pdbid in line:
                    bound_ligand_pdb.append(line)
                else:
                    pass

            return bound_ligand_pdbid, bound_ligand_pdb
        except TypeError:
            print("No bound ligands found from RCSB!")

            return "", []
    else:
        hetero_resname: list[str] = []
        for line in pdb_contents:
            pdb_terms_dict = get_pdb_line_components(line)
            if pdb_terms_dict["RECORD TYPE"] == "HETATM":
                hetero_resname.append(str(pdb_terms_dict["RESIDUE NAME"]))
            elif pdb_terms_dict["RECORD TYPE"] == "ATOM":
                if pdb_terms_dict["RESIDUE NAME"] not in amino_acids_3_letter_codes:
                    hetero_resname.append(str(pdb_terms_dict["RESIDUE NAME"]))
        hetero_resname = list(set(hetero_resname))

        choices: Dict[int, str] = {}
        if len(hetero_resname) > 1:
            i = 1
            for resname in hetero_resname:
                choices[i] = resname
                i += 1

            for key, value in choices.items():
                print(f"{key}. {value}")

            bound_ligand_pdbid = choices[int(input("Choose Bound Ligand: "))]

            bound_ligand_pdb: list[str] = []

            for line in pdb_contents:
                if "HETATM" in line or "ATOM" in line and bound_ligand_pdbid in line:
                    bound_ligand_pdb.append(line)
                else:
                    pass

            dehydrogenated_ligand: list[str] = bound_ligand_pdb | delete_hydrogens

            return bound_ligand_pdbid, dehydrogenated_ligand

        elif len(hetero_resname) == 1:
            bound_ligand_pdbid = hetero_resname[0]

            bound_ligand_pdb: list[str] = []

            for line in pdb_contents:
                if "HETATM" in line or "ATOM" in line and bound_ligand_pdbid in line:
                    bound_ligand_pdb.append(line)
                else:
                    pass

            dehydrogenated_ligand: list[str] = bound_ligand_pdb | delete_hydrogens

            return bound_ligand_pdbid, dehydrogenated_ligand
        else:
            print("No Ligands Found!")
            return "", []


@Pipe
def delete_waters(pdb_contents: list[str]) -> list[str]:
    """
    Removes water molecules from a PDB file.
    """

    return [line for line in pdb_contents if "HOH" not in line]
