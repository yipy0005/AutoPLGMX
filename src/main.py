import argparse
import shutil
import subprocess
import sys

# from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Tuple

from check_pdb import check_missing_pdb_info, get_pdb_line_components
from gmx_ops import (  # convert_lig_to_gro,
    MDP,
    combine_prot_lig_gro,
    edit_mdp_files,
    edit_topol_top,
    fix_lig_bond_order,
    generate_lig_files,
    generate_prot_top,
    hydrogenate_lig,
)
from md import (
    AddIons,
    Minimization,
    NPT_Equilibration,
    NVT_Equilibration,
    Production_MD,
    Solvation,
)
from postprocessing import remove_pbc, rot_trans_fit
from protprep import (
    adding_ter,
    assign_chain_ids,
    delete_hetatms,
    delete_hydrogens,
    delete_waters,
    remove_alternate_positions,
)
from utils import Commandline_Operation, check_pdb_type, download_pdb

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


def obtain_pdb(pdbid: str, pdb_file: Path) -> Tuple[list[str], str]:
    pdb_contents = download_pdb(pdbid, pdb_file)
    pdb_filename = pdbid if pdb_file == Path("") else Path((pdb_file).name).stem
    return pdb_contents, pdb_filename


def clean_pdb(pdb_type: str, pdb_contents: list[str]) -> list[str]:

    return (
        (
            pdb_contents
            | remove_alternate_positions
            | delete_hydrogens
            | delete_hetatms
            | assign_chain_ids
            | adding_ter
        )
        if pdb_type in {"PROTEIN", "COMPLEX"}
        else pdb_contents | delete_hydrogens
    )  # type: ignore


def mdp_type(pdb_type: str, simulation_type: str, ligid: str) -> Any:
    if simulation_type == "ions":
        mdp = MDP(
            sim_type=simulation_type,
            integrator="steep",
            nsteps=50000,
            nstlist=1,
            rlist=1.0,
            coulombtype="cutoff",
            rcoulomb=1.0,
            rvdw=1.0,
        )

        return mdp

    elif simulation_type == "em":
        mdp = MDP(
            sim_type=simulation_type,
            integrator="steep",
            nsteps=50000,
            nstlist=1,
        )

        return mdp

    elif simulation_type == "nvt" and pdb_type == "COMPLEX":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=500000,
            # nsteps=50000,
            nstenergy=500,
            nstlog=500,
            nstxout_compressed=500,
            continuation="no",
            tcoupl="V-rescale",
            tc_grps=f"Protein_{ligid} Water_and_ions",
            tau_t=0.1,
            pcoupl="no",
            gen_vel="yes",
        )

        return mdp

    elif simulation_type == "nvt" and pdb_type == "PROTEIN":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=500000,
            # nsteps=50000,
            nstenergy=500,
            nstlog=500,
            nstxout_compressed=500,
            continuation="no",
            tcoupl="V-rescale",
            tc_grps="Protein Water_and_ions",
            tau_t=0.1,
            pcoupl="no",
            gen_vel="yes",
        )

        return mdp

    elif simulation_type == "nvt" and pdb_type == "LIGAND":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=500000,
            # nsteps=50000,
            nstenergy=500,
            nstlog=500,
            nstxout_compressed=500,
            continuation="no",
            tcoupl="V-rescale",
            tc_grps="System",
            tau_t=0.1,
            pcoupl="no",
            gen_vel="yes",
        )

        return mdp

    elif simulation_type == "npt" and pdb_type == "COMPLEX":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=2500000,
            # nsteps=50000,
            nstenergy=500,
            nstlog=500,
            nstxout_compressed=500,
            tc_grps=f"Protein_{ligid} Water_and_ions",
            pcoupl="C-rescale",
        )

        return mdp

    elif simulation_type == "npt" and pdb_type == "PROTEIN":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=2500000,
            # nsteps=50000,
            nstenergy=500,
            nstlog=500,
            nstxout_compressed=500,
            tc_grps="Protein Water_and_ions",
            pcoupl="C-rescale",
        )

        return mdp

    elif simulation_type == "npt" and pdb_type == "LIGAND":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=2500000,
            # nsteps=50000,
            nstenergy=500,
            nstlog=500,
            nstxout_compressed=500,
            tc_grps="System",
            pcoupl="C-rescale",
        )

        return mdp

    elif simulation_type == "md" and pdb_type == "COMPLEX":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=5000000,
            # nsteps=100000,
            tc_grps=f"Protein_{ligid} Water_and_ions",
        )

        return mdp

    elif simulation_type == "md" and pdb_type == "PROTEIN":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=500000000,
            # nsteps=100000,
            tc_grps="Protein Water_and_ions",
            # freezegrps="Protein",
            # freezedim="Y Y Y",
        )

        return mdp

    elif simulation_type == "md" and pdb_type == "LIGAND":
        mdp = MDP(
            sim_type=simulation_type,
            nsteps=5000000,
            # nsteps=100000,
            tc_grps="System",
        )

        return mdp


def calc_freeBE(production_md_output: Path, fit: Any, root_dir: str):
    # Calculate Binding Free Energy
    Path(f"{root_dir}/energy").mkdir(parents=True, exist_ok=True)

    echo = subprocess.run(
        "(echo 'q')",
        shell=True,
        check=True,
        capture_output=True,
    )

    gen_index = Commandline_Operation(
        params=["make_ndx"],
        input_files={"-f": Path(f"{root_dir}/md.gro")},
        output_files={"-o": Path(f"{root_dir}/all_index.ndx")},
        input_for_subprocess=echo.stdout.decode("utf-8"),
    )

    gen_index.run_gmx_cmdline()

    cmd = (
        "gmx_MMPBSA --clean | "
        f"gmx_MMPBSA -O -i {Path(root_dir).parent}/src/mmpbsa.in -cs {production_md_output}"
        f" -ci {str(gen_index.output_files['-o'])}"
        f" -cg 1 13 -ct {fit.output_files['-o']} -cp {root_dir}/topol.top "
        f"-o {root_dir}/energy/FINAL_RESULTS_MMPBSA.dat -eo {root_dir}/energy/FINAL_RESULTS_MMPBSA.csv -nogui"
    )
    subprocess.run(cmd, shell=True)


def main(pdb_file: Path, gpu_id: int, output_dir: Path) -> Any:
    root_dir: str = f"{output_dir}/{pdb_file.name.split('_')[0]}"

    # Check PDB
    # If PDB is a protein-ligand/complex file, split them up into protein and ligand PDB files separately.
    # Otherwise, remain as it is.
    with open(pdb_file, "r") as pdb_input:
        pdb_contents: list[str] = [
            line.strip("\n") for line in pdb_input.readlines()
        ] | delete_waters

    # Check PDB for missing information
    pdb_type: str = check_pdb_type(pdb_contents)

    problematic_lines = check_missing_pdb_info(pdb_contents)
    if len(problematic_lines) > 0:
        print("Missing Data/Data in Wrong Column in PDB File:\n")
        for key, value in problematic_lines.items():
            print(f"{key}\t(Missing: {value})")
        sys.exit("\nPlease check your PDB! Exiting...")

    pdb_type_dict: Dict[str, list[str]] = {}
    protein_pdb_contents: list[str] = []
    ligand_pdb_contents: list[str] = []
    if pdb_type == "COMPLEX":
        pdb_type_dict["COMPLEX"] = pdb_contents

        ligand_atom_numbers: list[int] = []
        for line in pdb_contents:
            pdb_dict_terms = get_pdb_line_components(line)
            if (
                line.startswith("ATOM")
                or line.startswith("HETATM")
                and str(pdb_dict_terms["RESIDUE NAME"])
                not in amino_acids_3_letter_codes
            ):
                ligand_atom_numbers.append(int(pdb_dict_terms["ATOM SERIAL NUMBER"]))

            if (
                line.startswith("CONECT")
                and int(pdb_dict_terms["ATOM SERIAL NUMBER"]) in ligand_atom_numbers
            ):
                ligand_pdb_contents.append(line)
            elif (
                line.startswith("ATOM")
                and str(pdb_dict_terms["RESIDUE NAME"]) in amino_acids_3_letter_codes
            ):
                protein_pdb_contents.append(line)
            elif (
                line.startswith("ATOM")
                or line.startswith("HETATM")
                and str(pdb_dict_terms["RESIDUE NAME"])
                not in amino_acids_3_letter_codes
            ):
                ligand_pdb_contents.append(line)

    elif pdb_type == "LIGAND":
        pdb_type_dict["COMPLEX"] = []

        ligand_atom_numbers: list[int] = []
        for line in pdb_contents:
            ligand_pdb_contents.append(line)

    elif pdb_type == "PROTEIN":
        pdb_type_dict["COMPLEX"] = []
        protein_pdb_contents = pdb_contents
    pdb_type_dict["PROTEIN"] = protein_pdb_contents
    pdb_type_dict["LIGAND"] = ligand_pdb_contents
    if ligand_pdb_contents:
        ligand_id: str = str(
            (get_pdb_line_components(ligand_pdb_contents[0]))["RESIDUE NAME"]
        )
    else:
        ligand_id = ""

    for system_type, pdb_contents in pdb_type_dict.items():
        if len(pdb_contents) != 0:
            with open(f"{root_dir}/{system_type}_cleaned.pdb", "w") as cleaned_pdb_file:
                for line in clean_pdb(system_type, pdb_contents):
                    cleaned_pdb_file.write(line + "\n")

    # Generate various GROMACS required files for MD
    if len(pdb_type_dict["PROTEIN"]) > 0:
        generate_prot_top("PROTEIN", root_dir)

    # Process and Generate Ligand Files
    if len(pdb_type_dict["LIGAND"]) > 0:
        hydrogenate_lig(root_dir=root_dir, ligand_id=ligand_id)
        fix_lig_bond_order(root_dir=root_dir, ligand_id=ligand_id)
        generate_lig_files(
            root_dir=root_dir,
            ligand_pdb_contents=ligand_pdb_contents,
            ligand_id=ligand_id,
        )

    # Combine {protein}_processed.gro and {ligand}.gro
    print(pdb_type)
    if pdb_type == "COMPLEX":
        with open(f"{root_dir}/complex.gro", "w") as complex_f:
            for line in combine_prot_lig_gro(root_dir=root_dir, ligand_id=ligand_id):
                complex_f.write(f"{line}")
    elif pdb_type == "LIGAND":
        shutil.copy(
            Path(f"{root_dir}/{ligand_id}_fix.acpype/{ligand_id}_fix_GMX.gro"),
            Path(f"{root_dir}/complex.gro"),
        )

        shutil.copy(
            Path(f"{root_dir}/{ligand_id}_fix.acpype/{ligand_id}_fix_GMX.top"),
            Path(f"{root_dir}/topol.top"),
        )
        shutil.copy(
            Path(f"{root_dir}/{ligand_id}_fix.acpype/{ligand_id}_fix_GMX.itp"),
            Path(f"{root_dir}"),
        )
        shutil.copy(
            Path(f"{root_dir}/{ligand_id}_fix.acpype/posre_{ligand_id}_fix.itp"),
            Path(f"{root_dir}"),
        )

    elif pdb_type == "PROTEIN":
        shutil.copy(
            Path(f"{root_dir}/{pdb_type}_processed.gro"),
            Path(f"{root_dir}/complex.gro"),
        )
    edit_topol_top(pdb_type, Path(f"{root_dir}/topol.top"), ligand_id)  # Edit topol.top

    solvation_output = Solvation(
        Path(f"{root_dir}/complex.gro"), root_dir
    ).solvate_system()  # Solvate the Complex

    # Create MDP files
    for sim_type in ["ions", "em", "nvt", "npt", "md"]:
        mdp = mdp_type(
            simulation_type=sim_type,
            pdb_type=pdb_type,
            ligid=ligand_id,
        )

        write_mdp = mdp.write_mdp_file()
        with open(Path(f"{root_dir}/{sim_type}.mdp"), "w") as sim_mdp_file:
            for key, value in write_mdp.items():
                if mdp.tc_grps != "System" and key in ["tau_t", "ref_t"]:
                    sim_mdp_file.write(f"{key} = {value} {value}\n")
                else:
                    sim_mdp_file.write(f"{key} = {value}\n")

    # Add ions
    ions_output = AddIons(
        ions_mdp_file=Path(f"{root_dir}/ions.mdp"),
        solvated_gro_file=solvation_output[0].output_files["-o"],
        pdb_type=pdb_type,
        root_dir=root_dir,
    ).add_ions_to_system()

    edit_mdp_files(pdb_type, Path(f"{root_dir}/topol.top"), ligand_id)  # Edit MDP files

    # Energy minimzation
    Minimization(
        minimization_mdp_file=Path(f"{root_dir}/em.mdp"),
        solvated_gro_file_with_ions=ions_output[0].output_files["-o"],
        root_dir=root_dir,
        # additional_gmx_args=additional_gmx_args,
    ).preprocess()

    subprocess.run(
        (
            f"OMP_NUM_THREADS=5 GMX_GPU_ID={gpu_id} nixGL gmx mdrun -v -pin on -pinoffset {gpu_id * 16}"
            f" -ntmpi 3 -deffnm em"
        ),
        shell=True,
    )

    # NVT Equilibration
    NVT_Equilibration(
        minimization_gro_file=Path("em.gro"),
        nvt_mdp_file=Path("nvt.mdp"),
        pdb_type=pdb_type,
        ligid=ligand_id,
        root_dir=root_dir,
        # additional_gmx_args=additional_gmx_args,
    ).preprocess()

    subprocess.run(
        (
            f"OMP_NUM_THREADS=5 GMX_GPU_ID={gpu_id} nixGL gmx mdrun -v -pin on"
            f" -pinoffset {gpu_id * 16} -ntmpi 3 -nb gpu -pme gpu -npme 1 -deffnm nvt"
        ),
        shell=True,
    )

    # NPT Equilibration
    NPT_Equilibration(
        nvt_gro_file=Path(f"{root_dir}/nvt.gro"),
        npt_mdp_file=Path(f"{root_dir}/npt.mdp"),
        pdb_type=pdb_type,
        ligid=ligand_id,
        root_dir=root_dir,
        # additional_gmx_args=additional_gmx_args,
    ).preprocess()

    subprocess.run(
        (
            f"OMP_NUM_THREADS=5 GMX_GPU_ID={gpu_id} nixGL gmx mdrun -v -pin on"
            f" -pinoffset {gpu_id * 16} -ntmpi 3 -nb gpu -pme gpu -npme 1 -deffnm npt"
        ),
        shell=True,
    )

    # Production MD
    Production_MD(
        npt_gro_file=Path(f"{root_dir}/npt.gro"),
        md_mdp_file=Path(f"{root_dir}/md.mdp"),
        pdb_type=pdb_type,
        ligid=ligand_id,
        root_dir=root_dir,
        # additional_gmx_args=additional_gmx_args,
    ).preprocess()

    subprocess.run(
        (
            f"OMP_NUM_THREADS=5 GMX_GPU_ID={gpu_id} nixGL gmx mdrun -v -pin on"
            f" -pinoffset {gpu_id * 16} -ntmpi 3 -nb gpu -pme gpu -npme 1 -deffnm md"
        ),
        shell=True,
    )

    # Post-processing Trajectory
    no_pbc = remove_pbc(
        Path(f"{root_dir}/md.tpr"), Path(f"{root_dir}/md.xtc"), root_dir
    )
    fit = rot_trans_fit(no_pbc.input_files["-s"], no_pbc.output_files["-o"], root_dir)

    if pdb_type == "COMPLEX":
        calc_freeBE(
            Path(f"{root_dir}/md.tpr"), fit, root_dir
        )  # Calculate Binding Free Energy

    # Move all output files into directory
    # datetime_now = (
    #     str(datetime.now()).replace(" ", "-").replace(":", "-").replace(".", "-")
    # )
    Path(f"./{pdb_file.stem}/{pdb_type}").mkdir(parents=True, exist_ok=True)
    subprocess.run(
        (
            f"mv {pdb_file} *_cleaned.pdb gmx_MMPBSA.log md* *.itp UNL* *.top \\#* *.gro em*"
            f" *.ndx ions.* npt* nvt* LIG* complex.pdb.diff_chain complex.pdb.same_chain"
            f" ligand.pdb test_complex.pdb _GMXMMPBSA* *.mmxsa COM* REC* energy {ligand_id}*"
            f" ./{pdb_file.stem}/{pdb_type}"
        ),
        shell=True,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--pdb_file", type=Path, required=True)
    # parser.add_argument("--additional_gmx_args", type=str, default="")
    parser.add_argument("--gpu_id", type=int, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)

    args = parser.parse_args()

    main(
        pdb_file=args.pdb_file,
        # additional_gmx_args=args.additional_gmx_args,
        gpu_id=args.gpu_id,
        output_dir=args.output_dir,
    )
