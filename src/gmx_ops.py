import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

from utils import Commandline_Operation, calc_net_charge, read_topol


def edit_mdp_files(pdb_type: str, topol_file: Path, ligid: str) -> None:
    topol_contents: Dict[str, list[str]] = read_topol(pdb_type, topol_file, ligid)
    topol_lines_containing_ions: list[str] = []
    for line in topol_contents["MOLECULES"]:
        if "NA" in line or "CL" in line:
            topol_lines_containing_ions.append(line)

    for file in os.listdir("."):
        if len(topol_lines_containing_ions) == 0:
            if Path(file).suffix == ".mdp":
                subprocess.run(f"sed -i 's/Water_and_ions/Water/g' {file}", shell=True)


def generate_prot_top(pdb_type: str, root_dir: str) -> subprocess.CompletedProcess[str]:
    echo = subprocess.run(
        "(echo '5' && echo '1')",
        shell=True,
        check=True,
        capture_output=True,
    )

    gen_prot_top = Commandline_Operation(
        params=["pdb2gmx", "-posrefc", "1000", "-merge", "all"],
        input_files={"-f": Path(f"{root_dir}/{pdb_type}_cleaned.pdb")},
        output_files={"-o": Path(f"{root_dir}/{pdb_type}_processed.gro")},
        input_for_subprocess=echo.stdout.decode("utf-8"),
    )

    return gen_prot_top.run_gmx_cmdline()


def hydrogenate_lig(
    ligand_id: str, root_dir: str, lig_pdb_path: Path = Path("LIGAND_cleaned.pdb")
) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        f"obabel -ipdb {root_dir}/{lig_pdb_path} -O {root_dir}/{ligand_id}.mol2 -h",
        shell=True,
    )


def fix_lig_bond_order(
    ligand_id: str, root_dir: str
) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        (
            f"perl {Path(root_dir).parent}/src/sort_mol2_bonds.pl {root_dir}/{ligand_id}.mol2"
            f" {root_dir}/{ligand_id}_fix.mol2"
        ),
        shell=True,
    )


def generate_lig_files(
    ligand_pdb_contents: list[str], ligand_id: str, root_dir: str
) -> subprocess.CompletedProcess[bytes]:
    net_charge = calc_net_charge(ligand_pdb_contents)
    return subprocess.run(
        f"acpype -i {root_dir}/{ligand_id}_fix.mol2 -n {net_charge} -f", shell=True
    )


def combine_prot_lig_gro(
    root_dir: str, ligand_id: str, protein_pdb: str = "PROTEIN"
) -> list[str]:

    complex_lines: list[str] = []
    with open(f"{root_dir}/{protein_pdb}_processed.gro", "r") as protein_f:
        protein_lines = protein_f.readlines()
    with open(
        f"{root_dir}/{ligand_id}_fix.acpype/{ligand_id}_fix_GMX.gro", "r"
    ) as ligand_f:
        ligand_lines = ligand_f.readlines()
    for line in protein_lines[2:-1]:
        complex_lines.append(line)
    for line in ligand_lines[2:-1]:
        complex_lines.append(line)
    complex_lines.insert(0, protein_lines[0])
    complex_lines.insert(1, f"{len(complex_lines) - 1}\n")
    complex_lines.append(protein_lines[-1])

    return complex_lines


def edit_topol_top(pdb_type: str, topol_file: Path, ligid: str):
    # Edit topol.top
    edited_topol_contents = read_topol(pdb_type, topol_file, ligid)
    with open("tmp.top", "w") as tmp_topol_file:
        topol_dict_keys = [
            "FORCEFIELD PARAMS",
            # "DEFAULTS",
            "MOLECULE TYPE",
            "ATOMS",
            "BONDS",
            "PAIRS",
            "ANGLES",
            "DIHEDRALS1",
            "DIHEDRALS2",
            "CMAP",
            "POSRE",
            "CHAIN TOPOLOGIES",
            "WATER TOPOLOGY",
            "WATER POSRE",
            "IONS TOPOLOGY",
            "SYSTEM",
            "MOLECULES",
        ]
        if ligid != "":
            topol_dict_keys.insert(11, "LIGAND POSRE")
            topol_dict_keys.insert(1, "LIGAND TOPOLOGY")

        for key in topol_dict_keys:
            for line in edited_topol_contents[key]:
                tmp_topol_file.write(f"{line}\n")
            if key != "MOLECULES":
                tmp_topol_file.write("\n")

    shutil.move(Path("tmp.top"), topol_file)


@dataclass
class MDP:
    # Type of Simulation
    sim_type: str = "md"
    # Preprocessing
    define: str = "-DPOSRES"
    # Run control
    integrator: str = "md"
    dt: float = 0.002
    nsteps: int = 500000
    # Energy minimization
    emtol: float = 1000.0
    emstep: float = 0.01
    # Output control
    nstlog: int = 5000
    nstenergy: int = 5000
    nstxout_compressed: int = 5000
    # Neighbour searching
    cutoff_scheme: str = "Verlet"
    nstlist: int = 20
    pbc: str = "xyz"
    ns_type: str = "grid"
    rlist: float = 1.2
    # Electrostatics
    coulombtype: str = "PME"
    rcoulomb: float = 1.2
    # Van der Waals
    vdwtype: str = "cutoff"
    vdw_modifier: str = "force-switch"
    rvdw_switch: float = 1.0
    rvdw: float = 1.2
    DispCorr: str = "no"
    # Ewald
    fourierspacing: float = 0.16
    pme_order: int = 4
    # Temperature coupling
    tcoupl: str = "nose-hoover"
    tc_grps: str = "Protein Non-Protein"
    tau_t: float = 0.4
    ref_t: float = 300
    # Pressure coupling
    pcoupl: str = "Parrinello-Rahman"
    pcoupltype: str = "isotropic"
    tau_p: float = 2.0
    compressibility: float = 4.5e-5
    ref_p: float = 1.0
    refcoord_scaling: str = "com"
    # Velocity generation
    gen_vel: str = "no"
    gen_temp: float = 300
    gen_seed: int = -1
    # Bonds
    constraints: str = "h-bonds"
    constraint_algorithm: str = "LINCS"
    continuation: str = "yes"
    lincs_order: int = 4
    lincs_iter: int = 1
    # Non-equilibrium MD
    freezegrps: str = ""
    freezedim: str = ""

    def write_mdp_file(self) -> Dict[str, Any]:

        sim_types = ["ions", "em", "nvt", "npt", "md"]

        default_params_for_ions_and_em: Dict[str, Any] = {
            "integrator": self.integrator,
            "nsteps": self.nsteps,
            "emtol": self.emtol,
            "emstep": self.emstep,
            "cutoff_scheme": self.cutoff_scheme,
            "nstlist": self.nstlist,
            "pbc": self.pbc,
            "ns_type": self.ns_type,
            "rlist": self.rlist,
            "coulombtype": self.coulombtype,
            "rcoulomb": self.rcoulomb,
            "rvdw": self.rvdw,
        }

        default_params_for_nvt_npt_md: Dict[str, Any] = {
            "integrator": self.integrator,
            "dt": self.dt,
            "nsteps": self.nsteps,
            "nstlog": self.nstlog,
            "nstenergy": self.nstenergy,
            "nstxout_compressed": self.nstxout_compressed,
            "cutoff_scheme": self.cutoff_scheme,
            "nstlist": self.nstlist,
            "pbc": self.pbc,
            "ns_type": self.ns_type,
            "rlist": self.rlist,
            "coulombtype": self.coulombtype,
            "rcoulomb": self.rcoulomb,
            "vdwtype": self.vdwtype,
            "vdw_modifier": self.vdw_modifier,
            "rvdw_switch": self.rvdw_switch,
            "rvdw": self.rvdw,
            "DispCorr": self.DispCorr,
            "fourierspacing": self.fourierspacing,
            "pme_order": self.pme_order,
            "tcoupl": self.tcoupl,
            "tc_grps": self.tc_grps,
            "tau_t": self.tau_t,
            "ref_t": self.ref_t,
            "pcoupl": self.pcoupl,
            "gen_vel": self.gen_vel,
            "constraints": self.constraints,
            "constraint_algorithm": self.constraint_algorithm,
            "continuation": self.continuation,
            "lincs_order": self.lincs_order,
            "lincs_iter": self.lincs_iter,
        }

        if self.sim_type in sim_types and self.sim_type == "ions":
            return default_params_for_ions_and_em

        elif self.sim_type in sim_types and self.sim_type == "em":

            additional_params: Dict[str, Any] = {
                "vdwtype": self.vdwtype,
                "vdw_modifier": self.vdw_modifier,
                "rvdw_switch": self.rvdw_switch,
                "DispCorr": self.DispCorr,
            }

            return default_params_for_ions_and_em | additional_params

        elif self.sim_type in sim_types and self.sim_type == "nvt":
            additional_params: Dict[str, Any] = {
                "define": self.define,
                "gen_temp": self.gen_temp,
                "gen_seed": self.gen_seed,
            }

            return default_params_for_nvt_npt_md | additional_params

        elif self.sim_type in sim_types and self.sim_type == "npt":
            additional_params: Dict[str, Any] = {
                "define": self.define,
                "pcoupltype": self.pcoupltype,
                "tau_p": self.tau_p,
                "compressibility": self.compressibility,
                "ref_p": self.ref_p,
                "refcoord_scaling": self.refcoord_scaling,
            }

            return default_params_for_nvt_npt_md | additional_params

        elif self.sim_type in sim_types and self.sim_type == "md":
            additional_params: Dict[str, Any] = {
                "pcoupltype": self.pcoupltype,
                "tau_p": self.tau_p,
                "compressibility": self.compressibility,
                "ref_p": self.ref_p,
                # "freezegrps": self.freezegrps,
                # "freezedim": self.freezedim,
            }

            return default_params_for_nvt_npt_md | additional_params

        return {}
