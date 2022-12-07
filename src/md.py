import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Tuple

from utils import Commandline_Operation


@dataclass
class Solvation:
    complex_file: Path
    root_dir: str

    def generate_solvent_box(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        gen_solv_box = Commandline_Operation(
            params=["editconf", "-bt", "cubic", "-c", "-d", "1.0"],
            input_files={"-f": self.complex_file},
            output_files={"-o": Path(f"{self.root_dir}/newbox.gro")},
        )

        return (gen_solv_box, gen_solv_box.run_gmx_cmdline())

    def solvate_system(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        solvent_box = self.generate_solvent_box()

        run_solvation = Commandline_Operation(
            params=["solvate"],
            input_files={
                "-cp": Path(solvent_box[0].output_files["-o"]),
                "-cs": Path("spc216.gro"),
            },
            output_files={
                "-p": Path(f"{self.root_dir}/topol.top"),
                "-o": Path(f"{self.root_dir}/solv.gro"),
            },
        )

        return (run_solvation, run_solvation.run_gmx_cmdline())


@dataclass
class AddIons:
    ions_mdp_file: Path
    solvated_gro_file: Path
    pdb_type: str
    root_dir: str

    def preprocess(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        gen_preproc_file = Commandline_Operation(
            params=["grompp"],
            input_files={
                "-f": self.ions_mdp_file,
                "-c": self.solvated_gro_file,
                "-p": Path(f"{self.root_dir}/topol.top"),
            },
            output_files={"-o": Path(f"{self.root_dir}/ions.tpr")},
        )

        return (gen_preproc_file, gen_preproc_file.run_gmx_cmdline())

    def add_ions_to_system(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        ions = self.preprocess()

        if self.pdb_type == "COMPLEX":
            echo = subprocess.run(
                "(echo '15')",
                shell=True,
                check=True,
                capture_output=True,
            )
        elif self.pdb_type == "PROTEIN":
            echo = subprocess.run(
                "(echo '13')",
                shell=True,
                check=True,
                capture_output=True,
            )
        else:
            echo = subprocess.run(
                "(echo '4')",
                shell=True,
                check=True,
                capture_output=True,
            )

        run_adding_ions = Commandline_Operation(
            params=["genion", "-pname", "NA", "-nname", "CL", "-neutral"],
            input_files={
                "-s": Path(ions[0].output_files["-o"]),
            },
            output_files={
                "-p": Path(f"{self.root_dir}/topol.top"),
                "-o": Path(f"{self.root_dir}/solv_ions.gro"),
            },
            input_for_subprocess=echo.stdout.decode("utf-8"),
        )

        return (run_adding_ions, run_adding_ions.run_gmx_cmdline())


@dataclass
class Minimization:
    minimization_mdp_file: Path
    solvated_gro_file_with_ions: Path
    root_dir: str
    # additional_gmx_args: str

    def preprocess(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        gen_preproc_file = Commandline_Operation(
            params=["grompp", "-maxwarn", "99"],
            input_files={
                "-f": self.minimization_mdp_file,
                "-c": self.solvated_gro_file_with_ions,
                "-p": Path(f"{self.root_dir}/topol.top"),
            },
            output_files={"-o": Path(f"{self.root_dir}/em.tpr")},
        )

        return (gen_preproc_file, gen_preproc_file.run_gmx_cmdline())

    def run_minimization(self) -> subprocess.CompletedProcess[str]:
        self.preprocess()

        run_minimization = Commandline_Operation(
            params=[
                "mdrun",
                "-v",
                "-deffnm",
                "em",
                "-pin",
                "on",
                "-pinoffset",
                # f"{int(((self.additional_gmx_args).split('='))[-1]) * 16}",
                "-ntmpi",
                "3",
            ],
            input_files={},
            output_files={},
            # additional_gmx_args=self.additional_gmx_args,
        )

        return run_minimization.run_gmx_cmdline()


@dataclass
class NVT_Equilibration:
    pdb_type: str
    minimization_gro_file: Path
    nvt_mdp_file: Path
    ligid: str
    root_dir: str
    # additional_gmx_args: str

    def gen_complex_index_file(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        echo = subprocess.run(
            "(echo '1 | 13' && echo 'q')",
            shell=True,
            check=True,
            capture_output=True,
        )

        index_file = Commandline_Operation(
            params=["make_ndx"],
            input_files={
                "-f": self.minimization_gro_file,
            },
            output_files={"-o": Path(f"{self.root_dir}/index.ndx")},
            input_for_subprocess=echo.stdout.decode("utf-8"),
        )

        return (index_file, index_file.run_gmx_cmdline())

    def preprocess(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:

        input_files: Dict[str, Path] = {
            "-f": self.nvt_mdp_file,
            "-c": self.minimization_gro_file,
            "-r": self.minimization_gro_file,
            "-p": Path(f"{self.root_dir}/topol.top"),
        }

        if self.pdb_type == "COMPLEX":
            complex_index_file = self.gen_complex_index_file()
            input_files["-n"] = Path(complex_index_file[0].output_files["-o"])

        gen_preproc_file = Commandline_Operation(
            params=["grompp", "-maxwarn", "99"],
            input_files=input_files,
            output_files={"-o": Path(f"{self.root_dir}/nvt.tpr")},
        )

        return (gen_preproc_file, gen_preproc_file.run_gmx_cmdline())

    def run_nvt_equil(self) -> subprocess.CompletedProcess[str]:
        self.preprocess()

        run_nvt_equilibration = Commandline_Operation(
            params=[
                "mdrun",
                "-v",
                "-deffnm",
                "nvt",
                "-pin",
                "on",
                "-pinoffset",
                # f"{int(((self.additional_gmx_args).split('='))[-1]) * 16}",
                "-ntmpi",
                "3",
                "-nb",
                "gpu",
                "-pme",
                "gpu",
                "-npme",
                "1",
            ],
            input_files={},
            output_files={},
            # additional_gmx_args=self.additional_gmx_args,
        )

        return run_nvt_equilibration.run_gmx_cmdline()


@dataclass
class NPT_Equilibration:
    pdb_type: str
    nvt_gro_file: Path
    npt_mdp_file: Path
    ligid: str
    root_dir: str
    # additional_gmx_args: str

    def preprocess(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        input_files: Dict[str, Path] = {
            "-f": self.npt_mdp_file,
            "-c": self.nvt_gro_file,
            "-r": self.nvt_gro_file,
            "-p": Path(f"{self.root_dir}/topol.top"),
        }

        if self.pdb_type == "COMPLEX":
            input_files["-n"] = Path(f"{self.root_dir}/index.ndx")

        gen_preproc_file = Commandline_Operation(
            params=["grompp", "-maxwarn", "99"],
            input_files=input_files,
            output_files={"-o": Path(f"{self.root_dir}/npt.tpr")},
        )

        return (gen_preproc_file, gen_preproc_file.run_gmx_cmdline())

    def run_npt_equil(self) -> subprocess.CompletedProcess[str]:
        self.preprocess()

        run_npt_equilibration = Commandline_Operation(
            params=[
                "mdrun",
                "-v",
                "-deffnm",
                "npt",
                "-pin",
                "on",
                "-pinoffset",
                # f"{int(((self.additional_gmx_args).split('='))[-1]) * 16}",
                "-ntmpi",
                "3",
                "-nb",
                "gpu",
                "-pme",
                "gpu",
                "-npme",
                "1",
            ],
            input_files={},
            output_files={},
            # additional_gmx_args=self.additional_gmx_args,
        )

        return run_npt_equilibration.run_gmx_cmdline()


@dataclass
class Production_MD:
    pdb_type: str
    npt_gro_file: Path
    md_mdp_file: Path
    ligid: str
    root_dir: str
    # additional_gmx_args: str

    def preprocess(self) -> Tuple[Any, subprocess.CompletedProcess[str]]:
        input_files: Dict[str, Path] = {
            "-f": self.md_mdp_file,
            "-c": self.npt_gro_file,
            "-r": self.npt_gro_file,
            "-p": Path(f"{self.root_dir}/topol.top"),
        }

        if self.pdb_type == "COMPLEX":
            input_files["-n"] = Path(f"{self.root_dir}/index.ndx")

        gen_preproc_file = Commandline_Operation(
            params=["grompp", "-maxwarn", "99"],
            input_files=input_files,
            output_files={"-o": Path(f"{self.root_dir}/md.tpr")},
        )

        return (gen_preproc_file, gen_preproc_file.run_gmx_cmdline())

    def run_prod_md(self) -> Tuple[Path, subprocess.CompletedProcess[str]]:
        prod_md_preproc = self.preprocess()

        run_production_md = Commandline_Operation(
            params=[
                "mdrun",
                "-v",
                "-deffnm",
                "md",
                "-pin",
                "on",
                "-pinoffset",
                # f"{int(((self.additional_gmx_args).split('='))[-1]) * 16}",
                "-ntmpi",
                "3",
                "-nb",
                "gpu",
                "-pme",
                "gpu",
                "-npme",
                "1",
            ],
            input_files={},
            output_files={},
            # additional_gmx_args=self.additional_gmx_args,
        )

        return (
            prod_md_preproc[0].output_files["-o"],
            run_production_md.run_gmx_cmdline(),
        )
