import subprocess
from pathlib import Path
from typing import Any

from utils import Commandline_Operation


def remove_pbc(tpr_input: Path, xtc_input: Path, root_dir: str) -> Any:
    # Remove PBC
    echo = subprocess.run(
        "(echo '1' && echo '0')",
        shell=True,
        check=True,
        capture_output=True,
    )

    remove_pbc = Commandline_Operation(
        params=["trjconv", "-pbc", "nojump", "-center"],
        input_files={
            "-s": tpr_input,
            "-f": xtc_input,
        },
        output_files={"-o": Path(f"{root_dir}/md_noPBC.xtc")},
        input_for_subprocess=echo.stdout.decode("utf-8"),
    )
    remove_pbc.run_gmx_cmdline()

    return remove_pbc


def rot_trans_fit(tpr_input: Path, xtc_input: Path, root_dir: str) -> Any:
    # Perform Rotational and Translational Fitting
    echo = subprocess.run(
        "(echo '4' && echo '0')",
        shell=True,
        check=True,
        capture_output=True,
    )

    fit = Commandline_Operation(
        params=["trjconv", "-fit", "rot+trans"],
        input_files={
            "-s": tpr_input,
            "-f": xtc_input,
        },
        output_files={"-o": Path(f"{root_dir}/md_fit.xtc")},
        input_for_subprocess=echo.stdout.decode("utf-8"),
    )
    fit.run_gmx_cmdline()

    return fit
