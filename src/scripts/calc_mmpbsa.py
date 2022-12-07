import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


@dataclass
class Commandline_Operation:
    params: List[str]
    input_files: Dict[str, Path]
    output_files: Dict[str, Path]
    input_for_subprocess: str = ""

    def run_gmx_cmdline(self) -> subprocess.CompletedProcess[str]:
        self.params.insert(0, "gmx")

        for file_group in [self.input_files, self.output_files]:
            for key, value in file_group.items():
                self.params.append(key)
                self.params.append(str(value))

        return subprocess.run(
            self.params, input=self.input_for_subprocess, encoding="utf-8"
        )


def calc_freeBE():
    # Calculate Binding Free Energy
    Path("energy").mkdir(parents=True, exist_ok=True)

    echo = subprocess.run(
        "(echo 'q')",
        shell=True,
        check=True,
        capture_output=True,
    )

    gen_index = Commandline_Operation(
        params=["make_ndx"],
        input_files={"-f": Path("md.gro")},
        output_files={"-o": Path("all_index.ndx")},
        input_for_subprocess=echo.stdout.decode("utf-8"),
    )

    gen_index.run_gmx_cmdline()

    cmd = (
        "gmx_MMPBSA --clean | gmx_MMPBSA -O -i mmpbsa.in -cs md.tpr -ci all_index.ndx -cg 1 13"
        " -ct md_fit.xtc -cp topol.top -o energy/FINAL_RESULTS_MMPBSA.dat"
        " -eo energy/FINAL_RESULTS_MMPBSA.csv -nogui"
    )
    subprocess.run(cmd, shell=True)


if __name__ == "__main__":
    calc_freeBE()
