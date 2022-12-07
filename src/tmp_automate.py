import argparse
import subprocess

# import time
from pathlib import Path
from typing import Any, List


def main(
    pdb_dir: Path, output_dir: Path, prefix_list: List[str], ligand_type: str
) -> Any:
    root = str(pdb_dir)
    for gpu_id, prefix in enumerate(prefix_list):
        subprocess.Popen(
            (
                f"mkdir {output_dir}/{prefix.split('_')[0]} && cd {output_dir}/{prefix.split('_')[0]}"
                " && python /home/ubuntu/automated_gromacs_md/src/main.py"
                f" --pdb_file {root}/{prefix}_{ligand_type}_complex.pdb --gpu_id={gpu_id}"
                f" --output_dir={output_dir}"
            ),
            shell=True,
        )
        # time.sleep(30)
        # subprocess.run(
        #     f"cp LIGAND_cleaned.pdb ./{prefix}_ligand.pdb",
        #     shell=True,
        # )
        # subprocess.run(f"python src/main.py --pdb_file {prefix}_ligand.pdb", shell=True)

        # if args.ligand_type == "active":
        #     subprocess.run(
        #         f"cp PROTEIN_cleaned.pdb ./{prefix}_protein.pdb",
        #         shell=True,
        #     )
        # subprocess.run(
        #     (
        #         f"python src/main.py --pdb_file {prefix}_protein.pdb"
        #         f" --additional_gmx_args={additional_gmx_args}"
        #     ),
        #     shell=True,
        # )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--pdb_dir", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--prefix_list", type=str, nargs="*", required=True)
    parser.add_argument(
        "--ligand_type", type=str, choices=["active", "inactive"], required=True
    )
    # parser.add_argument("--additional_gmx_args", type=str, default="")

    args = parser.parse_args()

    main(
        pdb_dir=args.pdb_dir,
        output_dir=args.output_dir,
        prefix_list=args.prefix_list,
        ligand_type=args.ligand_type,
        # additional_gmx_args=args.additional_gmx_args,
    )
