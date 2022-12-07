from simple_term_menu import TerminalMenu  # type: ignore
from .pdb_utils import download_pdb
from typing import List, Optional
from subprocess import run
import sys


def main():
    while True:
        terminal_menu = TerminalMenu(
            [
                "Download PDB",
                "Return to Main Menu",
                "Quit",
            ],
            multi_select=False,
        )
        terminal_menu.show()

        if terminal_menu.chosen_menu_entry == "Download PDB":
            pdb_id = input("Please enter the 4-letter PDB ID: ")
            with open(f"{pdb_id.upper()}.pdb", "w") as pdb_file:
                pdb_contents: Optional[List[str]] = download_pdb(pdb_id)
                if pdb_contents is not None:
                    for line in pdb_contents:
                        pdb_file.write(f"{line}\n")

            print(f"{pdb_id.upper()}.pdb downloaded!\n")

        elif terminal_menu.chosen_menu_entry == "Return to Main Menu":
            run("python main.py", shell=True)

        elif terminal_menu.chosen_menu_entry == "Quit":
            sys.exit("Goodbye!")


if __name__ == "__main__":
    main()
