from simple_term_menu import TerminalMenu  # type: ignore
from protein_preparation_wizard import main as ppw_menu
import sys


def main():
    terminal_menu = TerminalMenu(
        [
            "Protein Preparation Wizard",
            "Generate MDP Files",
            "Run Molecular Dynamics Simulation",
            "Quit",
        ],
        multi_select=False,
    )
    terminal_menu.show()

    if terminal_menu.chosen_menu_entry == "Protein Preparation Wizard":
        ppw_menu()

    elif terminal_menu.chosen_menu_entry == "Quit":
        sys.exit("Goodbye!")


if __name__ == "__main__":
    main()
