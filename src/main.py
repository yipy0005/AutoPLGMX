from simple_term_menu import TerminalMenu  # type: ignore


def main():
    terminal_menu = TerminalMenu(
        [
            "Protein Preparation Wizard",
            "Generate MDP Files",
            "Run Molecular Dynamics Simulation",
        ],
        multi_select=True,
        show_multi_select_hint=True,
    )
    terminal_menu.show()
    print("\nYou have chosen:\n")
    for chosen_item in terminal_menu.chosen_menu_entries:  # type: ignore
        print(chosen_item)


if __name__ == "__main__":
    main()
