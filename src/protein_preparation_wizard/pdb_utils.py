from typing import List, Optional
import requests


def download_pdb(pdb_code: str) -> Optional[List[str]]:
    """
    Download a PDB file from the web.

    Args:
        pdb_code (str): The 4-letter PDB ID

    Returns:
        Optional[List[str]]: Contents of the PDB file as a list of strings, otherwise None
    """
    
    pdb_url: str = f"http://www.rcsb.org/pdb/files/{pdb_code}.pdb"

    response = requests.get(pdb_url)

    if response.status_code == 200:
        return requests.get(pdb_url).text.split("\n")
    else:
        print("Error: PDB file not found")
