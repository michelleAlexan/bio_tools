#%%
from ete3 import NCBITaxa
from pathlib import Path
from typing import Iterable
import warnings

from bio_tools.utils.io import load_yaml

ncbi = NCBITaxa()

UP_TO_DATE_SPECIES_NAME_YAML = Path("/Users/michellealexander/projects/bio_tools/config/up_to_date_species_name.yaml")
#%%
def internet_on(timeout: float = 3.0) -> bool:
    """Check internet by resolving a hostname."""
    import socket
    try:
        socket.setdefaulttimeout(timeout)
        host = socket.gethostbyname("www.google.com")
        s = socket.create_connection((host, 80), timeout)
        s.close()
        return True
    except OSError:
        return False


def ensure_list(x):
    """
    Helper function. If species is a single value and not an iterable, convert to list. 
    """
    if isinstance(x, str):
        return [x]

    if isinstance(x, Iterable):
        return list(x)

    raise TypeError("Species must be a string or an iterable of strings.")

def map_species_to_correct_names(scientific_names: (list[str] )| (str), 
                                 up_to_date_species_name_yaml: Path):
    """
    Take a scientific name(s) that are misspelled 
        and take a path to a yaml file in which wrong names map to up-to-date (correct) names. 
    Return the same list with mapped values. 
    """
    up_to_date_species_name_dict = load_yaml(up_to_date_species_name_yaml)

    scientific_names = ensure_list(scientific_names)
    corrected_list = [up_to_date_species_name_dict[s] if s in up_to_date_species_name_dict.keys() else s for s in scientific_names]
    return corrected_list


def scientific_name_to_tax_id(species:(list[str] )| (str), 
                              update_ncbi_db: bool = False, 
                              up_to_date_species_name_yaml: Path | None = None
                            ) -> (list[int] | int):
    """
    Take a list of species, and return as list with corresponding taxon ids.
    
    For this, convert scientific name to taxonomic id using the ete3.NCBITaxa package. 
        If used for first time, a stable internet is required to fetch the taxonomy database.
            - > It will be saved in under '~/.etetoolkit/taxa.sqlite'
        If you haven't updated it for a while and want to fetch the up-to-date database, 
            - > set 'update_ncbi_db' parameter to True (may take 2 minutes)

    """

    def ete3_ncbi_approach(species):
        """
        Helper function. Fetch the species tax ids using the ete3.NCBITaxa approach. 
        """
        tax_ids = []
        for s in species:
            try:
                taxid = ncbi.get_name_translator([s])[s][0]
                tax_ids.append(taxid)
            except KeyError:
                warnings.warn(f"The species '{s}' was not found by ete3.NCBITaxa.get_name_translator and is thus discarded."
                            "Check for spelling mistakes or if species name is depricated.", UserWarning)
                continue

        return tax_ids
    

    # normalize to list
    species = ensure_list(species)

    # In case the list of scientific species name includes any typos or depricated species names, 
    # it is possible to manually create a yaml file containing the mapping:
    #    depricated_species_name -> up_to_date_species_name.
    if up_to_date_species_name_yaml:
        species = map_species_to_correct_names(species, up_to_date_species_name_yaml)


    # check if ete3 ncbi approach can be used 
    # this fetches the most up to date species - tax id mapping.
    if Path('~/.etetoolkit/taxa.sqlite').expanduser().exists():
        # update taxonomy database 
        if internet_on() and update_ncbi_db:
            ncbi.update_taxonomy_database()
        else:
            tax_ids = ete3_ncbi_approach(species)

    #todo: implement taxoniq offline method

    # return tax ids 
    if len(tax_ids) == 1:
        return tax_ids[0]
    else:
        return tax_ids


# TODO: TEST
def get_taxa_topology(taxa_ids: list[int]):
    """
    
    """
    ncbi = NCBITaxa()

    tree = ncbi.get_topology(taxa_ids, annotate=True)

    for node in tree.traverse():
        # Remove problematic characters from names
        if hasattr(node, "sci_name") and node.sci_name:
            node.name = node.sci_name.replace(":", "_").replace("(", "").replace(")", "").replace(" ", "_")
        elif not node.name:
            node.name = "NA"
    # Write a clean Newick file
    tree.write(outfile="plant_tree_clean.nwk", format=5)




#%%


# %%
