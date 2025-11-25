import requests, sys
import urllib.parse

def get_url(url, **kwargs):
  """
  Helper function for sending request
  """
  response = requests.get(url, **kwargs)

  if not response.ok:
    print(response.text)
    response.raise_for_status()
    sys.exit()

  return response


# %% -------------------------UNIPROT REQUESTS---------------------------------------

# Documentation: https://rest.uniprot.org/beta/docs/
UNIPROT_API = "https://rest.uniprot.org/beta"

# set of valid UniProt query fields 
VALID_FIELDS = {
    "gene",
    "taxonomy_id",
    "organism_id",
    "protein_name",
    "xref" # accession
}

class UniprotQueryBuilder:
    def __init__(self):
        self.parts = []

    def add_term(self, field, value):
        if field not in VALID_FIELDS:
            raise ValueError(f"Invalid field: {field}. Must be one of {list(VALID_FIELDS)}")
        query_part = f"({field}:{value})"
        self.parts.append(query_part)
        return self  

    def add_raw(self, raw_expression):
        """Allows adding raw custom query parts (must be UniProt-compatible)"""
        self.parts.append(f"({raw_expression})")
        return self

    def and_(self):
        self.parts.append("AND")
        return self

    def or_(self):
        self.parts.append("OR")
        return self

    def not_(self):
        self.parts.append("NOT")
        return self

    def build(self):
        return urllib.parse.quote(" ".join(self.parts))

    def __str__(self):
        return urllib.parse.unquote(self.build())


def advanced_uniprot_search(query, size=25):
    """
    Advance search with offset and size
    """
    url = f"{UNIPROT_API}/uniprotkb/search?query={query}&size={size}"
    r = get_url(url)
    data = r.json()
    
    return r, data

