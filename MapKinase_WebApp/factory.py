from kegg_api import KeggAPI
from wikipathways_api import WikiPathwaysAPI
from pathbank_api import PathBankAPI

def get_pathway_api(source):
    if source == "kegg":
        return KeggAPI()
    elif source == "wikipathways":
        return WikiPathwaysAPI()
    elif source == "pathbank":
        return PathBankAPI()
    else:
        raise ValueError(f"Unknown pathway source: {source}")