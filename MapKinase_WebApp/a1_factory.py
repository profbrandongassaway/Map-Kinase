from MapKinase_WebApp.a2_kegg_api import KeggAPI
from MapKinase_WebApp.a2_wikipathways_api import WikiPathwaysAPI

def get_pathway_api(source):
    if source == "kegg":
        return KeggAPI()
    elif source == "wikipathways":
        return WikiPathwaysAPI()
    else:
        raise ValueError(f"Unknown pathway source: {source}")
