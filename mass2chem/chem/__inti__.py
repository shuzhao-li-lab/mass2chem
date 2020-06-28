'''
This is a simplified version of pychemy copied from
https://github.com/shuzhao-li/pychemy.git

"simpleChem"?


Installation of openbabel is tricky with python3 binding. Not worthy the trouble?
Openbable is only used in 
inchi.py
from .inchi_converter import convert_inchi_to_formula

Will remove this and visualization elements (matplotlib, networkx) etc., 



To use PubChem API (replacing openbabel?):

https://bioinformatics.stackexchange.com/questions/10755/is-there-a-python-package-to-convert-inchi-to-molecular-structures

import requests
def get_smiles_from_inchikey(inchikey):
    r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON').json()
    return r['PropertyTable']['Properties'][0]['CanonicalSMILES']

inchikey = 'SGNXVBOIDPPRJJ-UHFFFAOYSA-N'
get_smiles_from_inchikey(inchikey)

'''