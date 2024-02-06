# parse HMDB ver 4, single XML files for all metabolites
# Retrieved from https://hmdb.ca/downloads
# Using Python 3 and lxml
# lxml is not in current Python 3 package. 
# It has similar API as ElementTree but faster. To install:
# pip install lxml

'''
(base) MLG-JGM467:HMDB-5 lish$ time python3 parse_hmdb5.py 
Found 217920 entries.

real	2m37.922s
user	2m1.660s
sys	0m26.192s

# Larger files will have to be partitioned, depending on memeory limit.

In [1]: from parse_hmdb4 import *                                                              

In [2]: myList = parse_a_file(infile, wanted)                                                  
Found 25411 entries.

In [3]: myList[0]                                                                              
Out[3]: 
{'accession': 'HMDB0000001',
 'name': '1-Methylhistidine',
 'chemical_formula': 'C7H11N3O2',
 'monisotopic_moleculate_weight': '',
 'iupac_name': '(2S)-2-amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid',
 'traditional_iupac': '1 methylhistidine',
 'cas_registry_number': '332-80-9',
 'smiles': 'CN1C=NC(C[C@H](N)C(O)=O)=C1',
 'inchi': 'InChI=1S/C7H11N3O2/c1-10-3-5(9-4-10)2-6(8)7(11)12/h3-4,6H,2,8H2,1H3,(H,11,12)/t6-/m0/s1',
 'inchikey': 'BRMWTNUJHUMWMS-LURJTMIESA-N',
 'pathways': '',
 'normal_concentrations': '',
 'abnormal_concentrations': '',
 'diseases': '',
 'drugbank_id': 'DB04151',
 'drugbank_metabolite_id': '',
 'phenol_explorer_compound_id': '',
 'phenol_explorer_metabolite_id': '',
 'foodb_id': 'FDB093588',
 'knapsack_id': '',
 'chemspider_id': '83153',
 'kegg_id': 'C01152',
 'biocyc_id': '',
 'bigg_id': '',
 'wikipidia': '',
 'nugowiki': '',
 'metagene': '',
 'metlin_id': '3741',
 'pubchem_compound_id': '92105',
 'het_id': '',
 'chebi_id': '50599',
 'protein_associations': ''}



In [40]: with open('HMDB0000031.xml', 'w') as O: 
    ...:     O.write(str(ET.tostring(x), 'utf-8')) 
    ...:                                                                                       


>>> prop_extract_dict(m)
{'logp': '0.051', 'logs': '0.67', 'solubility': '484 g/L', 'pka_strongest_acidic': '3.99', 'pka_strongest_basic': '-3.8', 'iupac': '(2S)-2-hydroxybutanoic acid', 'average_mass': '104.105', 'mono_mass': '104.047344118', 'smiles': 'CC[C@H](O)C(O)=O', 'formula': 'C4H8O3', 'inchi': 'InChI=1S/C4H8O3/c1-2-3(5)4(6)7/h3,5H,2H2,1H3,(H,6,7)/t3-/m0/s1', 'inchikey': 'AFENDNXGAFYKQO-VKHMYHEASA-N', 'polar_surface_area': '57.53', 'refractivity': '23.36', 'polarizability': '9.98', 'rotatable_bond_count': '2', 'acceptor_count': '3', 'donor_count': '2', 'physiological_charge': '-1', 'formal_charge': '0', 'number_of_rings': '0', 'bioavailability': 'Yes', 'rule_of_five': 'Yes', 'ghose_filter': 'No', 'veber_rule': 'No', 'mddr_like_rule': 'No'}


>>> deep_extract(m, '{http://www.hmdb.ca}' + 'class')
'Hydroxy acids and derivatives'

'''

from lxml import etree as ET 

infile = 'serum_metabolites.xml'

# fields to retrieve
wanted = ['accession', 'name', 'chemical_formula', 'monisotopic_molecular_weight', 'iupac_name', 
    'traditional_iupac', 'cas_registry_number', 'smiles', 'inchi', 'inchikey',]
     
prop_wanted = [ 'logp', 'number_of_rings',]

taxonomy_wanted = ['kingdom', 'super_class', 'class', 'sub_class']


def extract(obj_metabolite, x):
    try:
        y = obj_metabolite.find(x).text.strip()
        if y:
            return y
        else: return ''
    except AttributeError: 
        return ''

def deep_extract(obj_metabolite, x):
    try:
        y = obj_metabolite.find('.//' + x).text.strip()
        if y:
            return y
        else: return ''
    except AttributeError: 
        return ''


def prop_extract_dict(obj_metabolite, prop='predicted_properties', prefix='{http://www.hmdb.ca}'):
    result = {}
    for x in obj_metabolite.find('.//'+prefix + prop).getchildren():
        result[x.find(prefix+'kind').text] = x.find(prefix+'value').text
    return result



def extract_dict(obj_metabolite, wanted, prefix):
    result = {}
    for x in wanted:
        result[x] = extract(obj_metabolite, prefix+x)
    return result
    

def multi_extract_dict(obj_metabolite, wanted, taxonomy_wanted, prop_wanted, prefix):
	
    result = {}
    
    for x in wanted:
        result[x] = extract(obj_metabolite, prefix+x)
        
    for x in taxonomy_wanted:
	    result[x] = deep_extract(obj_metabolite, prefix+x)
        
    _d = prop_extract_dict(obj_metabolite, prop='predicted_properties', prefix=prefix)
    for x in prop_wanted:
	    result[x] = _d.get(x, '')
    
    return result



def parse_a_file(f, wanted, taxonomy_wanted, prop_wanted, prefix='{http://www.hmdb.ca}'):
    '''
    Input
    =====
    f: HMDB 5 XML file
    wanted: fields to retrieve
    prefix: prefix that is used in the XML file

    Return
    ======
    return as a list of dictionaries. 
    '''
    tree = ET.parse(f)
    root = tree.getroot()
    print("Found %d entries." %len(root))

    results = []
    for child in root.getchildren():
        results.append(
        multi_extract_dict(child, wanted, taxonomy_wanted, prop_wanted, prefix)
        )

    return results


def write_tsv(results, wanted, outfile):
    s = '\t'.join(wanted) + '\n'
    for R in results:
        s += '\t'.join([R[x] for x in wanted]) + '\n'
    with open(outfile, 'w') as O:
        O.write(s)



if __name__ == '__main__':
    # infile = 'serum_metabolites.xml'
    # infile= 'urine_metabolites.xml'
    # infile= 'feces_metabolites.xml'
    infile = 'hmdb_metabolites.xml'

    write_tsv(
        parse_a_file(infile, wanted, taxonomy_wanted, prop_wanted),
        wanted+taxonomy_wanted+prop_wanted, 
        "parsed_"+infile.replace(".xml", ".tsv")
    )
