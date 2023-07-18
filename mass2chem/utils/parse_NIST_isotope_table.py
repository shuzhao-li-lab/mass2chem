import requests as re
import json
import os 

NIST_TABLE_URL = 'https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all'

def fetch_NIST_isotope_table(url=NIST_TABLE_URL):
    r = re.get(url)
    return r.text.split("\n")

def parse_NIST_table_to_isotopes(table_lines, field_names_values=None):
    if field_names_values is None:
        field_names_values = {
            'Atomic Number': int,
            'Atomic Symbol': str,
            'Mass Number': int,
            'Relative Atomic Mass': float,
            'Isotopic Composition': float,
        }
    isotopes = []
    isotope = {}
    for line in table_lines:
        key = line.split(' = ')[0]
        if key in field_names_values:
            try:
                isotope[key] = field_names_values[key](line.rstrip().split("=")[-1].strip().replace('(','').replace(')','').replace('#', ''))
            except:
                isotope[key] = None
        if set(field_names_values.keys()).issubset(set(isotope.keys())):
            isotopes.append(isotope)
            isotope = {}
    return isotopes

def isotopes_to_mass_tables(isotope_list, inject_subatomic=True):
    element_masses = {}
    isotope_masses = {}
    for isotope_dict in isotope_list:
        element = isotope_dict['Atomic Symbol']
        mass_number = isotope_dict['Mass Number']
        abundance = isotope_dict['Isotopic Composition']
        mass = isotope_dict['Relative Atomic Mass']
        if not abundance:
            abundance = 0
        if element == 'D' or element == 'T':
            isotope_masses[element] = mass
            isotope_masses['[' + str(mass_number) + 'H]'] = mass
        else:
            if element not in element_masses:
                element_masses[element] = (mass, abundance)
            else:
                if abundance > element_masses[element][1]:
                    element_masses[element] = (mass, abundance)
            isotope_masses['[' + str(mass_number) + element + ']'] = mass
    subatomic = {
        'e': 0.000549, 
        'PROTON': 1.00727646677
    }
    element_masses = {k: v[0] for k, v in element_masses.items()}
    all_masses = {}
    all_masses.update(element_masses)
    all_masses.update(isotope_masses)
    all_masses.update(subatomic)
    return {
        "elements": element_masses,
        "isotopes": isotope_masses,
        "subatomic": subatomic,
        "all": all_masses 
    }

def dump_mass_isotope_mass_table(mass_tables, path=None):
    if path is None:
        path = os.path.join(os.path.abspath(os.path.dirname(__file__)), '../source_data/NIST_isotope_data.json')
    json.dump(mass_tables, open(path, 'w+'), indent=4)
    
def retrieve_dump_NIST_masses(path=None):
    table_text = fetch_NIST_isotope_table()
    isotopes = parse_NIST_table_to_isotopes(table_text)
    mass_tables = isotopes_to_mass_tables(isotopes)
    dump_mass_isotope_mass_table(mass_tables, path)

if __name__ == '__main__':
    retrieve_dump_NIST_masses()
