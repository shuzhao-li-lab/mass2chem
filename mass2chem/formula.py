'''
Handling chemical formulae.

This should switch to a more dedicated chemical tool, 
handling structural formats too.

Primary goal in the Mummichog suite is to compute:
1. associations btw chemical formulae and masses
2. Isotopes
3. In source reactions, including adducts, neutral loss, fragments
4. Biochemical reactions
'''

import re
from collections import namedtuple
import warnings
import os
import json

"""
Isotope information sourced from:

https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses

Coursey, J.S., Schwab, D.J., Tsai, J.J., and Dragoset, R.A. (2015), 
Atomic Weights and Isotopic Compositions (version 4.1). 
[Online] Available: http://physics.nist.gov/Comp [2023, July, 18]. 
National Institute of Standards and Technology, Gaithersburg, MD.

And subatomic masses sourced from:

https://physics.nist.gov/cgi-bin/cuu/Value?mpu or meu

Per NIST website:

NIST web pages are provided as a public service by the National Institute of Standards and Technology (NIST). 
With the exception of material marked as copyrighted, information presented on NIST sites are considered 
public information and may be distributed or copied. Use of appropriate byline/photo/image credits is requested.
"""

Ion = namedtuple('Ion', ['mz', 'ion', 'delta_formula'])
NIST_mass_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'source_data/NIST_isotope_data.json')
print(NIST_mass_path)
atom_mass_dict = json.load(open(NIST_mass_path))['all']
PROTON = atom_mass_dict["PROTON"]
ELECTRON = atom_mass_dict['e']

def parse_chemformula_dict(x):  # return to its original form
    '''This does not deal with nested groups in chemical formula, or isotopes.
    Formula from HMDB 3.5 are compatible.
    Example
    -------
    parse_chemformula_dict('C3H6N2O2')
    {'C': 3, 'H': 6, 'N': 2, 'O': 2}
    Note
    ----
    If formula contains negative number of elements, it will have no warnings but get a errenous dictionary (e.g., 'C2H12O-5' will return {'C': 2, 'H': 12, 'O': 1})
    Be careful of those circumstances.
    '''
    p = re.findall(r'([A-Z][a-z]*)(\d*)', x)
    d = {}
    for pair in p: d[pair[0]] = int( pair[1] or 1 )
    return d



def parse_chemformula_dict_comprehensive(x, remove_additive_cpd = True, handle_negative_number = True):
    '''
    A more inclusive parsing function (compare to `parse_chemformula_dict`) to handle some rare cases of formula parsing
    This also deals with nested groups in chemical formula. And the number of repetitive elements present in the formula will be sumed (e.g., CH2O2H24O4)
    Default: removing additive compounds (e.g., [13C]4H12N2·2HCl will turn into [13C]4H12N2)
    Default: Handle negative number in parsing dictionary (e.g.,, C-2H-2O-2)
    Handles common isotopes.
    Formula from HMDB 3.5 are compatible.

    Example
    -------
    parse_chemformula_dict_comprehensive('C-2H12O5H-4·H2O', remove_additive_cpd=True, handle_negative_number=True)
    {'C': -2, 'H': 8, 'O': 5}
    '''
    if remove_additive_cpd == True: # remove additive compound is true
        if '·' in x:
            x = x.split('·')[0]
    if handle_negative_number == True: # handle formula which has negative value; If the formula somehow has `-` not indicating negative number, this function will not work
        p = re.findall(r'(\[\d*[A-Z][a-z]*\]|[A-Z][a-z]*)(-?\d*)', x)
    else: 
        p = re.findall(r'([A-Z][a-z]*)(\d*)', x)
    k = set([x[0] for x in p]) # handle repetitive element (nested group)
    d = {k:int(0) for k in k}
    for pair in p: d[pair[0]] += int( pair[1] or 1 )
    return d


def add_formula_dict(dict1, dict2):
    '''
    Addition of two formulae as dictionaries.
    This allows calculating formula after a reaction, as dict2 can contain substraction of elements.
    Not as good as using real chemical structures, just a simple approximation.

    '''
    new, result = {}, {}
    for k in set(dict1.keys()).union( set(dict2.keys()) ):
        if k in dict1 and k in dict2:
            new[k] = dict1[k] + dict2[k]
        elif k in dict1:
            new[k] = dict1[k]
        else:
            new[k] = dict2[k]
    # check validity, as no negative number of element is allow in product
    _VALID = True
    for k in new:
        if new[k] < 0:
            _VALID = False
        elif new[k] > 0:     # remove if element is 0
            result[k] = new[k] 
    if _VALID:
        return result
    else:
        return None


def calculate_formula_mass(formula, massDict=atom_mass_dict):
    '''
    Modified from MG's calculate_mass. Rounding is not useful during calculation. 
    '''
    formula_dict = parse_chemformula_dict(formula)
    _m = 0
    for k,v in formula_dict.items():
        _m += v * massDict[k]
    return _m


def calculate_mass(formula_dict, decimal_places=6, massDict=atom_mass_dict):
    '''
    Calculate mass from the formula dictionary, with rounded decimals. 
    Modified for MG. 
    '''
    _m = 0
    for k,v in formula_dict.items():
        _m += v * massDict[k]
    _m = round(_m,decimal_places)
    return _m
    

def check_elemental_subset(Fm1, Fm2):
    '''
    Check if Fm2 a subunit of Fm1; Dictionary representation of chem formula.
    Not used now.
    '''
    Fm2_in_Fm1 = True
    Elements1 = Fm1.keys()
    Elements2 = Fm2.keys()
    if [x for x in Elements2 if x not in Elements1]:    #element not in Fm1
        return False
    else:
        for e in Elements2:
            if Fm2[e] > Fm1[e]: Fm2_in_Fm1 = False
        return Fm2_in_Fm1


def dict_to_hill_formula(d):
    '''
    Generate simple formula from dictionary.
    Hill notion: First C, H and other elements in alphabetical order.
    Carbon isotope (C13) is considered, but not other isotopic elements for now.

    Input 
    -----
    formula dictionary, e.g. {'C': 3, 'H': 6, 'N': 2, 'O': 2} 

    Return
    ------
    Formula in Hill notion, e.g. 'C3H6N2O2'

    '''
    def _str_element_number(int_x):
        # convert int_x to str and '' for int_x is 1. Not checking validity. 
        if int_x == 1:
            return ''
        else:
            return str(int_x)

    elements = list(d.keys())
    elements.sort()
    new_order = [x for x in elements if x not in ['C', 'H', '(C13)']]
    if 'H' in d:
        new_order = ['H'] + new_order
    if '(C13)' in d:
        new_order = ['(C13)'] + new_order
    if 'C' in d:
        new_order = ['C'] + new_order

    formula = ''
    for x in new_order:
        formula += x + _str_element_number(d[x])
    return formula


def __get_adduct_list__(mw, mode, primary_only):
    '''
    Return list of adduct or isotope, [(mz, ion, formula_change_dictionary), ...].
    Next step compute_adducts_formulae will validate the ion.
    Enforcing primary ions because others are not valid without them;
    but logic implemented in search functions, not here.
    '''

    if mode == 'pos': 
        primaryList = [
            (mw - ELECTRON, 'M[1+]', {}),                              # e not changing formula
            (mw + PROTON, 'M+H[1+]', {'H': 1}),                        #  
            (mw + 21.9820 + PROTON, 'M+Na[1+]', {'Na': 1}),            # Na = 21.9820 + PROTON = 22.9893
            (mw + 18.0106 + PROTON, 'M+H2O+H[1+]', {'H': 3, 'O':1}),   #  
            (mw + 18.033823, 'M+NH4[1+]', {'H': 4, 'N': 1}),
            ]
        if primary_only:
            return [Ion(*L) for L in primaryList]

        else:
            extendedList = [
                # Use {'C':-1} to force check that C exists in formula
                (mw +1.0034 - ELECTRON, 'M(C13)[1+]', {'C':-1, '(C13)':1}), 
                (mw +1.0034 + PROTON, 'M(C13)+H[1+]', {'C':-1, '(C13)':1, 'H':1}),
                
                (mw/2 + PROTON, 'M+2H[2+]', {'H': 2}),
                (mw/3 + PROTON, 'M+3H[3+]', {'H': 3}),

                (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]', {'C':-1, '(C13)':1, 'H':2}),
                (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]', {'C':-1, '(C13)':1, 'H':3}),
                # Sulfur:  32S (95.02%), 33S (0.75%), 34S (4.21%), .
                (mw +1.9958 + PROTON, 'M(S34)+H[1+]', {'S':-1, '(S34)':1, 'H':1}),
                # Chlorine has two stable isotopes, 35Cl (75.77%) and 37Cl (24.23%)
                (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]', {'Cl':-1, '(Cl37)':1, 'H':1}),
                
                (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]', {'H':1, 'Na':1}),
                (mw + 37.9555 + PROTON, 'M+K[1+]', {'K':1}),         # K = 37.9555 + PROTON = 38.9628
                (mw + 57.9586 + PROTON, 'M+NaCl[1+]', {'Na': 1, 'Cl': 1}), 
                
                (mw - 18.0106 + PROTON, 'M-H2O+H[1+]', {'H': -1, 'O': -1}), 
                (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]', {'H': -3, 'O': -2}),
                (mw - 17.0265 + PROTON, 'M-NH3+H[1+]', {'H': -2, 'N': -1}),
                (mw - 27.9950 + PROTON, 'M-CO+H[1+]', {'H': 1, 'C': -1, 'O': -1}),
                (mw - 43.9898 + PROTON, 'M-CO2+H[1+]', {'H': 1, 'C': -1, 'O': -2}),
                (mw - 46.0054 + PROTON, 'M-HCOOH+H[1+]', {'H': -1, 'C': -1, 'O': -2}),
                (mw + 67.9874 + PROTON, 'M+HCOONa[1+]', {'H': 1, 'C': 1, 'O': 2, 'Na': 1}),
                (mw - 72.0211 + PROTON, 'M-C3H4O2+H[1+]', {'H': -3, 'C': -3, 'O': -2}),
                (mw + 83.9613 + PROTON, 'M+HCOOK[1+]', {'H': 1, 'C': 1, 'O': 2, 'K': 1}),

                # removed, too infrequent
                # (mw - 67.9874 + PROTON, 'M-HCOONa+H[1+]', 'HCO2Na'),
                # (mw - 83.9613 + PROTON, 'M-HCOOK+H[1+]', 'HCO2K'),
                ]
            return [Ion(*L) for L in primaryList + extendedList]

    elif mode == 'neg':
        primaryList = [
            (mw - PROTON, 'M-H[-]', {'H': -1}),       
            (mw + ELECTRON, 'M[-]', {}), 
            (mw - 18.0106 - PROTON, 'M-H2O-H[-]', {'H': -3, 'O':-1}),   
            (mw + 34.9689, 'M+Cl[-]', {'Cl':1}),
            ]
        if primary_only:
            return [Ion(*L) for L in primaryList]

        else:
            extendedList = [
                (mw + 1.0034 - PROTON, 'M(C13)-H[-]', {'C':-1, '(C13)':1, 'H':-1}),
                (mw + 36.9659, 'M+Cl37[-]', {'(Cl37)':1}),
                (mw + 1.9972 - PROTON, 'M(Cl37)-H[-]', {'Cl':-1, '(Cl37)':1, 'H':-1}),
                (mw + 1.9958 - PROTON, 'M(S34)-H[-]', {'S':-1, '(S34)':1, 'H':-1}),

                (mw/2 - PROTON, 'M-2H[2-]', {'H': -2}),
                (mw + 22.9893 - 2*PROTON, 'M+Na-2H[-]', {'H':-2, 'Na':1}),
                (mw + 38.9628 - 2*PROTON, 'M+K-2H[-]', {'H':-2, 'K':1}),
                
                (mw + 78.9183, 'M+Br[-]', {'Br':1}),
                (mw + 80.9163, 'M+Br81[-]', {'(Br81)':1}),

                (mw + 40.01926853323, 'M+ACN-H[-]', {'C':2, 'H':2, 'N':1}),
                (mw + 1.007825 + 12 + 2*15.99491, 'M+HCOO[-]', {'C':1, 'H':1, 'O':2}),
                (mw + 59.013295, 'M+CH3COO[-]', {'C':2, 'H':3, 'O':2}),
                (mw - PROTON + 15.99491, 'M-H+O[-]', {'H': -1, 'O':1}),
                ]
            return [Ion(*L) for L in primaryList + extendedList]


def compute_adducts_formulae(mw, neutral_formula,  mode='pos', primary_only=False):
    '''
    Calculating isotopes and adducts, return m/z and formulae.

    A dictionary of element changes is used in the calculation.
    Another function add_formula_dict checks validity of required numbers of atoms.
    Electron is not considered in formula calculation, but may affected mass calculation.

    The isotopes, neutral losses and adducts can have multiplications,
    but multiple rounds of search will take care of those.

    Reaction notion in dictionary: 
    {'C': 3, 'H': 6, 'N': 2, 'O': 2} changed by {'C': -1, 'H': 1, 'O': -1} 
        = {'C': 2, 'H': 5, 'N': 2, 'O': 1}

    Return
    -------
    List of adducts: e.g. [(58.53894096677, 'M+2H[2+]', result_formula), ...,]

    '''
    dict_neutral_formula = parse_chemformula_dict(neutral_formula)
    adducts2get = []
    for x in __get_adduct_list__(mw, mode, primary_only):
        # ['mz', 'ion', 'delta_formula']
        result = add_formula_dict(dict_neutral_formula, x.delta_formula)
        if result:
            adducts2get.append( [x.mz, x.ion, dict_to_hill_formula(result)] )

    return adducts2get


def calculate_formula_diff(FM_D1, FM_D2):
    '''
    Calculate the differences between two formula dictionaries (FM_D1 - FM_D2), and return a formula dictionary documenting the differences
    
    '''
    unique_elements = set(list(FM_D1.keys()) + list(FM_D2.keys()))
    
    for e in unique_elements:      # update the dictionary with keys that present in another formula
        if e not in FM_D1:
            FM_D1.update({e:0})
        if e not in FM_D2:
            FM_D2.update({e:0})
    
    diff_dict = {key: FM_D1[key] - FM_D2.get(key, 0) for key in FM_D1}
    diff_dict = {k:v for k,v in diff_dict.items() if v !=0} # remove those with zero
    if not all([x > 0 for x in diff_dict.values()]) or not all([x < 0 for x in diff_dict.values()]):
        warnings.warn('the differences between the two formulas are not all positive or all negative values')
    
    return(diff_dict)

def adjust_salt_formula(F, export_str = True):
    '''
    Formula in the commercially provided products usually contain salts, but ions detected mass spec are not in salt form.
    To convert to compatible format for adduct calculation, salts (e.g., K, Na) need to replace by Hydrogens.
    '''
    exchange_dict = {
    'K': {'H': 1},
    'Na': {'H': 1},
    'Ca': {'H': 2},
    'Mg': {'H': 2}}
    F_D = parse_chemformula_dict(F, remove_additive_cpd = True)
    for k,v in exchange_dict.items():
        if k in F_D:
            for i in range(F_D[k]):
                F_D = add_formula_dict(F_D,v)
            F_D.pop(k,None)
    if export_str:
        res = dict_to_hill_formula(F_D)
    else:
        res = F_D
    return res