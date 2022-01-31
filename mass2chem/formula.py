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

PROTON = 1.00727646677
electron = 0.000549

Ion = namedtuple('Ion', ['mz', 'ion', 'delta_formula'])


def parse_chemformula_dict(x):
    '''This does not deal with nested groups in chemical formula, or isotopes.
    Formula from HMDB 3.5 are compatible.

    Example
    -------
    parse_chemformula_dict('C3H6N2O2')
    {'C': 3, 'H': 6, 'N': 2, 'O': 2}
    '''
    p = re.findall(r'([A-Z][a-z]*)(\d*)', x)
    d = {}
    for pair in p: d[pair[0]] = int( pair[1] or 1 )
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
            (mw - electron, 'M[1+]', {}),                              # e not changing formula
            (mw + PROTON, 'M+H[1+]', {'H': 1}),                        #  
            (mw + 21.9820 + PROTON, 'M+Na[1+]', {'Na': 1}),            # Na = 21.9820 + PROTON = 22.9893
            (mw + 18.0106 + PROTON, 'M+H2O+H[1+]', {'H': 3, 'O':1}),   #  
            ]
        if primary_only:
            return [Ion(*L) for L in primaryList]

        else:
            extendedList = [
                # Use {'C':-1} to force check that C exists in formula
                (mw +1.0034 - electron, 'M(C13)[1+]', {'C':-1, '(C13)':1}), 
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
                (mw + 18.033823, 'M+NH4[1+]', {'H': 4, 'N': 1}),
                
                # Not checking existance of chemical groups here, but should do so when switching to structure based calculations
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
            (mw + electron, 'M[-]', {}), 
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
