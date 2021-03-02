'''
Under development.

Ultimately, adducts should be computed per compound based on chemical formula.
The Azimuth DB will host the info and provide compiled releases.

To also consider list in
https://github.com/stanstrup/commonMZ


'''

import re

# Pychemy is incluced as `chem`, a stripped down version.
from .chem.molmass import Formula

PROTON = 1.00727646677

# currency metabolites of ubiquitous presence
currency = ['C00001', 'C00080', 'C00007', 'C00006', 'C00005', 'C00003',
            'C00004', 'C00002', 'C00013', 'C00008', 'C00009', 'C00011', 
            'G11113', '',
            'H2O', 'H+', 'Oxygen', 'NADP+', 'NADPH', 'NAD+', 'NADH', 'ATP', 
            'Pyrophosphate', 'ADP', 'Orthophosphate', 'CO2',]

  
primary_ions = ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]']

wanted_adduct_list = {
    'pos_default': ['M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 
                    'M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]'
                    ],

    'generic_positive': ['M[1+]','M+H[1+]','M+2H[2+]','M+3H[3+]','M(C13)+H[1+]','M(C13)+2H[2+]',
                    'M(C13)+3H[3+]','M(S34)+H[1+]','M(Cl37)+H[1+]','M+Na[1+]','M+H+Na[2+]','M+K[1+]',
                    'M+H2O+H[1+]','M-H2O+H[1+]','M-H4O2+H[1+]','M-NH3+H[1+]','M-CO+H[1+]',
                    'M-CO2+H[1+]','M-HCOOH+H[1+]','M+HCOONa[1+]','M-HCOONa+H[1+]','M+NaCl[1+]',
                    'M-C3H4O2+H[1+]','M+HCOOK[1+]','M-HCOOK+H[1+]',
                    ],
    'negative': ['M-H[-]','M-2H[2-]','M(C13)-H[-]','M(S34)-H[-]','M(Cl37)-H[-]',
                    'M+Na-2H[-]','M+K-2H[-]','M-H2O-H[-]','M+Cl[-]','M+Cl37[-]',
                    'M+Br[-]','M+Br81[-]','M+ACN-H[-]','M+HCOO[-]','M+CH3COO[-]','M-H+O[-]'
                    ],

    # to add options

    }




def formula2mass( x):
    return Formula( x ).isotope.mass



def parse_chemformula(x):
    '''This does not deal with nested groups in chemical formula.
    Formula from HMDB 3.5 are compatible.
    '''
    p = re.findall(r'([A-Z][a-z]*)(\d*)', x)
    d = {}
    for pair in p: d[pair[0]] = int( pair[1] or 1 )
    return d

def check_sub(Fm1, Fm2):
    '''Check if Fm2 a subunit of Fm1; Dictionary representation of chem formula
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


def compute_adducts(mw, cFormula):
    '''
    # new function calculating adducts
    # A few rules here:
    # C13, S34, Cl37 pos only applicable if the atom is in formula
    # and others in `required subgroup`, as the last item in the tuples
    # Better way is to based on chemical structure - future direction.


    Initially, 
     the most frequent derivatives under positive mode are adopted from
        Brown et al. Analyst, 2009, 134, 1322-1332.

    # test
    print(dict_formula['C05587'])
    print(model.dict_cpds_def['C05587'])
    print(model.dict_cpds_mass['C05587'])

    print (compute_adducts( model.dict_cpds_mass['C05587'], dict_formula['C05587'] ))

    '''
    # PROTON = 1.00727646677
    
    addList = [(mw, 'M[1+]', ''), 
         (mw + PROTON, 'M+H[1+]', ''),
         (mw/2 + PROTON, 'M+2H[2+]', ''),
         (mw/3 + PROTON, 'M+3H[3+]', ''),
         (mw +1.0034 + PROTON, 'M(C13)+H[1+]', 'C'),
         (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]', 'C'),
         (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]', 'C'),
         (mw +1.9958 + PROTON, 'M(S34)+H[1+]', 'S'),
         (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]', 'Cl'),
         (mw + 21.9820 + PROTON, 'M+Na[1+]', ''),        # Na = 21.9820 + PROTON = 22.9893
         (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]', ''),
         (mw + 37.9555 + PROTON, 'M+K[1+]', ''),         # K = 37.9555 + PROTON = 38.9628
         (mw + 18.0106 + PROTON, 'M+H2O+H[1+]', ''), 
         (mw - 18.0106 + PROTON, 'M-H2O+H[1+]', 'H2O'), 
         (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]', 'H4O2'),
         (mw - 17.0265 + PROTON, 'M-NH3+H[1+]', 'NH3'),
         (mw - 27.9950 + PROTON, 'M-CO+H[1+]', 'CO'),
         (mw - 43.9898 + PROTON, 'M-CO2+H[1+]', 'CO2'),
         (mw - 46.0054 + PROTON, 'M-HCOOH+H[1+]', 'H2CO2'),
         (mw + 67.9874 + PROTON, 'M+HCOONa[1+]', ''),
         (mw - 67.9874 + PROTON, 'M-HCOONa+H[1+]', 'HCO2Na'),
         (mw + 57.9586 + PROTON, 'M+NaCl[1+]', ''), 
         (mw - 72.0211 + PROTON, 'M-C3H4O2+H[1+]', 'C3H4O2'),
         (mw + 83.9613 + PROTON, 'M+HCOOK[1+]', ''),
         (mw - 83.9613 + PROTON, 'M-HCOOK+H[1+]', 'HCO2K'),
         ] + [
             (mw - PROTON, 'M-H[-]', ''),
        (mw/2 - PROTON, 'M-2H[2-]', ''),
        (mw + 1.0034 - PROTON, 'M(C13)-H[-]', 'C'),
        (mw + 1.9958 - PROTON, 'M(S34)-H[-]', 'S'),
        (mw + 1.9972 - PROTON, 'M(Cl37)-H[-]', 'Cl'),
        (mw + 22.9893 - 2*PROTON, 'M+Na-2H[-]', ''),
        (mw + 38.9628 - 2*PROTON, 'M+K-2H[-]', ''),
        (mw - 18.0106 - PROTON, 'M-H2O-H[-]', 'H2O'),
        (mw + 34.9689, 'M+Cl[-]', ''),
        (mw + 36.9659, 'M+Cl37[-]', ''),
        (mw + 78.9183, 'M+Br[-]', ''),
        (mw + 80.9163, 'M+Br81[-]', ''),
        (mw + 2*12 + 3*1.007825 + 14.00307 - PROTON, 'M+ACN-H[-]', ''),
        (mw + 1.007825 + 12 + 2*15.99491, 'M+HCOO[-]', ''),
        (mw + 3*1.007825 + 2*12 + 2*15.99491, 'M+CH3COO[-]', ''),
        (mw - PROTON + 15.99491, 'M-H+O[-]', ''),
        ]

    dict_cFormula = parse_chemformula(cFormula)
    mydict = {}
    for x in addList:
        if check_sub(dict_cFormula, parse_chemformula(x[2])):
            mydict[x[1]] = x[0]

    return mydict


#
# from mummichog
#



def adduct_function(mw, mode):
    '''
    return a list of derivatives/adducts according to operation mode.
    The most frequent derivatives under positive mode are adopted from
        Brown et al. Analyst, 2009, 134, 1322-1332.
    'dpj_positive' is a customized version for DPJ lab.
    Negative mode is an empirical compilation.
    
    Some derivatives are not possible for some compounds, 
    subject to future upgrade.
    
    
    Paul Benton sent a list used in XCMSonline.
    
    This is to be replaced by pre-computed tables based on chemical formula.
    
    '''
    if mode == 'dpj_positive':
        return [(mw, 'M[1+]'), 
                (mw + PROTON, 'M+H[1+]'),
                (mw/2 + PROTON, 'M+2H[2+]'),
                (mw +1.0034 + PROTON, 'M(C13)+H[1+]'),
                (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]'),
                (mw +1.9958 + PROTON, 'M(S34)+H[1+]'),
                (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]'),
                (mw + 21.9820 + PROTON, 'M+Na[1+]'), 
                (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]'),
                (mw + 37.9555 + PROTON, 'M+K[1+]'), 
                (mw + 67.9874 + PROTON, 'M+HCOONa[1+]'),
                (mw + 83.9613 + PROTON, 'M+HCOOK[1+]'),
                ]
    
    elif mode == 'generic_positive':
        return [(mw, 'M[1+]'), 
                (mw + PROTON, 'M+H[1+]'),
                (mw/2 + PROTON, 'M+2H[2+]'),
                (mw/3 + PROTON, 'M+3H[3+]'),
                (mw +1.0034 + PROTON, 'M(C13)+H[1+]'),
                (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]'),
                (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]'),
                (mw +1.9958 + PROTON, 'M(S34)+H[1+]'),
                (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]'),
                #
                (mw + 21.9820 + PROTON, 'M+Na[1+]'),    # Na = 21.9820 + PROTON = 22.9893
                (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]'),
                (mw + 37.9555 + PROTON, 'M+K[1+]'),     # K = 37.9555 + PROTON = 38.9628
                #
                (mw + 18.0106 + PROTON, 'M+H2O+H[1+]'), 
                (mw - 18.0106 + PROTON, 'M-H2O+H[1+]'), 
                (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]'),
                (mw - 17.0265 + PROTON, 'M-NH3+H[1+]'),
                (mw - 27.9950 + PROTON, 'M-CO+H[1+]'),
                (mw - 43.9898 + PROTON, 'M-CO2+H[1+]'),
                (mw - 46.0054 + PROTON, 'M-HCOOH+H[1+]'),
                (mw + 67.9874 + PROTON, 'M+HCOONa[1+]'),
                (mw - 67.9874 + PROTON, 'M-HCOONa+H[1+]'),
                #
                (mw + 57.9586 + PROTON, 'M+NaCl[1+]'), 
                (mw - 72.0211 + PROTON, 'M-C3H4O2+H[1+]'),
                (mw + 83.9613 + PROTON, 'M+HCOOK[1+]'),
                (mw - 83.9613 + PROTON, 'M-HCOOK+H[1+]'),
                ]
    
    elif mode == 'park':
        return [(mw + PROTON, 'M+H[1+]'),
                (mw + 21.9820 + PROTON, 'M+Na[1+]'), 
                #(mw + 18.0106 + PROTON, 'M+H2O+H[1+]'), 
                (mw - 18.0106 + PROTON, 'M-H2O+H[1+]'), 
                (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]'),
                ]
        
    elif mode == 'negative':
        return [(mw - PROTON, 'M-H[-]'),
               (mw/2 - PROTON, 'M-2H[2-]'),
               (mw + 1.0034 - PROTON, 'M(C13)-H[-]'),
               (mw + 1.9958 - PROTON, 'M(S34)-H[-]'),
               (mw + 1.9972 - PROTON, 'M(Cl37)-H[-]'),
               #
               (mw + 22.9893 - 2*PROTON, 'M+Na-2H[-]'),
               (mw + 38.9628 - 2*PROTON, 'M+K-2H[-]'),
               (mw - 18.0106 - PROTON, 'M-H2O-H[-]'),
               #
               (mw + 34.9689, 'M+Cl[-]'),
               (mw + 36.9659, 'M+Cl37[-]'),
               (mw + 78.9183, 'M+Br[-]'),
               (mw + 80.9163, 'M+Br81[-]'),
               (mw + 2*12 + 3*1.007825 + 14.00307 - PROTON, 'M+ACN-H[-]'),
               (mw + 1.007825 + 12 + 2*15.99491, 'M+HCOO[-]'),
               (mw + 3*1.007825 + 2*12 + 2*15.99491, 'M+CH3COO[-]'),
               (mw - PROTON + 15.99491, 'M-H+O[-]'),
               ]

    elif mode == 'neutral':
        print ("Neutral mode of instrumentation is not supported.")
        return []
    
    else:
        print ("Unrecognized mode of instrumentation.")
        return []

# not used
# accuracy of the MS instrument
def mz_tolerance(mz, instrument):
    '''
    instrument is ppm in the input parameter. 
    Using flat ppm now.
    Future experimental function to add instrument specific calibration, but no point right now.
    '''
    try:
        instrument = int(instrument)
        return 0.000001 * instrument * mz
    
    except ValueError:
        if instrument == 'FTMS':
            return max(0.00001*mz, 3*2**(1.98816*np.log2(mz) - 26.1945))
        
        elif instrument == 'ORBITRAP':
            # needs further calibration
            # return max(0.00001*mz, 2*2**(1.45325*np.log2(mz) - 20.8554))
            return 0.000010*mz
        
        else:
            return 0.000010*mz
  


