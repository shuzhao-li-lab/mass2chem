'''
DB based, step wise, formula mass centered MS1 annotation

use case 1:
    input feature table
    output:
        feature table with flags of contaminants, empCpd and corresponding ions
        nonredundant table with one feature (selected by highest intensity) per empCpd, and empCpd annotations

use case 2:
    single sample annotation, to faciliate correspondence in data preprocessing


To-do:
# to import library of authentic chemical standards

# core data structures
# from metDataModel import derived

Old data structure in mummichog v2:
"C05769": {"formula": "C36H34N4O8", "mw": 654.269, "name": "Coproporphyrin I", 
"adducts": {"M+2H[2+]": 328.14177646677, "M+Br81[-]": 735.1853, "M-H2O+H[1+]": 637.2656764667701, "M-C3H4O2+H[1+]": 583.25517646677, "M-HCOOH+H[1+]": 609.27087646677, "M-CO+H[1+]": 627.28127646677, 
        "M+K[1+]": 693.23177646677, "M+Cl[-]": 689.2379, "M+Na-2H[-]": 674.23644706646, "M-CO2+H[1+]": 611.28647646677, "M+Na[1+]": 677.25827646677, "M-2H[2-]": 326.12722353323, "M+H[1+]": 655.27627646677, 
        "M-H4O2+H[1+]": 619.25507646677, "M(C13)-H[-]": 654.26512353323, "M+HCOONa[1+]": 723.26367646677, "M(C13)+2H[2+]": 328.64347646677004, "M+HCOOK[1+]": 739.23757646677, "M+HCOO[-]": 699.266645, 
        "M(C13)+3H[3+]": 219.43134313343666, "M-H[-]": 653.26172353323, "M+ACN-H[-]": 694.2882685332299, "M+Cl37[-]": 691.2349, "M-H2O-H[-]": 635.25112353323, "M+Br[-]": 733.1873, 
        "M+3H[3+]": 219.09694313343667, "M+CH3COO[-]": 713.282295, "M(C13)+H[1+]": 656.2796764667701, "M[1+]": 654.269, "M-NH3+H[1+]": 638.2497764667701, "M+NaCl[1+]": 713.2348764667701, 
        "M+H+Na[2+]": 339.13277646677, "M+H2O+H[1+]": 673.28687646677, "M-H+O[-]": 669.25663353323, "M+K-2H[-]": 690.20994706646}}, 

Use DB1, DB2, (DB3, ...) for tiered searches.
Start with DB1 for searching primary ions, because other adducts/isotopes can only exist if a primary ion does.
Precompute selectivity within each DB. Use in indexed format.

Use indexed_DB            
example DB_1[99] = 
            [['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]], 
            ['C4H5NS_99.013721', 99.01372099999999, 'C4H5NS', 0.9967247109333564, [('C4H5NS_99.01427', 'M[1+]')]], 
            ['C3H8NaO2_99.041706', 99.04170646677, 'C3H8NaO2', 0.999999999998776, [('C3H8O2_76.05243', 'M+Na[1+]')]], ...]

empCpds will be organized based on matched DB entries after formula_mass annotation.


'''

# import re
import numpy as np
# from scipy.optimize import curve_fit
from scipy.stats import norm as normal_distribution

# from metDataModel.core import import Feature, EmpiricalCompound

# databases
from .lib.common_mass import mass_signatures
from .lib.massDicts import massDict_hmdb_serum, massDict_hmdb
from .lib.PubChemLite import PubChemLite

from .lib.LCMS_contaminants import contaminants_pos, contaminants_neg

from .io import tFeature, read_feature_tuples
from .formula import compute_adducts_formulae

# from .calibrate import mass_calibrate


# ad hoc, will repalce with better lists for calibration, for each of pos and neg.
calibration_mass_dict_pos = {'C4H9N3O2_131.069477': 132.07675346677001, 'C6H11NO2_129.078979': 130.08625546677, 'C5H9NO4_147.053158': 148.06043446677, 'C6H10O6_178.047738': 179.05501446677002, 'C9H11NO3_181.073893': 182.08116946677, 'C9H11NO2_165.078979': 166.08625546677, 'C9H8O3_164.047344': 165.05462046677002, 'C23H45NO4_399.334859': 400.34213546677, 'C2H7NO3S_125.014664': 126.02194046676999, 'C5H7NO3_129.042593': 130.04986946677002, 'C5H4N4O3_168.02834': 169.03561646677, 'C5H8O5_148.037173': 149.04444946677, 'C5H12O5_152.068473': 153.07574946677002, 'C9H8O2_148.052429': 149.05970546677, 'C5H10N2O3_146.069142': 147.07641846677, 'C17H33NO4_315.240959': 316.24823546677, 'C5H11NO2S_149.051049': 150.05832546677001, 'C5H10O2_102.06808': 103.07535646676999, 'C15H29NO4_287.209658': 288.21693446677, 'C8H17NO2_159.125929': 160.13320546677002, 'H2O4S_97.967379': 98.97465546676999, 'C6H13NO5_179.079373': 180.08664946677, 'C16H22O4_278.151809': 279.15908546677, 'C20H28O2_300.20893': 301.21620646677, 'C7H8N4O2_180.064726': 181.07200246677002, 'C7H6O2_122.036779': 123.04405546676999, 'C14H22N2O3_266.163043': 267.17031946677, 'C19H16O4_308.104859': 309.11213546677, 'C21H39NO4_369.287909': 370.29518546677, 'C8H6O4_166.026609': 167.03388546677002, 'C18H35NO_281.271865': 282.27914146677, 'C19H37NO4_343.272259': 344.27953546677, 'C26H52NO7P_521.34814': 522.35541646677, 'C7H8N2O2_152.058578': 153.06585446677002, 'C9H17NO2_171.125929': 172.13320546677002, 'C25H47NO4_425.350509': 426.35778546677, 'C8H10O3_154.062994': 155.07027046677, 'C25H45NO4_423.334859': 424.34213546677, 'C22H46NO7P_467.301189': 468.30846546677003, 'C23H48NO7P_481.316839': 482.32411546677, 'C24H48NO7P_493.316839': 494.32411546677, 'C26H54NO7P_523.36379': 524.37106646677, 'C26H48NO7P_517.316839': 518.32411546677, 'C28H52NO7P_545.34814': 546.35541646677, 'C28H50NO7P_543.332489': 544.33976546677, 'C24H50NO6P_479.337575': 480.34485146677, 'C25H52NO7P_509.34814': 510.35541646677, 'C21H38O4_354.27701': 355.28428646677, 'C17H31NO4_313.225308': 314.23258446677, 'C15H27NO4_285.194008': 286.20128446677, 'C19H35NO4_341.256609': 342.26388546677003, 'C22H37NO3_363.277344': 364.28462046677004, 'C19H22ClN5O_371.151288': 372.15856446677003, 'C24H30O6_414.204239': 415.21151546677, 'C21H25N_291.1987': 292.20597646677, 'C15H23NO2_249.172879': 250.18015546677, 
    'C17H27NO4_309.194008': 310.20128446677, 'C21H25ClN2O3_388.15537': 389.16264646677, 'C16H18O10_370.089997': 371.09727346677}


MASS_RANGE = (50, 2000)
# fraction of total retention time, or of ranks of retention time
# used to determine coelution of ions ad hoc
RETENTION_TIME_TOLERANCE_FRAC = 0.02

FEATURE_REGISTRY = {}
EMPCPD_REGISTRY = {}


def make_formula_mass_id(formula, mass):
    return formula + '_' + str(round(mass,6))

def list_search_formula_mass_db(list_query_mz, indexed_DB, limit_ppm=10):
    return [search_formula_mass_db(query_mz, indexed_DB, limit_ppm) for query_mz in list_query_mz]

def bin_by_median(List_of_tuples, func_tolerance):
    '''
    Not perfect because left side may deviate out of tolerance, but LC-MS data always have enough gaps for separation.
    Will add kernel density method for grouping m/z features.
    List_of_tuples: [(value, object), (value, object), ...], to be separated into bins by values (either rt or mz).
                    objects have attribute of sample_name if to align elsewhere.
    return: [seprated bins], each as a list in the same input format. If all falls in same bin, i.e. [List_of_tuples].
    '''
    new = [[List_of_tuples[0], ], ]
    for X in List_of_tuples[1:]:
        if X[0]-np.median([ii[0] for ii in new[-1]]) < func_tolerance(X[0]):       # median moving with list change
            new[-1].append(X)
        else:
            new.append([X])
    PL = []
    for L in new:
        PL.append([X[1] for X in L])
    return PL


def assemble_empcpds_by_formula(List_of_features):
    pass


def assemble_empcpds_denovo(List_of_features):
    pass


def search_formula_mass_db(query_mz, indexed_DB, limit_ppm=10):
    '''
    Find best matched formula_mass in indexed_DB for query_mz within ppm limit.
    
    indexed_DB: # [formula_mass, m/z, charged_formula, selectivity, [relations to neutral formula]]
            example DB_1[99] = 
            [['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]], 
            ['C4H5NS_99.013721', 99.01372099999999, 'C4H5NS', 0.9967247109333564, [('C4H5NS_99.01427', 'M[1+]')]], 
            ['C3H8NaO2_99.041706', 99.04170646677, 'C3H8NaO2', 0.999999999998776, [('C3H8O2_76.05243', 'M+Na[1+]')]], ...]

    return
    ------
    Closest match as an entry,
             (match_ppm, [formula_mass, m/z, charged_formula, selectivity, [relations to neutral formula]])
    None if out of ppm limit.
    '''
    _delta = query_mz * limit_ppm * 0.000001
    _low_lim, _high_lim = query_mz - _delta, query_mz + _delta
    result = []
    for ii in range(int(_low_lim), int(_high_lim)+1):
        if ii in indexed_DB:
            for F in indexed_DB[ii]:
                # F[1] is m/z, fixed DB format
                result.append( (abs(query_mz-F[1]), F) )

    result.sort()
    if result and result[0][0] < _delta:
        return result[0][1]
    else:
        return None

def create_extended_DB(list_empCpd_seeds, mode='pos'):
    '''
    This creates a new DB of common isotopes and adducts based on a given seed list.
    This is using ions specified in .formula.compute_adducts_formulae.
    This can be used on what's found after getting primary ions in DB1, 
    because isotopes/adducts depend on the existence of a primary ion.
    Otherwise, the presence of a primary ion will be enforced at empCpd level.

    input
    -----
    list_empCpd_seeds: list of entries as [formula_mass, ...]. Not full empCpd class.
    
    return
    ------
    indexed_DB in the format of each entry like:
    [formula_mass, m/z, charged_formula, selectivity, [relations to neutral formula]]

    '''
    new_DB = {}
    for x in set(list_empCpd_seeds):             # remove redundancy
        neutral_formula, mw = x.split("_")       # [('C2H4O3_76.016044',..
        adducts = compute_adducts_formulae(float(mw), neutral_formula, mode)    
                                                 # e.g. [(58.53894096677, 'M+2H[2+]', result_formula), ...,]
        for A in adducts:
            new = [ make_formula_mass_id(A[2], A[0]), A[0], A[2], None, [(x, A[1]),] ]
            ii = int(A[0])
            if ii in new_DB:
                new_DB[ii].append(new)
            else:
                new_DB[ii] = [new]

    return new_DB


def annotate_formula_mass(list_query_mz, indexed_DB, check_mass_accuracy=False, limit_ppm=10):
    '''
    Annotate input m/z list using formula based mass in the indexed_DB,
    which should be generic, not assuming separation of primary ions and their common isotopes and adducts.
    Search efficiency will be improved by the construction strategy of the indexed_DB,
    not dependent on the algorithm in this function. E.g., more frequent features should be in first DB to search.

    If check_mass_accuracy is True, the ppm deviations of each matched primary ion are fitted to a normal distribution,
    to estimate mass_accuracy and ppm_std.

    input
    -----
    input_list: m/z values.
    indexed_DB: per entry - [formula_mass, m/z, charged_formula, ion, selectivity, [relations to neutral formula]]
                This can be a pre-loaded DB, or hot DB in asari.

    # mode='pos', 

    return
    ------
    {mass_accuracy: estimated based on difference from formula value
    ppm_std:         estimated based on precision distribution
    annotated_list: [matched DB entry or None, ...]
    corrected_list_query_mz: m/z values will be corrected if they deviate from formula values by mass_accuracy > 5 ppm.
    }
    '''
    result = list_search_formula_mass_db(list_query_mz, indexed_DB, limit_ppm)
    mass_accuracy, ppm_std, corrected_list_query_mz = None, None, None
    N = len(list_query_mz)
    if check_mass_accuracy:
        list_ppm_errors = []
        for ii in range(N):
            if result[ii]:
                list_ppm_errors.append( 1000000 * (list_query_mz[ii] - result[ii][1])/result[ii][1] )

        mass_accuracy, ppm_std = normal_distribution.fit(list_ppm_errors)
        if abs(mass_accuracy) > 5:   # this is considered significant mass shift, requiring m/z correction for all 
            list_query_mz = [x-x*0.000001*mass_accuracy for x in list_query_mz]
            # redo search because we may find different matches after mass correction
            # Also change search ppm if self.__mass_stdev__ > 5
            result = list_search_formula_mass_db(list_query_mz, indexed_DB, max(10, 2*ppm_std))
            corrected_list_query_mz = list_query_mz
            # update ppm_std because it may be used for downstream search parameters
            list_ppm_errors = []
            for ii in range(N):
                if result[ii]:
                    list_ppm_errors.append( 1000000 * (list_query_mz[ii] - result[ii][1])/result[ii][1] )
            _mu, ppm_std= normal_distribution.fit(list_ppm_errors)

        elif ppm_std > 5:      # no mass correction but need redo search using higher ppm for unmatched mz
            result = list_search_formula_mass_db(list_query_mz, indexed_DB, 2*ppm_std)

    return {
        'mass_accuracy': mass_accuracy,
        'ppm_std': ppm_std,
        'annotated_list': result,
        'corrected_list_query_mz': corrected_list_query_mz,
        }


# -----------------------------------------------------------------------------
#
# Prototype code with short-handed data structure, 
# still used for particular cases, backward compatibility only during dev
#
# -----------------------------------------------------------------------------

def do_two_step_annotation(indexed_features, anchorDict, mode='pos',  
                                            MASS_RANGE=MASS_RANGE, ppm=2):
    '''
    Wrapper of annotate_anchorList and extend_annotate_anchorList.
    '''
    _empCpds, _unmatched_features = annotate_anchorList(
                                            indexed_features, anchorDict, ppm)
    _empCpds, _unmatched_features = extend_annotate_anchorList( 
                                            _empCpds, _unmatched_features, mode,  MASS_RANGE, ppm)
    print("Got matched empCpds: ", len(_empCpds))

    return _empCpds, _unmatched_features


def annotate_anchorList(indexed_features, anchorList, ppm=2):
    '''
    Annotate indexed_features by an anchor set of formula_mass values.
    
    Input
    -----
    indexed_features: {int_mz: [{'id': '', 'mz': 0, 'calibrated_mz': 0, 'rtime': 0, ...}, ...],
                                ...}
    anchorList: {'C7H11N3O2_169.085127': [(170.09240346677, '+H[1+]'), (192.074345, 'Na'), (169.085127, 'M')], ...}
    ppm: Based on m/z only, but for calibrated m/z, this should be very stringent, thus ppm default at 2.

    Return
    -------
    matched_empCpds: List of empCpds, e.g. below.
    unmatched_features: Features not matched to anchorList, unindexed
    
    Example
    -------
    matched_empCpds, unmatched_features = annotate_anchorList(my_cal_flaged, anchorList, 2)

    726 702 3939

    matched_empCpds[33]

    {'list_of_features': [{'id': '"FT0394"',
    'mz': 116.070646922025,
    'rtime': 36.3316230773926,
    'intensities': [],
    'calibrated_mz': 116.07051083925319,
    'potential_contaminants': [],
    'ion': '+H[1+]'},
    {'id': '"FT0670"',
    'mz': 138.052651169964,
    'rtime': 36.1285610198975,
    'intensities': [],
    'calibrated_mz': 138.052503250255,
    'potential_contaminants': [],
    'ion': 'Na'},
    {'id': '"FT0379"',
    'mz': 115.063235073422,
    'rtime': 31.9573593139648,
    'intensities': [],
    'calibrated_mz': 115.06309953312442,
    'potential_contaminants': [],
    'ion': 'M'}],
    'interim_id': 'C5H9NO2_115.063329',
    'neutral_formula': 'C5H9NO2',
    'neutral_base_mass': 115.063329}

    len(matched_empCpds), len(unmatched_features)

    (504, 3939)

    # this example shows 504 empCpds are found in 4641 features.
    # The empCpds include 702 Features, of which 24 have multiple matches.
    '''

    matched_empCpds, unmatched_features, matched_features = [], [], []
    for anchor_compound in list(anchorList.items()):
        _m = __anchor_to_empCpd__(indexed_features, anchor_compound, ppm)
        if _m:
            matched_empCpds.append(_m)

    for F in matched_empCpds: 
        matched_features += [x['id'] for x in F['list_of_features']]

    for k,v in indexed_features.items():
        for F in v:
            if F['id'] not in matched_features:
                unmatched_features.append(F)

    print("Remaining unmatched featuress: ", len(unmatched_features))
    return matched_empCpds, unmatched_features


def extend_annotate_anchorList(list_empCpds, unmatched_features, mode='pos',  MASS_RANGE=MASS_RANGE, ppm=2):
    '''
    Annotate unmatched_features by extended adducts.
    
    Input
    -----
    list_empCpds: List of anchor matched empCpds
    unmatched_features: [{'id': '', 'mz': 0, 'calibrated_mz': 0, 'rtime': 0, ...}, ...],
    mode: ionization mode, pos or neg.
    ppm: Based on m/z only, but for calibrated m/z, this should be very stringent, thus ppm default at 2.

    Return
    -------
    matched_empCpds: List of updated empCpds, 
    updated_unmatched_features: Features not in matched_empCpds, unindexed
    '''

    unmatched_features = index_features(unmatched_features)
    return __extend_search__(list_empCpds, unmatched_features, mode,  MASS_RANGE, ppm)





def flag_contaminants(indexed_features, contaminant_list, ppm=2):
    '''
    Flag features if they match to known contaminants;
    set F['potential_contaminants'] = [contaminants, ].
    Based on m/z only, but for calibrated m/z, this should be very stringent, thus ppm default at 2.

    contaminant_list is dependent on ion mode, pos or neg. There's some redundancy in the lists.
    # mass, good_name, name, formula, ion_form, possible origin
    In [14]: contaminants_pos[:2]                                                                                             
    Out[14]: 
    [['537.8790134',
    'C2H4O2_[M6-H6+Fe3+O]+_537.879013',
    'Acetic Acid',
    'C2H4O2',
    '[M6-H6+Fe3+O]+',
    'Solvent'],
    ['555.8895784',
    'C2H4O2_[M6-H6+H2O+Fe3+O]+_555.889578',
    'Acetic Acid',
    'C2H4O2',
    '[M6-H6+H2O+Fe3+O]+',
    'Solvent']]

    This can be faster, as now looping unindexed contaminant_list

    Input
    -----
    indexed_features: {int_mz: [{'id': '', 'mz': 0, 'calibrated_mz': 0, 'rtime': 0, ...}, ...],
                                ...}

    '''
    def check_contaminant(mz, contaminant_list, ppm):
        matched = []
        for T in contaminant_list:
            if abs(T[0] - mz)/mz < 0.000001*ppm:
                matched.append(T)
        return matched

    for k,v in indexed_features.items():
        for F in v:
            F['potential_contaminants'] = check_contaminant(F['calibrated_mz'], contaminant_list, ppm)

    return indexed_features


def make_anchordict(massDict_db, mode='pos', MASS_RANGE=MASS_RANGE):
    '''
    To build anchor mz list based on keys like 'C7H11N3O2_169.085127'.
    This core list is used to separate ions from a feature list to reduce search space.
    The matched results will be organized by the same key of formula_mass.

    The core list is designed for mandating the presence of one of the primary ions, 
    without which other adducts won't be considered.

    Input
    -----
    massDict_db: using formula_mass as identifiers, e.g. 
                'C7H11N3O2_169.085127':[['HMDB0000001', '1-Methylhistidine'], ['HMDB0000479', '3-Methylhistidine']]

    Return
    ------
    {'C7H11N3O2_169.085127': [(170.09240346677, '+H[1+]'), (192.074345, 'Na'), (169.085127, 'M')], ...}
    '''
    posList = [(1.00727646677, 'M+H[1+]'), (22.989218, 'Na'), (-0.0005, 'M*'), 
                (19.01787646677, 'M+H2O+H[1+]')]
    # ? 'M-2H[2-]'
    negList = [(-1.00727646677, 'M-H[-]'), (34.969402, 'M+Cl-'), (0.0005, 'M+e[-]'), 
                (-19.01787646677, 'M-H2O-H[-]')]

    if mode == 'pos': useList = posList
    elif mode == 'neg': useList = negList
    else: print ("ionization mode is either pos or neg.")

    new = {}
    for k in massDict_db:
        kmass = float(k.split('_')[1])
        new[k] = [(x+kmass, y) for x,y in useList if MASS_RANGE[0] < x+kmass < MASS_RANGE[1]]
    return new

def search_feature_pair(F1, F2, mass_signatures, ppm=2):
    '''
    # unidirectional, F2-F1 only
    # msg is like (1.00727646677, '+H[1+]'),

    return
    ------
    list of matched: [('"FT1104"', '"FT1630"', 43.96389, '2Na-2H'), ...]
    '''
    _matched = []
    for msg in mass_signatures:
        if abs( F2['calibrated_mz'] - F1['calibrated_mz'] - msg[0] ) / F2['calibrated_mz'] < 0.000001*ppm:
            _matched.append( (F1['id'], F2['id'], msg[0], msg[1]) )
            
    return _matched


# -----------------------------------------------------------------------------
#
# Functions below are for internal use. Don't touch.
#
# -----------------------------------------------------------------------------


def __find_mz_indexed_features__(query_mz, indexed_features, limit_ppm=2):
    # to find all matches of a m/z in indexed_features, under limit_ppm
    matched = []
    _delta = query_mz * limit_ppm * 0.000001
    low_lim, high_lim = query_mz - _delta, query_mz + _delta
    for ii in range(int(low_lim), int(high_lim)+1):
        for F in indexed_features[ii]:
            if low_lim < F['calibrated_mz'] < high_lim:
                matched.append(F)

    return matched

def __anchor_to_empCpd__(indexed_features, anchor_compound, ppm):
    '''
    Annotate indexed_features into empCpd.
    EmpiricalCompound = {
        "neutral_base_mass": 0,
        "neutral_formula": [],
        "interim_id": '', 
        "list_of_features": [
                {'id': '', 'ion': 'M-H[1-]', 'm/z': 169.0013, 'rtime': 55, ..},
                {},
    }
    
    Input
    -----
    indexed_features: {int_mz: [{'id': '', 'mz': 0, 'calibrated_mz': 0, 'rtime': 0, ...}, ...],
                                ...}
    anchor_compound: ('C7H11N3O2_169.085127', [(170.09240346677, '+H[1+]'), (192.074345, 'Na'), (169.085127, 'M')])
    ppm: Based on m/z only, but for calibrated m/z, this should be very stringent, thus ppm default at 2.

    Return
    -------
    empCpd if matched; otherwise {}
    ''' 
    _matched = False
    tmp = {
        "list_of_features": [],
    }
    for v in anchor_compound[1]:
        found_features = __find_mz_indexed_features__(v[0], indexed_features, ppm)
        if found_features:
            _matched= True
            for F in found_features: F['ion'] = v[1]
            tmp['list_of_features'] += found_features

    if _matched:
        tmp['interim_id'] = anchor_compound[0]
        neutral_formula, neutral_base_mass = anchor_compound[0].split('_')
        tmp['neutral_formula'], tmp['neutral_base_mass'] = neutral_formula, float(neutral_base_mass)
        return tmp
    else:
        return None



def __extend_search__(list_empCpds, indexed_unmatched_features, 
                            mode='pos',  MASS_RANGE=MASS_RANGE, ppm=2):
    '''
    Annotate indexed_features by extended adducts.
    
    Input
    -----
    list_empCpds: List of anchor matched empCpds
    indexed_unmatched_features: indexed_features as {int_mz: [{'id': '', 'mz': 0, 'calibrated_mz': 0, 'rtime': 0, ...}, ...],
                                ...}
    mode: ionization mode, pos or neg.
    ppm: Based on m/z only, but for calibrated m/z, this should be very stringent, thus ppm default at 2.

    Return
    -------
    matched_empCpds: List of updated empCpds, 
    updated_unmatched_features: Features not in matched_empCpds, unindexed

    '''
    _tmp_matched = []
    for E in list_empCpds:
        for target_mz, note, result_formula in compute_adducts_formulae(E['neutral_base_mass'], E['neutral_formula'], mode) :
            if MASS_RANGE[0] < target_mz < MASS_RANGE[1]:
                found_features = __find_mz_indexed_features__(target_mz, indexed_unmatched_features, ppm)
                for F in found_features: 
                    _tmp_matched.append(F['id'])
                    # format?
                    F['ion'] = (result_formula, note)
                    E['list_of_features'].append(F)

    updated_unmatched_features = []
    for k,v in indexed_unmatched_features.items():
        for F in v:
            if F['id'] not in _tmp_matched:
                updated_unmatched_features.append(F)

    print("Remaining unmatched featuress: ", len(updated_unmatched_features))
    return list_empCpds, updated_unmatched_features







#
# -----------------------------------------------------------------------------
#
# not used now
# 

def find_isotopic_pairs(list_peaks, mztree, isotopic_pair, mz_tolerance_ppm=5, rt_tolerance_scans=5):
    '''
    - to update

    To find naturally occuring isotopic pairs in LC-MS peaks, 
    using strict coelution criteria and loose abundance ratio.
    This is done within a sample, so that all elution time is precisely matched by scan numbers.
    This avoids any issues in alignment across samples.
    Rules: elution apex should be within 5 scans; overlap of elution peaks minimal 50% of the smaller peak;
            abundance ratio loosely applied, leaving room for measuring errors.
    
    Abundance of natually occuring isotopic elements:
    12C~99%, 13C ~ 1%
    14N ~ 99.64%, 15N ~ 0.36%
    Sulfur:  32S (95.02%), 33S (0.75%), 34S (4.21%), 
    Chlorine has two stable isotopes, 35Cl (75.77%) and 37Cl (24.23%).
    1H and 31P are basically ~100% in nature.

    Input
    =====
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
    'left_base': 648, 'right_base': 655, 'id_number': 555},
                 ...]
    mztree is constructed on list_peaks.
    isotopic_pair: (mass_difference, notion, abundance_ratio). E.g. (1.0034, 'M(C13)', 0.2). - to update
    
    Return
    ======
    list of pairs of peaks.
    '''
    pairs = []
    (mass_difference, relation, abundance_ratio_min, abundance_ratio_max) = isotopic_pair
    for P1 in list_peaks:
        tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
        for P2 in tmp:
            if abs(P1['apex']-P2['apex']) <= rt_tolerance_scans and is_coeluted(P1, P2) \
                                             and abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                pairs.append((P1, P2))

    return pairs


def init_eTrees_by_isotopic_signatures(isotopic_signatures, iso_type, peak_list):
    '''
    An example of isotopic_signatures: 
    [(296, 'anchor'), (339, 'M(13C)'), (332, 'M(15N)'), (333, 'M(15N)')], using peak numbers. 
    

    '''
    eTrees = []
    (P1, P2) = iso_pairs
    for pp in iso_pairs:
        ET = epdTree(peak_list)
        ET.anchor_peak = P1['id_number']
        ET.isotope_peaks = [P1['id_number'], P2['id_number']]
        ET.peak_relationships.append(( P1['id_number'], P2['id_number'], iso_type ))
    return eTrees

def merge_isotopic_eTrees(list_eTrees):
    '''
    For each pair of eTrees after init_eTrees_by_isopairs, merge if at least one overlap peak is anchor. 
    Ignor if overlap is all at endpoints. I.e. trees can branch, but not converge.
    '''
    new = []
    peak2etree = {}
    for ET in list_eTrees:
        for P in ET.isotope_peaks:
            if P in peak2etree:
                peak2etree[P].append()

    pass



#
# -----------------------------------------------------------------------------
#



def search_formula_mass_db(query_mz, indexed_DB, limit_ppm=10):
    '''
    Find best matched formula_mass in indexed_DB for query_mz within ppm limit.
    
    indexed_DB: # [formula_mass, m/z, charged_formula, selectivity, [relations to neutral formula]]
            example DB_1[99] = 
            [['C2H4NaO3_99.00532', 99.00532046676999, 'C2H4NaO3', 0.9999999571942034, [('C2H4O3_76.016044', 'M+Na[1+]')]], 
            ['C4H5NS_99.013721', 99.01372099999999, 'C4H5NS', 0.9967247109333564, [('C4H5NS_99.01427', 'M[1+]')]], 
            ['C3H8NaO2_99.041706', 99.04170646677, 'C3H8NaO2', 0.999999999998776, [('C3H8O2_76.05243', 'M+Na[1+]')]], ...]

    return
    ------
    Closest match as an entry,
             (match_ppm, [formula_mass, m/z, charged_formula, selectivity, [relations to neutral formula]])
    None if out of ppm limit.
    '''
    _delta = query_mz * limit_ppm * 0.000001
    _low_lim, _high_lim = query_mz - _delta, query_mz + _delta
    result = []
    for ii in range(int(_low_lim), int(_high_lim)+1):
        if ii in indexed_DB:
            for F in indexed_DB[ii]:
                # F[1] is m/z, fixed DB format
                result.append( (abs(query_mz-F[1]), F) )

    result.sort()
    if result and result[0][0] < _delta:
        return result[0][1]
    else:
        return None



def annotate_formula_mass(list_query_mz, indexed_DB, check_mass_accuracy=False, limit_ppm=10):
    '''
    Annotate input m/z list using formula based mass in the indexed_DB,
    which should be generic, not assuming separation of primary ions and their common isotopes and adducts.
    Search efficiency will be improved by the construction strategy of the indexed_DB,
    not dependent on the algorithm in this function. E.g., more frequent features should be in first DB to search.

    If check_mass_accuracy is True, the ppm deviations of each matched primary ion are fitted to a normal distribution,
    to estimate mass_accuracy and ppm_std.

    input
    -----
    input_list: m/z values.
    indexed_DB: per entry - [formula_mass, m/z, charged_formula, ion, selectivity, [relations to neutral formula]]
                This can be a pre-loaded DB, or hot DB in asari.

    # mode='pos', 

    return
    ------
    {mass_accuracy: estimated based on difference from formula value
    ppm_std:         estimated based on precision distribution
    annotated_list: [matched DB entry or None, ...]
    corrected_list_query_mz: m/z values will be corrected if they deviate from formula values by mass_accuracy > 5 ppm.
    }
    '''
    result = list_search_formula_mass_db(list_query_mz, indexed_DB, limit_ppm)
    mass_accuracy, ppm_std, corrected_list_query_mz = None, None, None
    N = len(list_query_mz)
    if check_mass_accuracy:
        list_ppm_errors = []
        for ii in range(N):
            if result[ii]:
                list_ppm_errors.append( 1000000 * (list_query_mz[ii] - result[ii][1])/result[ii][1] )

        mass_accuracy, ppm_std = normal_distribution.fit(list_ppm_errors)
        if abs(mass_accuracy) > 5:   # this is considered significant mass shift, requiring m/z correction for all 
            list_query_mz = [x-x*0.000001*mass_accuracy for x in list_query_mz]
            # redo search because we may find different matches after mass correction
            # Also change search ppm if self.__mass_stdev__ > 5
            result = list_search_formula_mass_db(list_query_mz, indexed_DB, max(10, 2*ppm_std))
            corrected_list_query_mz = list_query_mz
            # update ppm_std because it may be used for downstream search parameters
            list_ppm_errors = []
            for ii in range(N):
                if result[ii]:
                    list_ppm_errors.append( 1000000 * (list_query_mz[ii] - result[ii][1])/result[ii][1] )
            _mu, ppm_std= normal_distribution.fit(list_ppm_errors)

        elif ppm_std > 5:      # no mass correction but need redo search using higher ppm for unmatched mz
            result = list_search_formula_mass_db(list_query_mz, indexed_DB, 2*ppm_std)

    return {
        'mass_accuracy': mass_accuracy,
        'ppm_std': ppm_std,
        'annotated_list': result,
        'corrected_list_query_mz': corrected_list_query_mz,
        }


