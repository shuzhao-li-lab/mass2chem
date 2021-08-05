"""
Grouping features into EmpiricalCompounds, as part of annotation process.
This is ideally done by annotation methods combined with pre-processing.
E.g. similar peak shapes can be taken into account.
Here is a fallback option. When only a feature table is supplied.

Main data structures:

Feature = {'id': '', 'mz': 0, 'rtime': 0, 
            'calibrated_mz': 0, 'calibrated_rtime': 0}

EmpiricalCompound = {
    "neutral_base_mass": 0,
    "neutral_formula": [],
    "list_of_features": [
            {'feature id': '', 'ion': 'M-H[1-]', 'm/z': 169.0013, 'rtime': 55, ..},
            {},
}

E.g. 
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


list(massDict_hmdb_serum.items())[:3] = 
    [('C7H11N3O2_169.085127',
    [['HMDB0000001', '1-Methylhistidine'],
    ['HMDB0000479', '3-Methylhistidine']]),
    ('C3H10N2_74.084398',
    [['HMDB0000002', '1,3-Diaminopropane'],
    ['HMDB0013136', '1,2 Diaminopropane']]),
    ('C4H6O3_102.031694',
    [['HMDB0000005', '2-Ketobutyric acid'],
    ['HMDB0000060', 'Acetoacetic acid']])]

Indexed_features are used internally for search speed.

SL 2021-08-04
"""

from ..massDicts import massDict_hmdb_serum
# this is list of common m/z diff values
from ..common_mass import mass_signatures

import re
import numpy as np
from scipy.optimize import curve_fit


calibration_mass_dict_pos = {'C4H9N3O2_131.069477': 132.07675346677001, 'C6H11NO2_129.078979': 130.08625546677, 'C5H9NO4_147.053158': 148.06043446677, 'C6H10O6_178.047738': 179.05501446677002, 'C9H11NO3_181.073893': 182.08116946677, 'C9H11NO2_165.078979': 166.08625546677, 'C9H8O3_164.047344': 165.05462046677002, 'C23H45NO4_399.334859': 400.34213546677, 'C2H7NO3S_125.014664': 126.02194046676999, 'C5H7NO3_129.042593': 130.04986946677002, 'C5H4N4O3_168.02834': 169.03561646677, 'C5H8O5_148.037173': 149.04444946677, 'C5H12O5_152.068473': 153.07574946677002, 'C9H8O2_148.052429': 149.05970546677, 'C5H10N2O3_146.069142': 147.07641846677, 'C17H33NO4_315.240959': 316.24823546677, 'C5H11NO2S_149.051049': 150.05832546677001, 'C5H10O2_102.06808': 103.07535646676999, 'C15H29NO4_287.209658': 288.21693446677, 'C8H17NO2_159.125929': 160.13320546677002, 'H2O4S_97.967379': 98.97465546676999, 'C6H13NO5_179.079373': 180.08664946677, 'C16H22O4_278.151809': 279.15908546677, 'C20H28O2_300.20893': 301.21620646677, 'C7H8N4O2_180.064726': 181.07200246677002, 'C7H6O2_122.036779': 123.04405546676999, 'C14H22N2O3_266.163043': 267.17031946677, 'C19H16O4_308.104859': 309.11213546677, 'C21H39NO4_369.287909': 370.29518546677, 'C8H6O4_166.026609': 167.03388546677002, 'C18H35NO_281.271865': 282.27914146677, 'C19H37NO4_343.272259': 344.27953546677, 'C26H52NO7P_521.34814': 522.35541646677, 'C7H8N2O2_152.058578': 153.06585446677002, 'C9H17NO2_171.125929': 172.13320546677002, 'C25H47NO4_425.350509': 426.35778546677, 'C8H10O3_154.062994': 155.07027046677, 'C25H45NO4_423.334859': 424.34213546677, 'C22H46NO7P_467.301189': 468.30846546677003, 'C23H48NO7P_481.316839': 482.32411546677, 'C24H48NO7P_493.316839': 494.32411546677, 'C26H54NO7P_523.36379': 524.37106646677, 'C26H48NO7P_517.316839': 518.32411546677, 'C28H52NO7P_545.34814': 546.35541646677, 'C28H50NO7P_543.332489': 544.33976546677, 'C24H50NO6P_479.337575': 480.34485146677, 'C25H52NO7P_509.34814': 510.35541646677, 'C21H38O4_354.27701': 355.28428646677, 'C17H31NO4_313.225308': 314.23258446677, 'C15H27NO4_285.194008': 286.20128446677, 'C19H35NO4_341.256609': 342.26388546677003, 'C22H37NO3_363.277344': 364.28462046677004, 'C19H22ClN5O_371.151288': 372.15856446677003, 'C24H30O6_414.204239': 415.21151546677, 'C21H25N_291.1987': 292.20597646677, 'C15H23NO2_249.172879': 250.18015546677, 
    'C17H27NO4_309.194008': 310.20128446677, 'C21H25ClN2O3_388.15537': 389.16264646677, 'C16H18O10_370.089997': 371.09727346677}


MASS_RANGE = (50, 2000)

# fraction of total retention time, or of ranks of retention time
# used to determine coelution of ions ad hoc
RETENTION_TIME_TOLERANCE_FRAC = 0.02    


def __read_features__(feature_table, 
                        id_col=False, mz_col=1, rtime_col=2, 
                        intensity_cols=(4,10), delimiter="\t"):
    '''
    Read a text feature table into a list of features.
    Internal use.

    To-do: should switch to pd.read_table??

    Input
    -----
    feature_table: Tab delimited feature table. First line as header.
                    Recommended col 0 for ID, col 1 for m/z, col 2 for rtime.
    id_col: column for id. If feature ID is not given, row_number is used as ID.
    mz_col: column for m/z.
    rtime_col: column for retention time.
    intensity_cols: range of columns for intensity values.

    Return
    ------
    List of features: [{'id': '', 'mz': 0, 'rtime': 0, 
                        intensities: [], 'representative_intensity': 0, ...}, 
                        ...], 
                        where representative_intensity is mean value.
    '''
    featureLines = open(feature_table).read().splitlines()
    header = featureLines[0].split(delimiter)
    num_features = len(featureLines)-1
    # sanity check
    print("table headers ordered: ", header[mz_col], header[rtime_col])
    print("Read %d feature lines" %num_features)
    L = []
    for ii in range(1, num_features+1):
        if featureLines[ii].strip():
            a = featureLines[ii].split(delimiter)
            if isinstance(id_col, int):         # feature id specified
                iid = a[id_col]
            else:
                iid = 'row'+str(ii)
            xstart, xend = intensity_cols
            intensities = [float(x) for x in a[xstart: xend]]
            L.append({
                'id': iid, 'mz': float(a[mz_col]), 'rtime': float(a[rtime_col]),
                'intensities': intensities,
                'representative_intensity': np.mean(intensities),
            })
    return L


def __index_features__(list_of_features, range_low=20, range_high=2000):
    '''
    Indexing a list of features to speed up search.
    A filter/Exception should be added for the mass range in production.
    Internal use.

    Input
    -----
    List of features, [{'id': '', 'mz': 0, 'rtime': 0, ...}, ...]
    range_low, range_high: m/z range, default [20, 2000]. 

    Return
    ------
    Dictionary of features, indexed using integer of m/z.
    '''
    L = {}
    for ii in range(range_low, range_high+1): L[ii] = []
    for F in list_of_features: L[int(F['mz'])].append(F)
    return L


def mass_calibrate(indexed_features, calibration_mass_dict, limit_ppm=25):
    '''
    Use formula based reference mass list to calibrate m/z values in list_of_features.
    Assuming a large number of compounds in ref_mass_dict exist in the list_of_features,
    which should be the case in common biological matrices.
    A regression model is built for (theoretical m/z) = a * (measured m/z) + b,
    then the model is applied to all features to calculate a calibrated mz.

    To-do: add 99% interval of ppm from fitting, and use for the whole expt.
    
    Input
    -----
    indexed_features: List of Feature instances, from the same Experiment.
    calibration_mass_dict: list of known formulae and corresponding m/z values. 
                E.g. in positive ESI {'C7H11N3O2_169.085127': 170.092403, 
                'C3H10N2_74.084398': 75.091674, ...}
    
    limit_ppm: 25. Expecting a high-res instrument performs within 25 ppm, 
                a should be 0.000025 within 1, and b should be < 0.01.
                Otherwise, raise flag.

    Return
    ------
    (a, b): coefficients in calibration model.
    updated_indexed_features: {int_mz: [{'id': '', 'mz': 0, 'calibrated_mz': 0, 'rtime': 0, ...}, ...],
                                ...}
    '''
    def find_closest(query_mz, indexed_features, limit_ppm):
        # to find closest match of a theoretical m/z in indexed_features, under limit_ppm
        relevant_features = indexed_features[int(query_mz)-1] + indexed_features[int(query_mz)] + indexed_features[int(query_mz)-1]
        result = [(abs(query_mz-F['mz']), F['mz'], query_mz) for F in relevant_features]
        result.sort()
        if result[0][0] < result[0][1] * limit_ppm * 0.000001:
            return result[0]
        else:
            return None

    def objective(x, a, b): return a * x + b

    matched = []
    for m in calibration_mass_dict.values():
        best_pair = find_closest(m, indexed_features, limit_ppm)
        if best_pair:
            matched.append(best_pair)

    num_matched = len(matched)
    if num_matched < 3:
        print("Aborted. Found %d matched mass features." %num_matched)

    else:
        if num_matched < 10:
            print("Warning. Using %d matched mass features, too few for regression." %num_matched)
        X, Y = [x[1] for x in matched], [x[2] for x in matched]
        popt, _ = curve_fit(objective, X, Y)
        a, b = popt
        if abs(b) > 0.01 or abs(a-1) > 0.000025:
            print(" -- WARNING -- extreme parameters are usually bad, ", a, b)

        updated_indexed_features = {}
        for k,v in indexed_features.items():
            updated_indexed_features[k] = []
            for F in v:
                F['calibrated_mz'] = a * F['mz'] + b
                updated_indexed_features[k].append(F)

        return (a,b), updated_indexed_features



def __flag_contaminants__(indexed_features, contaminant_list, ppm=2):
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
            F['potential_contaminants'] = check_contaminant(F['mz'], contaminant_list, ppm)

    return indexed_features


def __make_anchordict__(massDict_hmdb_serum, mode='pos', MASS_RANGE=(50, 2000)):
    '''
    To build anchor mz list based on keys like 'C7H11N3O2_169.085127'.
    This core list is used to separate ions from a feature list to reduce search space.
    The matched results will be organized by the same key of formula_mass.

    Return
    ------
    {'C7H11N3O2_169.085127': [(170.09240346677, '+H[1+]'), (192.074345, 'Na'), (169.085127, 'M')], ...}
    '''
    posList = [(1.00727646677, '+H[1+]'), (22.989218, 'Na'), (0, 'M')]
    negList = [(-1.00727646677, '-H[-]'), (34.969402, 'Cl-'), (0, 'M')]
    if mode == 'pos': useList = posList
    elif mode == 'neg': useList = negList
    else: print ("ionization mode is either pos or neg.")

    new = {}
    for k in massDict_hmdb_serum:
        kmass = float(k.split('_')[1])
        new[k] = [(x+kmass, y) for x,y in useList if MASS_RANGE[0] < x+kmass < MASS_RANGE[1]]
    return new

#
# Annotate, round 0 - flag potential contaminants
# round 1 - get anchor ions by formula matching; use HMDB serum for now
# round 2 - get empCdps corresponding to anchorList
# round 3 - extend major isotope and adducts, neutral loss, fragments (in source) for core empCdps
# round 4 - search unmatched_features in larger DB

# round 5 - extended adducts search

# round 6 - ab initio empCpd grouping

# 

# anchorList = __make_anchordict__(massDict_hmdb_serum, mode='pos')

def __find_mz_indexed_features__(query_mz, indexed_features, limit_ppm=2):
    # to find all matches of a m/z in indexed_features, under limit_ppm
    matched = []
    _delta = query_mz * limit_ppm * 0.000001
    low_lim, high_lim = query_mz - _delta, query_mz + _delta
    for ii in range(int(low_lim), int(high_lim)+1):
        for F in indexed_features[ii]:
            if low_lim < F['mz'] < high_lim:
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

    print(len(matched_features), len(set(matched_features)), len(unmatched_features))
    return matched_empCpds, unmatched_features


def extend_annotate_anchorList(list_empCpds, unmatched_features, mode='pos',  MASS_RANGE=(50, 2000), ppm=2):
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


    unmatched_features = __index_features__(unmatched_features)
    return __extend_search__(list_empCpds, unmatched_features, mode,  MASS_RANGE, ppm)

def __extend_search__(list_empCpds, indexed_unmatched_features, mode='pos',  MASS_RANGE=(50, 2000), ppm=2):
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
        for target_mz, ion in __compute_adducts__(E['neutral_base_mass'], E['neutral_formula'], mode) :
            if MASS_RANGE[0] < target_mz < MASS_RANGE[1]:
                found_features = __find_mz_indexed_features__(target_mz, indexed_unmatched_features, ppm)
                for F in found_features: 
                    _tmp_matched.append(F['id'])
                    F['ion'] = ion
                    E['list_of_features'].append(F)

    updated_unmatched_features = []
    for k,v in indexed_unmatched_features.items():
        for F in v:
            if F['id'] not in _tmp_matched:
                updated_unmatched_features.append(F)

    print(len(_tmp_matched), len(set(_tmp_matched)), len(updated_unmatched_features))
    return list_empCpds, updated_unmatched_features



def __compute_adducts__(mw, cFormula,  mode='pos'):
    '''
    Calculating isotopes and adducts.
    C13, S34, Cl37 pos only applicable if the atom is in formula,
    and others in `required subgroup`, as the last item in the tuples

    The isotopes, neutral losses and adducts can have multiplications,
    but the 3rd round search will take care of those.

    Return
    -------
    List of adducts: e.g. [(58.53894096677, 'M+2H[2+]'), (39.36171946677, 'M+3H[3+]'), (117.07400546676999, 'M(C13)+H[1+]'),]

    '''

    def __parse_chemformula__(x):
        '''This does not deal with nested groups in chemical formula.
        Formula from HMDB 3.5 are compatible.
        '''
        p = re.findall(r'([A-Z][a-z]*)(\d*)', x)
        d = {}
        for pair in p: d[pair[0]] = int( pair[1] or 1 )
        return d

    def __check_sub__(Fm1, Fm2):
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

    PROTON = 1.00727646677
    if mode == 'pos': 
        addList = [
            #(mw, 'M[1+]', ''),                              # already in anchor empCpd
            #(mw + PROTON, 'M+H[1+]', ''),
            (mw/2 + PROTON, 'M+2H[2+]', ''),
            (mw/3 + PROTON, 'M+3H[3+]', ''),
            (mw +1.0034 + PROTON, 'M(C13)+H[1+]', 'C'),
            (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]', 'C'),
            (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]', 'C'),
            (mw +1.9958 + PROTON, 'M(S34)+H[1+]', 'S'),
            (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]', 'Cl'),
            #(mw + 21.9820 + PROTON, 'M+Na[1+]', ''),        # Na = 21.9820 + PROTON = 22.9893
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
            ]

    elif mode == 'neg': 
        addList = [
            #(mw - PROTON, 'M-H[-]', ''),                   # already in anchor empCpd
            (mw/2 - PROTON, 'M-2H[2-]', ''),
            (mw + 1.0034 - PROTON, 'M(C13)-H[-]', 'C'),
            (mw + 1.9958 - PROTON, 'M(S34)-H[-]', 'S'),
            (mw + 1.9972 - PROTON, 'M(Cl37)-H[-]', 'Cl'),
            (mw + 22.9893 - 2*PROTON, 'M+Na-2H[-]', ''),
            (mw + 38.9628 - 2*PROTON, 'M+K-2H[-]', ''),
            (mw - 18.0106 - PROTON, 'M-H2O-H[-]', 'H2O'),
            #(mw + 34.9689, 'M+Cl[-]', ''),
            (mw + 36.9659, 'M+Cl37[-]', ''),
            (mw + 78.9183, 'M+Br[-]', ''),
            (mw + 80.9163, 'M+Br81[-]', ''),
            (mw + 2*12 + 3*1.007825 + 14.00307 - PROTON, 'M+ACN-H[-]', ''),
            (mw + 1.007825 + 12 + 2*15.99491, 'M+HCOO[-]', ''),
            (mw + 3*1.007825 + 2*12 + 2*15.99491, 'M+CH3COO[-]', ''),
            (mw - PROTON + 15.99491, 'M-H+O[-]', ''),
            ]

    dict_cFormula = __parse_chemformula__(cFormula)
    adducts2get = []
    for x in addList:
        if __check_sub__(dict_cFormula, __parse_chemformula__(x[2])):
            adducts2get.append(x[:2])

    return adducts2get










    


























        
def __is_coelution__(self, massFeature1, massFeature2):
    '''
    True if retention times are within a tolerance window in time or ranks.
    Not assuming massFeatures are sorted in this function.
    '''
    if abs(massFeature1.retention_time - massFeature2.retention_time) < self.rtime_tolerance or \
        abs(massFeature1.retention_time_rank - massFeature2.retention_time_rank) < self.rtime_tolerance_rank:
        return True
    else:
        return False


def __is_mz_match__(mz1, mz2, mz_tolerance):
    if abs(mz1 - mz2) < mz_tolerance:
        return True
    else:
        return False



