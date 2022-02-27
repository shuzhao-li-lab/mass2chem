'''
Two lines of search methods:
- A centurion_tree is an indexed dictionary of peaks/features.
- We also use DataFrame based vector operations for some cases (less often).

Emperical compounds are constructed by co-eluting isotopic and adduct patterns.
# isotopic_signatures: example as [(182, 'anchor'), (191, 'M(13C)'), (205, 'M(18O)')]
# consider Cl as adduct, both 35Cl and 37Cl are considered
# 34S and 37Cl are close, but should have diff ratios and are formula dependent

'15N/14N' and '33S/32S' are close for high m/z features, e.g. 0.002353/470 ~ 5 pmm.
'18O/16O' and 'M(13C),M(15N)' are close for high m/z features, e.g. 0.003855/770 ~ 5 pmm.

Comprehensive combinations of isotopes and adducts, vs seeded search:
all adducts can be combined with isotopes in theory, greatly depending on abundance in practice.
Controlled by coelution, it's high confidence in assigning these peaks to empCpds, 
but not 100% correct. 
Until we have a better strategy, step-wise searches are less prone to errors.

search_patterns: Use seed_empCpd_patterns. May update after seeing experimental statistics.
                    Avoid H2O because hard to distinguish from differences in chemical formulae.
                    The third item in tuples is to control abundance ratio, 
                    which is linient to leave room for high carbon numbers and measurement.

'''

# Mass difference and Abundance of natually occuring isotopic elements.
# (mz difference, notion, ratio low limit, ratio high limit)
# Ratio is loose to leave room for 1) multiple atoms, e.g. 40 C atoms lead to ~ 40*1% in abundance
# 2) measurement errors; not all peaks are in linear dynamic range
# 3) quantification by peak height is not same as peak area

isotopic_patterns = [
    # mass diff, isotopes, (intensity ratio constraint)
    (1.003355, '13C/12C', (0, 0.8)),      # 13C-12C, 12C~99%, 13C ~ 1%
    (0.997035, '15N/14N', (0, 0.2)),     # 15N-14N, 14N ~ 99.64%, 15N ~ 0.36%
    (2.004245, '18O/16O', (0, 0.2)),      # 18O-16O, 16O ~ 99.76, 16O ~ 0.2%
    (1.995796, '34S/32S', (0, 0.4)),      # 32S (95.02%), 33S (0.75%), 34S (4.21%)
    (0.999388, '33S/32S', (0, 0.1)),
    # double isotopes
    (2.00039, 'M(13C),M(15N)', (0, 0.2)),
    (2.999151, 'M(13C),M(34S)', (0, 0.4)),
    # double charged
    (0.5017, '13C/12C, double charged', (0, 0.8)),
    (0.4985, '15N/14N, double charged', (0, 0.2)),
]

common_adducts = {
    # mass diff, modification
    # not using (intensity ratio constraint), but it can be documented or learned
    'pos': [
        (1.0078, 'H'),
        (21.9820, 'Na/H'), # Na replacing H
        (10.991, 'Na/H, double charged'),
        (18.0106, '+H2O'), 
        (18.033823, '+NH4'),
        (37.9559, '39K/H'),
        (39.9540, '41K/H'),
        (41.026549, 'Acetonitrile'),
    ],
    'neg': [
        (1.0078, 'H'),
        (22.9893, 'Na'),
        (20.97474706646, '+Na-2H'),
        (18.0106, 'H2O'), 
        (34.9689, '35Cl'),
        (36.9659, '37Cl'),
        (40.01926853323, '+ACN-H'),
        (44.998201, 'COOH'),
        (59.013295, 'CH3COO'),
    ],
}

extended_adducts = {
    'pos': [],
    'neg': [],
}

# this does not include combinatorial values of isotopes and adducts
# Mostly used as seeds; better to do inclusive search after matched to formulae
seed_empCpd_patterns = {
    'pos': [(1.003355, '13C/12C', (0, 0.8)), (1.0078, 'H'), (21.9820, 'Na/H'), ],
    'neg': [(1.003355, '13C/12C', (0, 0.8)), (1.0078, 'H'), (20.97474706646, '+Na-2H'), (34.9689, '35Cl')],
}



#
# -----------------------------------------------------------------------------
#

def build_centurion_tree(list_peaks):
    '''
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...]
    Return a dictionary, indexing mzList by 100*mz bins.
    Because most high-resolution mass spectrometers measure well under 0.01 amu, 
    one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).
    list_mass_tracks has similar format as list_peaks.
    '''
    d = {}
    for p in list_peaks:
        cent = int(100 * p['mz'])
        if cent in d:
            d[cent].append(p)
        else:
            d[cent] = [p]
    return d


def build_peak_id_dict(list_peaks):
    d = {}
    for p in list_peaks:
        d[p['id_number']] = p
    return d

def __build_centurion_tree_mzlist(mzList):
    '''
    Return a dictionary, indexing mzList by 100*mz bins.
    Because most high-resolution mass spectrometers measure well under 0.01 amu, 
    one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).
    '''
    d = {}
    for ii in range(len(mzList)):
        cent = int(100*mzList[ii])
        if cent in d:
            d[cent].append((mzList[ii], ii))
        else:
            d[cent] = [(mzList[ii], ii)]
    return d

def find_all_matches_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=5):
    '''
    Return matched peaks in mz_centurion_tree.
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    results = []
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            if abs(peak['mz']-query_mz) < mz_tol:
                results.append(peak)
                
    return results

def find_best_match_centurion_indexed_list(query_mz, mz_centurion_tree, limit_ppm=2):
    '''
    Return matched indices in mz_centurion_tree (based on peak list).
    '''
    q = int(query_mz * 100)
    mz_tol = query_mz * limit_ppm * 0.000001
    result = (None, 999)
    for ii in (q-1, q, q+1):
        L = mz_centurion_tree.get(ii, [])
        for peak in L:
            _d = abs(peak['mz']-query_mz)
            if _d < min(result[1], mz_tol):     # enforce mz_tol here
                result = (peak, _d)
                
    return result[0]

def is_coeluted_by_overlap(P1, P2, rt_tolerance=10):
    '''
    coelution is defined by verlap more than half of the smaller peak, or apexes within rt_tolerance.
    If not enough parameters are given for peaks, fallback to is_coeluted_by_distance(P1, P2, rt_tolerance=10).
    
    Example peak format: {'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0,
    'left_base': 648, 'right_base': 655, 'id_number': 555}

    rt_tolerance: tolerance of retention time is used when peak boundaries are not given.
                Can be scan numbers or seconds.
                
    return: True or False
    '''
    _coeluted = False
    try:
        len1, len2 = P1['right_base'] - P1['left_base'], P2['right_base'] - P2['left_base']
        # overlap is the max L and min R
        overlap = min(P1['right_base'], P2['right_base']) - max(P1['left_base'], P2['left_base'])
        if overlap > 0.5 * min(len1, len2):
            _coeluted = True
        return _coeluted
    except:
        print("KeyError, fallback to function is_coeluted_by_distance.")
        return is_coeluted_by_distance(P1, P2, rt_tolerance=10)

def is_coeluted_by_distance(P1, P2, rt_tolerance=10):
    _coeluted = False
    if abs(P1['apex']-P2['apex']) <= rt_tolerance:
        _coeluted = True
    return _coeluted

def get_seed_empCpd_signatures(list_peaks, 
                    mztree, 
                    search_patterns = seed_empCpd_patterns['pos'],
                    mz_tolerance_ppm=5, 
                    isotope_rt_tolerance=5, coelution_rt_tolerance=10, 
                    is_coeluted = is_coeluted_by_overlap,
                    check_isotope_ratio = True,
                    ):
    '''
    To find a core group of ions that belong to an empirical compound (emp_cpd, or empCpd, ecd).
    This should be a combination of 12C/13C isotopes and common adducts;
    an empCpd can be constructed even if the 13C ion is missing.
    Other isotopes and adducts can be searched after establishing the seed empCpds.

    This is currently a heuristic approach to generate empCpds without too many errors.
    One may research more to polish an optimal method of assigning ions to empCpds, 
    where chemical formulae and learned patterns can be incorporated.

    Input
    =====
    list_peaks: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...]
    mztree: indexed list_peaks
    mz_tolerance_ppm: ppm tolerance in examining m/z patterns.
    seed_search_patterns: initial ions for constructing an empCpd. 
            This is very conservative and requires stricter isotope_rt_tolerance, because other ions can be searched later. 
            Can use self.seed_search_patterns, e.g. for pos ions:
            [(1.003355, '13C/12C', (0, 0.8)), (1.0078, 'H'), (21.9820, 'Na/H')]
    isotope_rt_tolerance: tolerance threshold for deviation in retetion time, stricter version for seed ions (isotopes). 
            Default unit is scan number, but can be arbitrary.
    coelution_rt_tolerance: tolerance threshold for deviation in retetion time of adducts (other than seeds).
    is_coeluted: coelution function

    Return
    ======
    list of lists of peak numbers that match search_patterns patterns, e.g.
    [ [(195, 'anchor'), (206, '13C/12C')],  ...]
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [ (P1['id_number'], 'anchor'), ]          # Nature is nice to have lowest mass for the most abundant 
        for _pair in search_patterns:
            (mass_difference, relation) = _pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if abs(P1['apex']-P2['apex']) <= isotope_rt_tolerance and is_coeluted(P1, P2, coelution_rt_tolerance):
                    if check_isotope_ratio and len(_pair) > 2:  # checking abundance ratio
                        (abundance_ratio_min, abundance_ratio_max) = _pair[2]
                        if abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                            matched.append( (P2['id_number'], relation) )
                    else:
                        matched.append( (P2['id_number'], relation) )
 
        if len(matched) > 1:
            # if len(matched) > 2: matched = merge_confused_peak_matches(matched)
            signatures.append(matched)

    return signatures


def extend_seed_empCpd_signature(seed_signatures, peak_dict, mztree, ext_search_patterns, 
                    mz_tolerance_ppm=5, 
                    coelution_rt_tolerance=10, 
                    is_coeluted = is_coeluted_by_overlap,
                    ):
    '''
    Extend a single signature by ext_search_patterns, including combinatorial values of isotopes and adducts.

    Input
    =====
    mztree: indexed list_peaks
    seed_signatures: e.g. [(182, 'anchor'), (191, '13C/12C'), (205, '18O/16O')]
    search_patterns: e.g. [(10.991, 'Na/H, double charged'), (18.0106, '+H2O'), (18.033823, '+NH4'),]
    example peak: {'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}

    Return
    ======
    One epd (empirical compound) based on seed peaks and ext_search_patterns, as a list of [(peak_id, relation), ...]
    '''
    matched = []
    for seed in seed_signatures:
        peak, relation = peak_dict[seed[0]], seed[1]
        for adduct in ext_search_patterns:
            tmp = find_all_matches_centurion_indexed_list( peak['mz'] + adduct[0], mztree, mz_tolerance_ppm )
            for P2 in tmp:
                if is_coeluted(peak, P2, coelution_rt_tolerance):
                    matched.append( (P2['id_number'], relation +','+ adduct[1]) )

    return matched


def find_isotopic_signatures(list_peaks, mztree, isotopic_patterns, mz_tolerance_ppm=5, rt_tolerance_scans=5):
    '''
    See find_isotopic_pairs. This extends to all related isotopic signatures.
    Allows ambiguous matches, which will be dealt with in empCpd statistics.
    But the peaks need to be unique. E.g. 
    '15N/14N' and '33S/32S', '18O/16O' and 'M(13C),M(15N)' are close for high m/z features.
    Matches like ('F5486', '15N/14N'), ('F5486', '33S/32S') need to merge on the peak/feature.

    Input
    =====
    isotopic_patterns = [(1.003355, '13C/12C', (0, 0.8)), ...], the third item is optional limits of abundance ratio.

    Return
    ======
    list of lists of peak numbers that match isotopic patterns, [ [], ... ]

    Example
    =======
    [ [(195, 'anchor'), (206, '13C/12C')], 
      [(182, 'anchor'), (191, '13C/12C'), (205, '18O/16O')],
      [(295, 'anchor'), (335, '13C/12C'), (368, 'M(13C),M(34S)')], ...]
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [ (P1['id_number'], 'anchor'), ]          # Nature is nice to have lowest mass for the most abundant 
        for isotopic_pair in isotopic_patterns:
            (mass_difference, relation) = isotopic_pair[:2]
            tmp = find_all_matches_centurion_indexed_list(P1['mz'] + mass_difference, mztree, mz_tolerance_ppm)
            for P2 in tmp:
                if abs(P1['apex']-P2['apex']) <= rt_tolerance_scans and is_coeluted_by_overlap(P1, P2):
                    if len(isotopic_pair) == 2:             # not checking abundance ratio
                        matched.append( (P2['id_number'], relation) )
                    else:
                        (abundance_ratio_min, abundance_ratio_max) = isotopic_pair[2]
                        # checking abundance ratio
                        if abundance_ratio_min*P1['height'] < P2['height'] < abundance_ratio_max*P1['height']:
                            matched.append( (P2['id_number'], relation) )
        if len(matched) > 1:
            # if len(matched) > 2: matched = merge_confused_peak_matches(matched)
            signatures.append(matched)

    return signatures


def find_adduct_signatures(list_peaks, mztree, adduct_patterns, mz_tolerance_ppm=5):
    '''
    Search adduct mass_diff in ceelution peaks, de novo. 
    Not requiring matched apex as in isotopic pairs, nor abundance ratio restriction.
    '''
    signatures = []
    for P1 in list_peaks:
        matched = [ (P1['id_number'], 'anchor'), ]
        for adduct in adduct_patterns:
            tmp = find_all_matches_centurion_indexed_list( P1['mz'] + adduct[0], mztree, mz_tolerance_ppm )
            for P2 in tmp:
                if is_coeluted_by_overlap(P1, P2):
                    matched.append( (P2['id_number'], adduct[1]) )

        if len(matched) > 1:
            signatures.append(matched)

    return signatures


def merge_confused_peak_matches(L):
    '''
    merge if same peak is matched to more than one isotopes or adducts.
    L example: [(295, 'anchor'), (335, '13C/12C'), (368, 'M(13C),M(34S)')]
    '''
    tmp, _d = [L[0]], {}
    for peak in L[1:]:
        if peak[0] in _d:
            _d[peak[0]] += ';' + peak[1]
        else:
            _d[peak[0]] = peak[1]
    for k,v in _d.items():
        tmp.append((k,v))
    return tmp


def find_mzdiff_pairs_from_masstracks(list_mass_tracks, list_mz_diff=[1.003355, 21.9820], mz_tolerance_ppm=5):
    '''
    Find all pairs in list_mass_tracks that match a pattern in list_mz_diff, and return their id_numbers as pairs.
    This function does not use coeluction (rtime) rules. 

    Input
    =====
    list_mass_tracks: [{ 'id_number': ii,  'mz': xic[0], 'rt_scan_numbers': xic[1],  'intensity': xic[2],  }, ...]
    list_mz_diff: defaul 1.003355, 21.9820 are 13C/12C, Na/H. Dependent on ionization mode.
    
    Return
    ======
    list of pairs of mass tracks numbers.
    '''
    pairs = []
    # list_mass_tracks has similar format as list_peaks.
    mztree = build_centurion_tree(list_mass_tracks)
    for mzdiff in list_mz_diff:
        for P1 in list_mass_tracks:
            P2 = find_best_match_centurion_indexed_list(P1['mz'] + mzdiff, mztree, mz_tolerance_ppm)
            if P2:
                pairs.append((P1['id_number'], P2['id_number']))

    return pairs


#
# -----------------------------------------------------------------------------
#

def search_formula_mass_dataframe(query_mz, DFDB, limit_ppm=10):
    '''
    return best match formula_mass in DFDB if under ppm limit.
    DFDB is using a Pandas DataFrame to house reference database.
    ppm is signed to capture the direction of mass shift.

    Not using idxmin,   #ii = DFDB.tmp.idxmin()
    # in pd.DF, index not necessarily integer; can be sequence if more than one match, but they are trying to fix in pandas dev version
    '''
    DFDB['tmp'] = abs(DFDB.mz - query_mz)
    ii = DFDB.tmp.values.argmin()
    ppm = 1000000 * (query_mz - DFDB.iloc[ii].mz)/query_mz
    if  abs(ppm) < limit_ppm:
        return (DFDB.iloc[ii].name, ppm)            # name is formula_mass
    else:
        return None

