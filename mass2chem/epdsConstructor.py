'''
Constructing empirical compounds de novo.

The ion patterns here are constructed differently from formula based calculations,
because formulae are not known here.

        Anchor peak is the most abundant isotopic peak. Modification peaks are from adducts, neutral loss and fragments.
        Requires at least one of isotopic_signatures and adduct_signatures. Use peak numbers only.
        If initiated by isotopic_signatures, do adduct search;
        elif initiated by adduct_signatures, assuming other primary isotopes are not seen.
        
        We build indexed centurion trees here to assist searches.
        Minor possibility that peak:empCpd is not N:1.

'''
from .search import *

class epdsConstructor:
    '''
    To organize a table of peaks/features into a list of empirical compounds
    (https://github.com/shuzhao-li/metDataModel).
    empCpds (or epds) can be initiated by signatures on isotopic relationships or adduct relationships.
    For low-intensity peaks, their isotopic counterparts may not be detectable.
    Only common adducts are considered at this step.
    This list of epds are used in m/z alignment/correspondence.
    For assigned peaks/features, will calculate selectivity/score later too.

    Input
    =====
    isotopic_patterns: [(1.003355, '13C/12C', (0, 0.2)), ...],
    peak_list: [{'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
    'left_base': 648, 'right_base': 655, 'id_number': 555},
                 ...]
    
    Example
    =======
    To get list of empCpds, run
    ECCON = epdsConstructor(list_peaks)
    list_empCpds = ECCON.peaks_to_epds()
    # e.g. {'id': 358, 'list_peaks': [(4215, 'anchor'), (4231, '13C/12C'), (4339, 'anchor,+NH4')]},

    To add -
    support of user input formats where rtime isn't precise or unavailable.

    '''
    def __init__(self, peak_list, mode='pos'):
        '''
        self.peak_dict: {peak_id: {'parent_masstrace_id': 1670, 'mz': 133.09702315984987, 'apex': 654, 'height': 14388.0, 
                    'left_base': 648, 'right_base': 655, 'id_number': 555}, ...}
        '''
        # will add parameters dict
        self.peak_list = peak_list
        self.peak_dict = build_peak_id_dict(self.peak_list)
        self.mode = mode
        self.seed_search_patterns = seed_empCpd_patterns[self.mode]
        self.ext_search_patterns = [x for x in common_adducts[self.mode] if x not in self.seed_search_patterns]

    def peaks_to_epdDict(self, 
                            seed_search_patterns, 
                            ext_search_patterns,
                            mz_tolerance_ppm=5, 
                            check_coelution=True,
                            isotope_rt_tolerance=5,     # default is scan number
                            coelution_rt_tolerance=10, 
                            coelution_function='overlap',
                            check_isotope_ratio = True):
        '''
        Wrapper function of constructing empCpds (empricial compounds).
        This complies with empCpd JSON convention, while self.peaks_to_epds produces short-handed lists.
        Input
        =====
        mz_tolerance_ppm: ppm tolerance in examining m/z patterns.
        seed_search_patterns: initial ions for constructing an empCpd. 
                This is very conservative and requires stricter isotope_rt_tolerance, because other ions can be searched later. 
                Can use self.seed_search_patterns, e.g. for pos ions:
                [(1.003355, '13C/12C', (0, 0.8)), (1.0078, 'H'), (21.9820, 'Na/H')]
        ext_search_patterns: 
                Additional ions to be searched to extend seed empCpds. Can use self.ext_search_patterns.
        check_coelution: True unless retention time is not available.
        isotope_rt_tolerance: tolerance threshold for deviation in retetion time, stricter version for seed ions (isotopes). 
                Default unit is scan number, but can be arbitrary.
        coelution_rt_tolerance: tolerance threshold for deviation in retetion time of adducts (other than seeds).
        coelution_function: default is 'overlap', checking overlap of peaks. Fallback is by distance, if any other value is given.

        Return
        ======
        list_empCpds, e.g. [{'interim_id': 12,
                        'neutral_formula_mass': None,
                        'neutral_formula': None,
                        'Database_referred': [],
                        'identity': [],
                        'MS1_pseudo_Spectra': [{'id_number': 'F94',
                            'mz': 98.97531935309317,
                            'apex': 700.0,
                            'ion_relation': 'anchor',
                            'parent_epd_id': 12},
                        {'id_number': 'F8161',
                            'mz': 99.9786892845517,
                            'apex': 701.0,
                            'ion_relation': '13C/12C',
                            'parent_epd_id': 12},
                        {'id_number': 'F317',
                            'mz': 116.98587698692174,
                            'apex': 700.0,
                            'ion_relation': 'anchor,+H2O',
                            'parent_epd_id': 12}],
                        'MS2_Spectra': []}, ...]
        '''
        return self.index_reformat_epds( 
                        self.peaks_to_epds(
                        seed_search_patterns, ext_search_patterns, mz_tolerance_ppm, 
                        check_coelution, isotope_rt_tolerance, coelution_rt_tolerance, coelution_function, check_isotope_ratio) 
                        )

    def peaks_to_epds(self, seed_search_patterns, ext_search_patterns,
                            mz_tolerance_ppm=5, 
                            check_coelution=True,       # use later
                            isotope_rt_tolerance=5, coelution_rt_tolerance=10, 
                            coelution_function='overlap',
                            check_isotope_ratio=True):
        '''
        Seeds by core set of 12C/13C isotopes and common adducts; then extend to common adducts.
        exhaustive_search is discouraged here, because it can be done better when formulae are supplied by a reference database.
        See parameter descriptions in peaks_to_epdDict.

        Return
        ======
        list_empCpds, [{'id': ii, 'list_peaks': [(peak_id, ion), (), ...],}, ...]
            The first peak is anchor ion.
        '''
        print("\n\nAnnotating empirical compounds on %d features/peaks, ..." %len(self.peak_list))
        if coelution_function == 'overlap':
            is_coeluted = is_coeluted_by_overlap
        else:
            is_coeluted = is_coeluted_by_distance
            # usage e.g. is_coeluted(P1, P2, coelution_rt_tolerance=10)

        list_empCpds = []
        epds = self.build_epds_2_steps(self.peak_list, seed_search_patterns, ext_search_patterns, 
                        mz_tolerance_ppm, 
                        isotope_rt_tolerance, coelution_rt_tolerance, is_coeluted, check_isotope_ratio)
        found2 = []
        for L in epds:
            found2 += [x[0] for x in L]
        found2 = set(found2)
        _NN2, ii = len(found2), -1
        for ii in range(len(epds)):
            list_empCpds.append(
                {'id': ii, 'list_peaks': epds[ii]}
            )
        print("epdsConstructor - numbers of seeded epds and included peaks: ", (ii+1, _NN2))
        return list_empCpds

    def build_epds_2_steps(self, peak_list, search_patterns, ext_search_patterns, mz_tolerance_ppm, 
                                isotope_rt_tolerance, coelution_rt_tolerance, is_coeluted, check_isotope_ratio):
        '''
        return list of epds by 2-step search: 


        '''

        subnetworks, peak_dict, edge_dict = peaks_to_networks(peak_list, 
                    isotope_search_patterns = isotope_search_patterns,
                    adduct_search_patterns = adduct_search_patterns,
                    mz_tolerance_ppm=mz_tolerance_ppm, 
                    rt_tolerance = rt_tolerance)



        epdTrees = []
        for subnet in subnetworks:
            E = empTree(subnet)
            epdTrees.append(
                E.to_tree()
            )







        return epdTrees

 


    def __old__build_epds_2_steps(self, peak_list, search_patterns, ext_search_patterns, mz_tolerance_ppm, 
                                isotope_rt_tolerance, coelution_rt_tolerance, is_coeluted, check_isotope_ratio):
        '''
        return list of epds by 2-step search: seeds first, then extended by another set of search pattern.
        search_patterns, e.g. [(1.003355, '13C/12C', (0, 0.8)), (10.991, 'Na/H, double charged'), (18.0106, '+H2O')]
        '''
        epds, found = [], []
        mztree = build_centurion_tree(peak_list)
        seeds = get_seed_empCpd_signatures(peak_list, mztree, search_patterns, mz_tolerance_ppm, 
                        isotope_rt_tolerance, coelution_rt_tolerance, is_coeluted, check_isotope_ratio)
        # [[(182, 'anchor'), (191, '13C/12C'), (205, 'H')], ...]
        for L in seeds:
            found += [x[0] for x in L]
        found = set(found)
        remaining_peaks = [P for P in peak_list if P['id_number'] not in found]
        mztree = build_centurion_tree(remaining_peaks)
        for G in seeds:
            epds.append(
                G + extend_seed_empCpd_signature(G, self.peak_dict, mztree, ext_search_patterns, 
                        mz_tolerance_ppm, coelution_rt_tolerance, is_coeluted)
            )
        return epds



    def index_reformat_epds(self, list_empCpds):
        '''
        Format list of empCpds to a dictionary with fields as defined in metDataModel.
        Be careful with format inconsistency.
        '''
        new = {}
        for E in list_empCpds:
            features = []
            for peak in E['list_peaks']:
                self.peak_dict[peak[0]]['ion_relation'] = peak[1]
                features.append( self.peak_dict[peak[0]] )
            new[E['id']] = {
                'interim_id': E['id'], 
                'neutral_formula_mass': None,
                'neutral_formula': None,
                'Database_referred': [],
                'identity': [],
                'MS1_pseudo_Spectra': features,
                'MS2_Spectra': [],
                }

        return new

