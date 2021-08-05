'''

use case:

input feature table
output:
    feature table with flags of contaminants, empCpd and corresponding ions
    nonredundant table with one feature (selected by highest intensity) per empCpd, and empCpd annotations




This needs to take step wise annotation.

Input data: list of features


# to import library of authentic chemical standards


# from metDataModel.derived import userData, metabolicModel
# from .io import *

# core data structures
# from metDataModel import derived


# from metDataModel.core import Feature, EmpiricalCompound

"C05769": {"formula": "C36H34N4O8", "mw": 654.269, "name": "Coproporphyrin I", "adducts": {"M+2H[2+]": 328.14177646677, "M+Br81[-]": 735.1853, "M-H2O+H[1+]": 637.2656764667701, "M-C3H4O2+H[1+]": 583.25517646677, "M-HCOOH+H[1+]": 609.27087646677, "M-CO+H[1+]": 627.28127646677, "M+K[1+]": 693.23177646677, "M+Cl[-]": 689.2379, "M+Na-2H[-]": 674.23644706646, "M-CO2+H[1+]": 611.28647646677, "M+Na[1+]": 677.25827646677, "M-2H[2-]": 326.12722353323, "M+H[1+]": 655.27627646677, "M-H4O2+H[1+]": 619.25507646677, "M(C13)-H[-]": 654.26512353323, "M+HCOONa[1+]": 723.26367646677, "M(C13)+2H[2+]": 328.64347646677004, "M+HCOOK[1+]": 739.23757646677, "M+HCOO[-]": 699.266645, "M(C13)+3H[3+]": 219.43134313343666, "M-H[-]": 653.26172353323, "M+ACN-H[-]": 694.2882685332299, "M+Cl37[-]": 691.2349, "M-H2O-H[-]": 635.25112353323, "M+Br[-]": 733.1873, "M+3H[3+]": 219.09694313343667, "M+CH3COO[-]": 713.282295, "M(C13)+H[1+]": 656.2796764667701, "M[1+]": 654.269, "M-NH3+H[1+]": 638.2497764667701, "M+NaCl[1+]": 713.2348764667701, "M+H+Na[2+]": 339.13277646677, "M+H2O+H[1+]": 673.28687646677, "M-H+O[-]": 669.25663353323, "M+K-2H[-]": 690.20994706646}}, 
"C05768": {"formula": "C36H44N4O8", "mw": 660.3159, "name": "Coproporphyrinogen I", "adducts": {"M+2H[2+]": 331.16522646677004, "M+Br81[-]": 741.2322, "M-H2O+H[1+]": 643.3125764667701, "M-C3H4O2+H[1+]": 589.30207646677, "M-HCOOH+H[1+]": 615.31777646677, "M-CO+H[1+]": 633.3281764667701, "M+K[1+]": 699.2786764667701, "M+Cl[-]": 695.2848, "M+Na-2H[-]": 680.28334706646, "M-CO2+H[1+]": 617.33337646677, "M+Na[1+]": 683.30517646677, "M-2H[2-]": 329.15067353323, "M+H[1+]": 661.3231764667701, "M-H4O2+H[1+]": 625.30197646677, "M(C13)-H[-]": 660.3120235332301, "M+HCOONa[1+]": 729.31057646677, "M(C13)+2H[2+]": 331.66692646677006, "M+HCOOK[1+]": 745.28447646677, "M+HCOO[-]": 705.3135450000001, "M(C13)+3H[3+]": 221.44697646677002, "M-H[-]": 659.30862353323, "M+ACN-H[-]": 700.33516853323, "M+Cl37[-]": 697.2818000000001, "M-H2O-H[-]": 641.2980235332301, "M+Br[-]": 739.2342000000001, "M+3H[3+]": 221.11257646677004, "M+CH3COO[-]": 719.329195, "M(C13)+H[1+]": 662.3265764667701, "M[1+]": 660.3159, "M-NH3+H[1+]": 644.2966764667701, "M+NaCl[1+]": 719.2817764667701, "M+H+Na[2+]": 342.15622646677, "M+H2O+H[1+]": 679.33377646677, "M-H+O[-]": 675.30353353323, "M+K-2H[-]": 696.2568470664601}}, 

'''



def compute_frequent_peaks(mw, formula, rule_table):
    '''
    '''
    pass



frequent_degenerate_peaks_pos = {
    (mw, 'M[1+]', ''), 
         (mw + PROTON, 'M+H[1+]', ''),
         (mw/2 + PROTON, 'M+2H[2+]', ''),

         (mw/3 + PROTON, 'M+3H[3+]', ''),
         (mw +1.0034 + PROTON, 'M(C13)+H[1+]', 'C'),
         (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]', 'C'),
         
         (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]', 'C'),

         (mw + 21.9820 + PROTON, 'M+Na[1+]', ''),  

         
         (mw +1.9958 + PROTON, 'M(S34)+H[1+]', 'S'),
         (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]', 'Cl'),

    'M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', }


frequent_degenerate_peaks_neg = {
    'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]'}


# to multiply with above
common_mz_diff = {}



mass_signatures = [
    (),
    ()
]









def grouping_mummichog(list_of_features, indexed_compounds):
    '''
    Compute list of empCpds from a list of features

    G1: core pairs of ions
    primary_ions = ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]']

    G2: search other ions based on G1
    G3: organize others, mostly singletons expected
    G4: match formula, thus indexed known compounds
    

    Input
    -----
    List of Feature instances, from the same Experiment.

    indexed_compounds:
    List of reference compounds, which have pre-computed ions.

    Return
    -------
    List of EmpiricalCompounds.

    '''
    Group1 = []
    list_of_features.sort(key=F.mz)
    indexed_features = index_features(list_of_features)
    relationships = {}
    for F in list_of_features:
        found = look_mass_relatives(F, indexed_features, relationships)
        Group1.append( found )

    # Group 2
    # extend to 







def group_features(list_of_features, algorithm='mummichog'):
    '''
    Grouping features to empCpds.

    This is mummichog algorithm, not using peak shapes, 
    which are not available if input data are feature table.
    Options can be added later when better connected with pre-processing.
    Will add other algorithms (e.g. CAMERA) later.

    Input
    -----
    List of Feature instances, from the same Experiment

    Return
    -------
    List of EmpiricalCompounds.

    Example
    -------
    

    '''
    pass



