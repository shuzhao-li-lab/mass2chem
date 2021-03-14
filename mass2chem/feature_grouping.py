"""
Grouping features into EmpiricalCompounds.

This is ideally done by annotation methods combined with pre-processing.
E.g. similar peak shapes can be taken into account.

Here is a fallback option. When only a feature table is supplied, 

"""

from metDataModel.core import Feature, EmpiricalCompound

# from metDataModel.derived import userData, metabolicModel

from .io import *



MASS_RANGE = (50, 2000)

# fraction of total retention time, or of ranks of retention time
# used to determine coelution of ions ad hoc
RETENTION_TIME_TOLERANCE_FRAC = 0.02    

# This will be cooridnated from common_mass?
wanted_mz_delta_pos = {
    'H+': 1.007276, 'C13': 1.0034, 'Na-H': 21.9820, '-H2O': 18.0106,
}




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


def text_to_empCpds(text):
    '''

    '''
    list_of_features = textfile_to_features()
    empCpds = group_features(list_of_features, algorithm='mummichog')
    return empCpds

def index_features(list_of_features):
    '''
    Index list_of_features to speed up search.

    Input
    -----
    List of metDataModel Feature instances.

    Return
    -------
    Dictionary of list of Features, indexed to integer of m/z.

    e.g. indexed_features[233] = [Feature3, Feature34, Feature556],
    where Feature3, Feature34, Feature556 have 233 < m/z < 234.
    '''
    indexed_features = {}
    for F in list_of_features:
        ii = int(F.mz)
        if ii in indexed_features:
            indexed_features[ii].append(F)
        else:
            indexed_features[ii] = [F]

    return indexed_features


def look_mass_relatives(
    feature0, 
    indexed_features, 
    mz_tolerance, 
    rtime_tolerance, 
    relationships):
    '''
    Look up features in indexed_features that have defined mass differences to feature0,
    and meet constraints of retention time.

    Input
    -----
    feature0: 
        feature of interest, metDataModel Feature instance.
    indexed_features: 
        Dictionary of list of Features, indexed to integer of m/z.
    relationships: 
        defined mass differences of interest, e.g. {'H+': 1.007276, 'C13': 1.0034, 'Na-H': 21.9820,}
    mz_tolerance: 
        m/z tolerance in matching.
    rtime_tolerance: 
        retention time tolerance in matching.

    Return
    -------
    List of matched Features, e.g. [('H+', Feature), ...]

    '''
    matched = []
    for r in relationships:
        mass_to_look = feature0.mz + relationships[r]
        # only need to look [floor-1, floor, floor+1]
        floor = int(mass_to_look)
        for F in indexed_features[floor-1] + indexed_features[floor] + indexed_features[floor+1]:
            if __is_mz_match__(F, feature0):
                if __is_coelution__(F, feature0):
                    matched.append( (r,F) )

    return matched
    


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



