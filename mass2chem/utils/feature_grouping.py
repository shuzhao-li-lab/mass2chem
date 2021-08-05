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
    [),
    ('C3H10N2_74.084398',
    [['HMDB0000002', '1,3-Diaminopropane'],
    ['HMDB0013136', '1,2 Diaminopropane']]),
    ('C4H6O3_102.031694',
    [['HMDB0000005', '2-Ketobutyric acid'],
    ['HMDB0000060', 'Acetoacetic acid']])]

Indexed_features are used internally for search speed.

SL 2021-08-04
"""





# -----------------------------------------------------------------------------
#
# not used now 
#
# -----------------------------------------------------------------------------



        
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


