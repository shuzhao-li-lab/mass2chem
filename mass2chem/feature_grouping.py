"""
Grouping features into EmpiricalCompounds.

This is ideally done by annotation methods combined with pre-processing.
E.g. similar peak shapes can be taken into account.

Here is a fallback option. When only a feature table is supplied, 

"""

from metDataModel.core import Feature, EmpiricalCompound

# from metDataModel.derived import userData, metabolicModel

from .io import *

def group_features(list_of_features, algorithm='mummichog'):
    '''
    Input
    -----
    List of Feature instances, from the same Experiment

    index features to integers, to save search time

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


def look_mass_relatives(feature1, indexed_features, relationships):
    '''
    Look up features in indexed_features that meet requirement of relationships to feature1.

    Input
    -----


    Return
    -------
    List of .

    '''
    pass


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



