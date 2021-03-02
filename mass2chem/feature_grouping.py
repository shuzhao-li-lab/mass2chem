"""
Grouping features into EmpiricalCompounds.

This is ideally done by annotation methods combined with pre-processing.
E.g. similar peak shapes can be taken into account.

Here is a fallback option. When only a feature table is supplied, 

"""

from metDataModel.core import Feature, EmpiricalCompound
from metDataModel.derived import userData, metabolicModel

def group_features(list_of_features):
    '''
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
