from collections import namedtuple
import numpy as np

# from metDataModel.core import Feature

# intensities is a list; mass_id links to MassTrace
tFeature = namedtuple('tFeature', ['feature_id', 'mass_id', 'mz', 'rtime', 'rt_min', 'rt_max', 
                                'peak_quality_max', 'peak_quality_median', 'number_peaks', 'perc_peaks',
                                'selectivity_combined', 'selectivity_mz', 'intensity_mean',
                                'intensities'])


def read_feature_tuples(feature_table, 
                        id_col=True, mz_col=1, rtime_col=2, 
                        intensity_cols=(3,4), delimiter="\t"):
    '''
    Read a text feature table into a list of features.

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
    List of namedtuple 'Features': [(Feature)...], 
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



def index_features(list_of_features, range_low=20, range_high=2000):
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


def read_features(feature_table, 
                        id_col=True, mz_col=1, rtime_col=2, 
                        intensity_cols=(3,4), delimiter="\t"):
    '''
    Read a text feature table into a list of features.
    Internal use.

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



def text_to_mcg_features(textValue, MASS_RANGE = (50, 2500), delimiter='\t'):
    '''
    Parse text input to List of features, in JSON peak formats. 
    The older versions used instances of metDataModel.Feature, which is better for complex operations.
    Here, simplicity is preferred.

    Input
    -----
    Text string of feature table, mummichog format, Column order is hard coded for now, 
    1st row as header = [mz, retention_time, p_value, statistic, CompoundID_from_user].
    Each row is a unique feature. Redundant entries should be check prior to this function.

    Return
    -------
    peak_list: [{''id_number': 555, 'mz': 133.09702315984987, 'apex': 654, 
    'height': 14388.0, 
    'left_base': 648, 'right_base': 655, },
                 ...]
    list of excluded lines in input text
    '''
    list_of_features, excluded_list = [], []
    lines = textValue.splitlines()
    _header_fields = lines[0].rstrip().split(delimiter)

    for ii in range(len( lines )-1):
        y = lines[ii+1].split('\t')
        CompoundID_from_user = ''
        if len(y) > 4: CompoundID_from_user = y[4]
        [mz, retention_time, p_value, statistic] = [float(x) for x in y[:4]]
        
        if MASS_RANGE[0] < mz < MASS_RANGE[1]:
            # row # human-friendly, numbering from 1
            F = {
                'id_number': 'row'+str(ii+1),
                'mz': mz,
                'retention_time': retention_time,
                'apex': retention_time,
                'statistic_score': statistic,
                'p_value': p_value,
                'CompoundID_from_user': CompoundID_from_user,
            }
            list_of_features.append( F )
        else:
            excluded_list.append( lines[ii+1] )

    return (list_of_features, excluded_list)


def textfile_to_mcg_features(infile):
    return text_to_mcg_features(open(infile).read())

