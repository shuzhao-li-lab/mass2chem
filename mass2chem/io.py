from metDataModel.core import Feature

MASS_RANGE = (50, 2500)


def text_to_features(textValue, delimiter='\t'):
    '''
    Parse text input to List of Feature instances.

    Input
    -----
    Text string of feature table, Column order is hard coded for now, 
    1st row as header = [mz, retention_time, p_value, statistic, CompoundID_from_user].
    Each row is a unique feature. Redundant entries should be check prior to this function.

    Return
    -------
    (List of metDataModel.Feature, list of excluded lines in input text)
    '''
    list_of_features, excluded_list = [], []
    lines = textValue.splitlines()
    header_fields = lines[0].rstrip().split(delimiter)

    for ii in range(len( lines )-1):
        y = lines[ii+1].split('\t')
        CompoundID_from_user = ''
        if len(y) > 4: CompoundID_from_user = y[4]
        [mz, retention_time, p_value, statistic] = [float(x) for x in y[:4]]
        
        if MASS_RANGE[0] < mz < MASS_RANGE[1]:
            # row # human-friendly, numbering from 1
            F = Feature(id='row'+str(ii+1))
            [F.mz, F.rtime] = [mz, retention_time]
            F.statistics['statistic_score'] = statistic
            F.statistics['p_value'] = p_value
            list_of_features.append( F )
        else:
            excluded_list.append( lines[ii+1] )

    return (list_of_features, excluded_list)


def textfile_to_features(infile):
    return text_to_features(open(inputFile).read())

