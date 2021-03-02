from metDataModel.core import Feature

def textfile_to_features(infile):
    return read_text_to_features(open(inputFile).read())

def read_text_to_features(textValue, delimiter='\t'):
    '''
    Column order is hard coded for now, as mz, retention_time, p_value, statistic, CompoundID_from_user
    '''
    #
    lines = __check_redundant__( textValue.splitlines() )
    self.header_fields = lines[0].rstrip().split(delimiter)
    excluded_list = []
    for ii in range(len( lines )-1):
        y = lines[ii+1].split('\t')
        
        CompoundID_from_user = ''
        if len(y) > 4: CompoundID_from_user = y[4]
        [mz, retention_time, p_value, statistic] = [float(x) for x in y[:4]]
        
        # row_number, mz, retention_time, p_value, statistic, CompoundID_from_user
        if MASS_RANGE[0] < mz < MASS_RANGE[1]:
            # row # human-friendly, numbering from 1
            self.ListOfMassFeatures.append( 
                MassFeature('row'+str(ii+1), mz, retention_time, p_value, statistic, CompoundID_from_user) 
                )
        else:
            excluded_list.append( (ii, mz, retention_time) )
    
    if excluded_list:
        print( "Excluding %d features out of m/z range %s." %(len(excluded_list), str(MASS_RANGE)) )

def __check_redundant__(L):
    redundant = len(L) - len(set(L))
    if redundant > 0:
        print( "Your input file contains %d redundant features." %(redundant) )
    return L
