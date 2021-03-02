'''
This needs to take step wise annotation.

Input data: list of features





PPM_tolerance = 0.000010
RTime_tolerance = 150       # seconds in retention time, usually a small number
                            # more lenient for diff instruments
                            # and possible diff void volume

def match1(a, b):
     # m/z match only
     if abs(a-b)/a < PPM_tolerance:
         return True
     else:
         return False

PROTON = 1.00727646677
SODIUM = 21.9820 + PROTON
H2O = 18.0106
adducts = []

for x in wanted:
     #pos
     #adducts.append( list(x) + ['M+H[1+]', x[2]+PROTON] ) 
     #adducts.append( list(x) + ['M+Na[1+]', x[2]+SODIUM] )
     #neg
     adducts.append( list(x) + ['M-H[1-]', x[2]-PROTON] ) 
     adducts.append( list(x) + ['M-H2O-H[-]', x[2]-PROTON-H2O] )
     adducts.append( list(x) + ['M-2H[2-]', x[2]/2-PROTON] )
        
# w is tmp; not confused when code blocks are run in wrong orders
w = open(new_file3).readlines()
s = 'Name\tID\tmono_mass\tion\ttheretical_mz\t' + w[0]

for line in w[1:]:
     mz = float(line.split('\t')[0])
     for a in adducts:
         if match1(mz, a[4]):
             s += '\t'.join([str(x) for x in a]) + '\t' + line
             
with open('microbial_matched_201705_C18neg.txt', 'w') as file:
    file.write(s)



confirmed = 'C18Neg_Orbitrap_IROA_Ken20170623.txt'

new = []
w = open(confirmed).read()
for line in w.splitlines()[1:636]:
    a = line.split('\t')
    # CNAME, M-H, Retention Time (C18), PARENT_CID, Preferred Adduct, m/z
    d = [a[ii] for ii in (3, 7, 8, 9, 19, 20)]
    if d[2].strip() != "*":
        # if alternative ion is preferred, it precedes default M-H
        # data format (Name, ID, ion, m/z, rtime)
        if d[4].strip():
            new.append( [d[0], d[3], d[4], d[5], d[2]] )
        else:
            new.append( [d[0], d[3], 'M-H', d[1], d[2]] )

confirmed = [[y.replace('*', '') for y in x] for x in new]
print("Got %d confirmed metabolites in library." %len(confirmed))

'''
# core data structures

from metDataModel import derived


# to import library of authentic chemical standards


#print(confirmed)

def match2(F1, F2):
    if abs(F1[0]-F2[0])/F1[0] < PPM_tolerance and abs(F1[1] - F2[1]) < RTime_tolerance:
        return True
    else:
        return False

def c18_annotate(infile, confirmed = confirmed):
    '''
    Last update, 2017-06-23
    '''
    STARTROW = 1
    col1, col2 = 0, 1   #assumming first two cols as m/z, rtime
    new = open(infile).readlines()
    s = 'Name\tID\tion\tlibrary_mz\tlib_retention_time\t' + new[0]
    
    for line in new[STARTROW: ]:
        x = line.split('\t')
        for y in confirmed:
            try:
                if match2( (float(x[col1]), float(x[col2])), (float(y[3]), float(y[4])) ):
                    s += '\t'.join(y) + '\t' + line
            except ValueError:
                pass    #ignoring irregular fields for now
                #print(y)
    
    with open('c18_annotate_' + infile, 'w') as file:
        file.write( s )
  
# standalone     
# c18_annotate(new_file3)



# metabolicNetwork

class DataMeetModel:
    '''
    # Key change to make in v3
    Move "annotation", i.e., listing empCpd, to separate package, mass2chem
    
    
    returns the tracking map btw massFeatures - EmpiricalCompounds - Compounds.

    Number of EmpiricalCompounds will be used to compute pathway enrichment, and control for module analysis.
    New in v2, but will move to a boutique database in v3.
    
    N:N matches:
    when a Compound matched to multiple MassFeatures, split by retention time to EmpiricalCompounds;
    when a Mass Feature matched to multiple Compounds, no need to do anything.
    
    ??Default primary ion is enforced, so that for an EmpiricalCompound, primary ion needs to exist before other ions.



    Key output:
    empCpd2Features = {empCpd: (), ...,}
    empCpd2Cpds = {empCpd: (), ...,}




    '''
    def __init__(self, metabolicModel, userData):
        '''
        # from ver 1 to ver 2, major change in .match()
        Trio structure of mapping
        (M.row_number, EmpiricalCompounds, Cpd)
        
        '''
        self.model = metabolicModel
        self.data = userData
        
        # retention time window for grouping, based on fraction of time or ranks
        self.rtime_tolerance = self.data.max_retention_time * RETENTION_TIME_TOLERANCE_FRAC
        self.rtime_tolerance_rank = len(self.data.ListOfMassFeatures) * RETENTION_TIME_TOLERANCE_FRAC
        
        # major data structures
        # web
        if self.data.web:
            wanted_ions = self.data.paradict['wanted_adduct_list']
        # local
        else:
            wanted_ions = wanted_adduct_list[ self.data.paradict['mode'] ]

        self.IonCpdTree = self.__build_cpdindex__(wanted_ions)

        self.rowDict = self.__build_rowindex__( self.data.ListOfMassFeatures )
        self.ListOfEmpiricalCompounds = self.get_ListOfEmpiricalCompounds()
        
        # this is the reference list
        self.mzrows = [M.row_number for M in self.data.ListOfMassFeatures]
        
        self.rowindex_to_EmpiricalCompounds = self.__make_rowindex_to_EmpiricalCompounds__()
        self.Compounds_to_EmpiricalCompounds = self.__index_Compounds_to_EmpiricalCompounds__()
        
        # this is the sig list
        self.significant_features = self.data.input_featurelist
        self.TrioList = self.batch_rowindex_EmpCpd_Cpd( self.significant_features )

    def __build_cpdindex__(self, wanted_ions):
        '''
        indexed Compound list, to speed up m/z matching.
        Limited to MASS_RANGE (default 50 ~ 2000 dalton).
        
        changing from adduct_function to wanted_adduct_list dictionary
        
        wanted_adduct_list['pos_default'] = ['M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 
                    'M+Na[1+]', 'M+H+Na[2+]', 'M+HCOONa[1+]'
                    ],
        
        # 
        >>> metabolicModels['human_model_mfn']['Compounds'].items()[92]
        ('C00217', {'formula': '', 'mw': 147.0532, 'name': 'D-Glutamate; D-Glutamic acid; D-Glutaminic acid; D-2-Aminoglutaric acid',
         'adducts': {'M+2H[2+]': 74.53387646677, 'M+Br81[-]': 227.9695, 'M-H2O+H[1+]': 130.04987646677, 
         'M-C3H4O2+H[1+]': 76.03937646677, 'M-HCOOH+H[1+]': 102.05507646676999, 'M-HCOONa+H[1+]': 80.07307646677, 
         'M+K[1+]': 186.01597646677, 'M+Cl[-]': 182.0221, 'M+Na-2H[-]': 167.02064706646001, 'M-CO2+H[1+]': 104.07067646677, 
         'M+Na[1+]': 170.04247646677, 'M+Br[-]': 225.9715, 'M(S34)-H[-]': 148.04172353323, 'M+H[1+]': 148.06047646677, 
         'M-H4O2+H[1+]': 112.03927646677, 'M(C13)-H[-]': 147.04932353323, 'M(Cl37)-H[-]': 148.04312353323, 'M+HCOONa[1+]': 216.04787646677, 'M(C13)+2H[2+]': 75.03557646677, 'M+HCOOK[1+]': 232.02177646677, 'M-CO+H[1+]': 120.06547646677, 'M+HCOO[-]': 192.050845, 'M(C13)+3H[3+]': 50.359409800103336, 'M(Cl37)+H[1+]': 150.05767646677, 'M-H[-]': 146.04592353323, 'M+ACN-H[-]': 187.07246853323, 'M+Cl37[-]': 184.0191, 'M-H2O-H[-]': 128.03532353322998, 'M(S34)+H[1+]': 150.05627646677002, 'M-HCOOK+H[1+]': 64.09917646677, 'M+3H[3+]': 50.025009800103334, 'M+CH3COO[-]': 206.066495, 'M(C13)+H[1+]': 149.06387646677, 'M[1+]': 147.0532, 'M-NH3+H[1+]': 131.03397646677, 'M+NaCl[1+]': 206.01907646677, 'M+H+Na[2+]': 85.52487646677, 'M+H2O+H[1+]': 166.07107646677002, 'M-H+O[-]': 162.04083353323, 'M+K-2H[-]': 182.99414706646002, 'M-2H[2-]': 72.51932353323001}})
        >>> len(metabolicModels['human_model_mfn']['Compounds'])
        3560
        '''

        IonCpdTree = []
        
        for ii in range(MASS_RANGE[1]+1): 
            IonCpdTree.append([])       #empty lists for anything below MASS_RANGE
            
        # iteritems vs items is contention of efficiency, but there's change btw Python 2 and Python 3...
        for c,d in self.model.Compounds.items():
            if d['mw']:                 #sanity check; bypass mistake in adducts type
                for ion,mass in d['adducts'].items():
                    if ion in wanted_ions and MASS_RANGE[0] < mass < MASS_RANGE[1]:
                        IonCpdTree[ int(mass) ].append( (c, ion, mass) )
                
        # tree: (compoundID, ion, mass), ion=match form; mass is theoretical
        return IonCpdTree


    def __build_rowindex__(self, ListOfMassFeatures):
        '''
        Index list of MassFeatures by row# in input data
        '''
        rowDict = {}
        for M in ListOfMassFeatures: rowDict[M.row_number] = M
        return rowDict


    def __match_all_to_all__(self):
        '''
        Major change of data structure here in version 2.
        In ver 1, matched m/z is stored in each Compound instance.
        Here, we produce mapping dictionaries for
            * mzFeatures to theoretical ions
            * Compounds to mzFeatures
        Then, 
            * EmpiricalCompounds are determined within Compound matched mzFeatures, considering retention time.
        
        
        '''
        self.__match_to_mzFeatures__()
        self.cpd2mzFeatures = self.index_Compounds_to_mzFeatures()
        return self.compound_to_EmpiricalCompounds()
        

    def __match_to_mzFeatures__(self):
        '''
        Fill mzFeatures with matched ions and compounds
        '''
        for M in self.data.ListOfMassFeatures:
            M.matched_Ions = self.__match_mz_ion__(M.mz, self.IonCpdTree)
        
        
    def index_Compounds_to_mzFeatures(self):
        '''
        compound ID - mzFeatures
        run after self.__match_to_mzFeatures__()
        L: (compoundID, ion, mass)
        cpd2mzFeatures[compoundID] = [(ion, mass, mzFeature), ...]
        '''
        cpd2mzFeatures = {}
        for M in self.data.ListOfMassFeatures:
            for L in M.matched_Ions:
                if L[0] in cpd2mzFeatures:
                    cpd2mzFeatures[L[0]].append( (L[1], L[2], M) )
                else:
                    cpd2mzFeatures[L[0]] = [(L[1], L[2], M)]
        
        print ("Got %d cpd2mzFeatures" %len(cpd2mzFeatures))
        return cpd2mzFeatures
        
        
    def __match_mz_ion__(self, mz, IonCpdTree):
        '''
        L: (compoundID, ion, mass)
        return ions matched to m/z
        '''
        floor = int(mz)
        matched = []
        mztol = mz_tolerance(mz, self.data.paradict['instrument'])
        for ii in [floor-1, floor, floor+1]:
            for L in IonCpdTree[ii]:
                if abs(L[2]-mz) < mztol:
                    matched.append( L )
                    
        return matched

    def compound_to_EmpiricalCompounds(self):
        '''
        EmpiricalCompounds are constructed in this function.
        First splitting features matching to same compound by retention time;
        then merging those matched to same m/z features.
        run after self.index_Compounds_to_mzFeatures()
        '''
        ListOfEmpiricalCompounds = []
        for k,v in self.cpd2mzFeatures.items():
            ListOfEmpiricalCompounds += self.__split_Compound__(k, v)      # getting inital instances of EmpiricalCompound
            
        print ("Got %d ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        
        # merge compounds that are not distinguished by analytical platform, e.g. isobaric
        return self.__merge_EmpiricalCompounds__( ListOfEmpiricalCompounds )
        
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

    def __split_Compound__(self, compoundID, list_match_mzFeatures):
        '''
        Determine EmpiricalCompounds among the ions matched to a Compound;
        return list of EmpiricalCompounds (not final, but initiated here).
        
        The retention time is grouped by tolerance value; 
        This method should be updated in the future.
        
        input data format:
        cpd2mzFeatures[compoundID] = list_match_mzFeatures = [(ion, mass, mzFeature), ...]
        
        '''
        # unpacked format: [retention_time, row_number, ion, mass, compoundID]
        all_mzFeatures = [(L[2].retention_time, L[2].row_number, L[0], L[1], compoundID) for L in list_match_mzFeatures]
        all_mzFeatures.sort()
        ECompounds = []
        tmp = [ all_mzFeatures[0] ]
        for ii in range(len(all_mzFeatures)-1):

            if self.__is_coelution__( self.rowDict[all_mzFeatures[ii+1][1]], self.rowDict[all_mzFeatures[ii][1]] ):
                tmp.append(
                            all_mzFeatures[ii+1] )
            else:
                ECompounds.append( EmpiricalCompound( tmp ) )
                tmp = [ all_mzFeatures[ii+1] ]
        
        ECompounds.append( EmpiricalCompound( tmp ) )
        return ECompounds


    def __merge_EmpiricalCompounds__(self, ListOfEmpiricalCompounds):
        '''
        If ion/mzFeatures are the same, merge EmpiricalCompounds
        EmpiricalCompounds.join() adds Compounds
        
        Because EmpiricalCompounds.str_row_ion uses mzFeatures sorted by row_number, this is 
        '''
        mydict = {}
        for L in ListOfEmpiricalCompounds:
            if L.str_row_ion in mydict:
                mydict[ L.str_row_ion ].join(L)
            else:
                mydict[ L.str_row_ion ]= L
        
        print ("Got %d merged ListOfEmpiricalCompounds" %len(mydict))
        return mydict.values()

    def __make_rowindex_to_EmpiricalCompounds__(self):
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.massfeature_rows:
                if m in mydict:
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict

    def __index_Compounds_to_EmpiricalCompounds__(self):
        '''
        Make dict cpd - EmpiricalCompounds
        '''
        mydict = {}
        for E in self.ListOfEmpiricalCompounds:
            for m in E.compounds:
                if m in mydict:
                    mydict[m].append(E)
                else:
                    mydict[m] = [E]
                    
        return mydict
        

    def batch_rowindex_EmpCpd_Cpd(self, list_features):
        '''
        Batch matching from row feature to Ecpds; Use trio data structure, (M.row_number, EmpiricalCompounds, Cpd).
        Will be used to map for both sig list and permutation lists.
        '''
        new = []
        for f in list_features:
            for E in self.rowindex_to_EmpiricalCompounds.get(f, []):
                for cpd in E.compounds:
                    new.append((f, E, cpd))
            
        return new

            
    def get_ListOfEmpiricalCompounds(self):
        '''
        Collect EmpiricalCompounds.
        Initiate EmpCpd attributes here.
        '''
        ListOfEmpiricalCompounds, ii = [], 1
        for EmpCpd in self.__match_all_to_all__():
            EmpCpd.evaluate()
            EmpCpd.EID = 'E' + str(ii)
            EmpCpd.get_mzFeature_of_highest_statistic( self.rowDict )
            ii += 1
            if self.data.paradict['force_primary_ion']:
                if EmpCpd.primary_ion_present:
                    ListOfEmpiricalCompounds.append(EmpCpd)
            else:
                ListOfEmpiricalCompounds.append(EmpCpd)
        
        print ("Got %d final ListOfEmpiricalCompounds" %len(ListOfEmpiricalCompounds))
        return ListOfEmpiricalCompounds

    
    def to_json(self):
        '''
        JSON export to be consumed by downstream functions

        empCpd2Cpds = {empCpd: (), ...,}

        Will update later in accordance to 
        https://github.com/shuzhao-li/metDataModel

        '''

        empCpd2Features, empCpd2Cpds = {}, {}
        for E in self.ListOfEmpiricalCompounds:
            empCpd2Features[E.EID] = E.massfeature_rows
            empCpd2Cpds[E.EID] = E.compounds

        return {
            'metabolic_model': self.model.version,
            'empCpd2Features': empCpd2Features,
            'empCpd2Cpds': empCpd2Cpds,
        }
