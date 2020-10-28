'''
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
  
    
c18_annotate(new_file3)


