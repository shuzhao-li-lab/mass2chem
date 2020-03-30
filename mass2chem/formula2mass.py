'''
# using pychemy module. It failed to install using pip because of `Open Babel` binding.
# But the mass calcuation function is not depdendent on `Open Babel`, 
# thus pychemy folder is copied here and it's usable as below

from pychemy.molmass import Formula

massdict_cpds2include = {}
unassigned = []
for x in cpds2include:
    if model.dict_cpds_mass.has_key(x):
        massdict_cpds2include[x] = model.dict_cpds_mass[x]
    elif dict_formula.has_key(x):
        try:
            massdict_cpds2include[x] = Formula( dict_formula[x] ).isotope.mass
        except:
            print(dict_formula[x])
            unassigned.append(x)
    else:
        unassigned.append(x)

print(len(massdict_cpds2include))
print("unassigned, ", len(unassigned))
'''

from pychemy.molmass import Formula

def formula2mass( x):
    return Formula( x ).isotope.mass

