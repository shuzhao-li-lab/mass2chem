import json
from functools import cache
import numpy as np
from scipy.stats import multinomial
from mass2chem.formula import parse_chemformula_dict
from itertools import product
import heapq

@cache
# this loads the isotope data file, need to consolidate to one representation
def isotopes(iso_path="/Users/mitchjo/Projects/ChemoselectiveDerivatizationTools/New/isotopes.json"):
    with open(iso_path) as isotope_fh:
        return json.load(isotope_fh)

@cache
# this builds and stores the multinomial object for use in calc_iso_nap
def build_multinomial(element, vector):
    return multinomial(sum(vector), isotopes()[element][2])

@cache
# calculate the nap of element and a given iso combination specified by vector
def calc_iso_nap(element, vector):
    return build_multinomial(element, vector).pmf(vector)

@cache
# calculate the mass of element and a given iso combination specified by vector
# offset lets us skip elements, useful for delta mass calcs
def calc_mass(element, vector, offset=0):
    return np.dot(isotopes()[element][1][offset:], vector[offset:])

@cache
# calculate the string representation of an element's isotopologues
def calc_iso(element, vector, offset=0):
    return '.'.join([iso + str(i) for iso, i in zip(isotopes()[element][0][offset:], vector[offset:]) if i > 0])

@cache
# calculate the string representation of an element's labels
def calc_lbl(element, vector):
    return calc_iso(element, vector, offset=1)

@cache
def __permute(num_isos, remaining_atoms, current_sum=0, current_perm=()):
    if num_isos == 1:
        return [current_perm + (remaining_atoms - current_sum,)]
    return [r for i in range(remaining_atoms - current_sum + 1) for r in __permute(num_isos - 1, remaining_atoms, current_sum + i, current_perm + (i,))]

def formula_to_components(formula):
    formula = parse_chemformula_dict(formula) if isinstance(formula, str) else formula
    components = {ele: [] for ele in formula}
    for ele, count in formula.items():
        for c_vec in __permute(len(isotopes()[ele][0]), count):
            components[ele].append({
                "nap": calc_iso_nap(ele, c_vec),
                "mass": calc_mass(ele, c_vec),
                "vector": c_vec
            })
        components[ele] = sorted(components[ele], key=lambda x: -x['nap'])
    return components

@cache
def generate_isotopologues(formula, min_NAP=0.01):
    formula = parse_chemformula_dict(formula)
    components = formula_to_components(formula)
    iso_heap = [[-np.prod([c[0]['nap'] for c in components.values()]), list([0 for _ in formula])]]
    heapq.heapify(iso_heap)
    in_heap = set()
    while iso_heap:
        _, top_solution = heapq.heappop(iso_heap)
        for i in range(len(top_solution)):
            new_solution = [x if j != i else x + 1 for j, x in enumerate(top_solution)]
            if new_solution[i] < len(components[list(formula.keys())[i]]) and tuple(new_solution) not in in_heap:
                naps = np.prod([components[ele][index]['nap'] for ele, index in zip(formula.keys(), new_solution)])
                if abs(naps) > min_NAP:
                    heapq.heappush(iso_heap, [-np.prod(naps), list(new_solution)])
                else:
                    pass
        if iso_heap and top_solution:
            yield compute_iso_properties([components[i][j]['vector'] for i, j in zip(formula.keys(), top_solution)], list(formula.keys()), 'unlabeled', 1)

    
def compute_iso_properties(combo_vector, element_list, label, weight=1):
    if combo_vector is not None:
        return {
            "NAP": {label: weight * np.prod([calc_iso_nap(e, ic) for e, ic in zip(element_list, combo_vector)]),},
            "MASS": np.sum([calc_mass(e, ic) for e, ic in zip(element_list, combo_vector)]),
            "ISO_LIST": ",".join([calc_iso(e, ic) for e, ic in zip(element_list, combo_vector)]),
            "combo": {label: (combo_vector, tuple(element_list))} 
        }

if __name__ == "__main__":
    for x in generate_isotopologues("C64H120O6"):
        if x:
            print(x)
