'''
This script implements a dynamic programming approach for de-novo L4 annotation of m/z values.

Given an m/z value and list of possible isotopes, it performs a top down recursive search implemented
using dynamic programming. The solutions are stored in a table that exists outside function scope for 
reuse in later problems. 

The output is the set of "components" or formula elements that consist of the solution. 

The approach is similar to, but distinct from other tools. 

First, there is no database of any kind. 
Second, the search space is much more flexible. We can have limits of elements and element
    combinations using a simple rule mechanism.
Third, with some computational improvements, we may be able to search of all chemical space
    for many sample types. 

The small size of the dynamic programming solution has implications if certain criteria are met. For 
newer AMD processors with very large L3 caches, the entirety of the backprop table can fit in CPU 
cache, possibly making this code extremely fast despite not having the best big O complexity. Second
re-writing this code in C++


'''

import numpy as np
from collections import defaultdict, Counter


class BackpropagateDP:
    def __init__(self, components, elements, max_mass, min_NAP, round, extra_rules=None):
        """
        Initialize the backpropagation-based evaluator.

        Parameters:
        - components: List of components, each being a dictionary with 'mz_delta' and other properties.
        - max_solutions: Maximum number of steps to consider for solutions.
        """
        if extra_rules is None:
            extra_rules = {}
        self.components = np.array(components)
        self.elements = elements
        # split components into groups by element that are mutually exclusive, consider only one set per step. 
        # put a null component in each set to ensure that we can exclude them though. 

        self.component_sets = {}
        for i, comp in enumerate(components):
            if comp['name'] not in self.component_sets:
                self.component_sets[comp['name']] = set()
            self.component_sets[comp['name']].add(i)

        self.max_solutions = len(self.elements)
        self.dp_cache = [{} for i in range(self.max_solutions + 1)]  # Caches intermediate results
        self.valences = np.array([x['valence'] for x in self.components], dtype=np.int16)
        self.naps = np.array([x['nap'] for x in self.components], dtype=np.float16)
        self.max_mass = max_mass
        self.min_NAP = min_NAP
        self.round = round
        self.extra_rules = extra_rules


    def backpropagate(self, query_mass, step=None):
        # first round the query mass, this makes it have finite number of states. 
        # memory consumption is therefore now ~ O(max_mz * 10**ROUND)
        if step is None:
            step = len(self.elements)
        query_mass = round(query_mass, self.round)

        # if we have the value, return it
        if query_mass in self.dp_cache:
            return self.dp_cache[step][query_mass]

        # out of bounds, too small or too big
        if query_mass < -1*(10**-self.round) or query_mass > self.max_mass:
            self.dp_cache[step][query_mass] = []
            return self.dp_cache[step][query_mass]

        # lets define base cases and what needs to be returned
        if abs(query_mass) <= 1*(10**-self.round):
            # here the query mass is within the target value, return an empty initial solution
            self.dp_cache[step][query_mass] = [frozenset()]
            return self.dp_cache[step][query_mass]
        
        if step == 0:
            return []
            #self.dp_cache[query_mass] = []
            #return self.dp_cache[query_mass]

        # now, we didn't have it cached and it wasn't zero and we still have table left to traverse
        
        self.dp_cache[step][query_mass] = set()
        for comp_i, component in [(i,x) for (i,x) in enumerate(self.components) if i in self.component_sets[self.elements[step - 1]]]:
            used = set()

            # solution is too massive
            if query_mass < component['mz_delta'] - (10**-self.round):
                continue
            
            MIN_NAP_eff = self.min_NAP / component['nap']
            for previous_path in self.backpropagate(query_mass - component['mz_delta'], step - 1):
                cpath = self.components[list(previous_path)]
                if (frozenset(previous_path) in used) \
                    or (component['name'] in [c['name'] for c in cpath]) \
                    or (np.prod(self.naps[list(previous_path)]) < MIN_NAP_eff) \
                    or (np.sum(self.valences[list(previous_path)]) + component['valence'] + 2 < 0)\
                    or any(Counter(tag for c in cpath for tag in c["tags"])[tag] > self.extra_rules.get(tag, float('inf')) for tag in self.extra_rules):
                    continue
                used.add(frozenset(previous_path))
                self.dp_cache[step][query_mass].add(frozenset(previous_path | {comp_i}))
        return self.dp_cache[step][query_mass]

    def query(self, query_mass):
        """
        Query paths leading to a specific mass with lazy evaluation.

        Parameters:
        - query_mass: The target mass to compute paths for.
        - max_steps: Maximum steps allowed (default: self.max_solutions).

        Returns:
        - List of paths leading to the query_mass.
        """
        return self.backpropagate(query_mass)
    
class MassExplainer:
    def __init__(self, csv_file, elements=["12C", "16O", "14N", "35Cl", "1H"], max_mass=1000, min_NAP=0.4, round=3):
        """
        Initialize the MassExplainer with isotopologue data.

        Args:
            csv_file (str): Path to the CSV file containing isotopologue data.
            max_solutions (int): Maximum number of components in a solution.
        """
        self.max_mass = max_mass
        self.min_NAP = min_NAP
        self.round = round
        self.components = self._parse_csv(csv_file, elements)
        self.lazy_dp = BackpropagateDP(self.components, elements, max_mass, min_NAP, round) #, max_solutions=len(ELEMENTS))

    def _parse_csv(self, file_path, to_extract):
        """
        Parse the CSV file and calculate isotope mass deltas relative to the most abundant isotopologue.

        Args:
            file_path (str): Path to the CSV file.

        Returns:
            list: List of components with relative mass differences.
        """
        isotopes = {}
        # this is hacky - fix this!!!
        for ele in to_extract:
            with open(file_path, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    name, mz_delta, prob, formula, valence, max_count,tags, _ = line.strip().split(',')
                    if name == ele:
                        mz_delta = float(mz_delta)
                        prob = float(prob)
                        valence = int(valence)
                        max_count = int(max_count)
                        tags = tags.split("|") if tags else []
                        if formula not in isotopes:
                            isotopes[formula] = []
                        isotopes[formula].append({'name': name, 'mz_delta': mz_delta, 'prob': prob, 'valence': valence, 'max': max_count, 'tags': tags})
        # Normalize isotopes to the most abundant (highest probability)
        components = []
        for formula, isotope_list in isotopes.items():
            isotope_list.sort(key=lambda x: -x['prob'])  # Sort by descending probability
            isotope = isotope_list[0]
            components.append({
                        'name': isotope['name'],
                        'name2': isotope['name'] + "_0",
                        'mz_delta': 0,
                        'nap': 1,
                        'formula': {},
                        'tags': [],
                        'valence': 0
                    })
            for num_iso in range(1, isotope['max']):
                if isotope['prob'] ** num_iso > self.min_NAP and float(isotope['mz_delta']) * num_iso < self.max_mass:
                    components.append({
                        'name': isotope['name'],
                        'name2': isotope['name'] + "_" + str(num_iso),
                        'mz_delta': float(isotope['mz_delta']) * num_iso,
                        'nap': isotope['prob'] ** num_iso,
                        'formula': {formula: num_iso},
                        'tags': isotope['tags'] * num_iso if isotope['tags'] else [],
                        'valence': isotope['valence'] * num_iso})
        return components

    def explains(self, query_mass):
        """
        Determine if the mz_value can be explained by a linear combination of components.

        Args:
            mz_value (float): Observed m/z value.
            ppm_tolerance (float): Error tolerance in parts-per-million.

        Returns:
            list: Solutions that explain the mz_value within the given tolerance.
        """
        from mass2chem.formula import dict_to_hill_formula
        results = set()
        for path in self.lazy_dp.query(query_mass):
            path_formula = {}
            path_nap = 1.0
            path_mass = 0.0
            for comp_i in path:
                path_formula.update(self.components[comp_i]['formula'])
                path_nap *= self.components[comp_i]['nap']
                path_mass += self.components[comp_i]['mz_delta']
            if path_formula:
                results.add((dict_to_hill_formula(path_formula), path_nap, path_mass))
        return {'query_mass': query_mass, 
                'possible_m0': [{
                'formula': r[0],
                'nap': r[1],
                'mass': r[2]
            } for r in results]}




# Example Usage
if __name__ == "__main__":
    import random
    import time
    X, Y = [], []
    # Initialize the MassExplainer with the provided CSV file
    explainer = MassExplainer("./isotopes.csv")
    for i in range(0, 1000):
        q = i 
        t1 = time.time()
        c = explainer.explains(q)
        t2 = time.time()

        if c['possible_m0']:
            print(q)
            print(c)
            print(t2 - t1)
        X.append(q)
        Y.append(t2-t1)

    import matplotlib.pyplot as plt
    plt.scatter(X, Y)
    plt.show()
