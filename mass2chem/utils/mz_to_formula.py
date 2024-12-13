import numpy as np
from collections import defaultdict

# values are rounded to this decimal place to enable states for dynamic programming.
ROUND = 4

# minNAP of solutions, assuming peak is monoisotopic and tallest assumes NAP of 0.5. 
MIN_NAP = 0.4

# limit bounds on how many of each element to expect
ELEMENT_RANGE = {
    "C": 40,
    "H": 130,
    "O": 10,
    "N": 5,
    "Cl": 6,
    "Br": 6,
}

#how many hydrogen bonds we expect per element max
VALENCE_TABLE = {
    "C": 3,
    "O": 1,
    "N": 2,
    "H": -1,
    "Cl": -1,
    "Br": -1,
}

# these can limit the space in a more complex way
TAGS = {
    "Cl": ["halogen"],
    "N": ["halogen"]
}

COUNT_RULES = {
    "halogen": 6
}

class BackpropagateDP:
    def __init__(self, components, max_solutions):
        """
        Initialize the backpropagation-based evaluator.

        Parameters:
        - components: List of components, each being a dictionary with 'mz_delta' and other properties.
        - max_solutions: Maximum number of steps to consider for solutions.
        """
        self.components = components
        self.max_solutions = max_solutions
        self.dp_cache = {}  # Caches intermediate results

    def backpropagate(self, query_mass, step, error, ROUND=ROUND, MIN_NAP=MIN_NAP):
        # first round the query mass, this makes it have finite number of states. 
        # memory consumption is therefore now ~ O(max_mz * 10**ROUND)
        query_mass = round(query_mass, ROUND)

        # if we have the value already, just return it
        if query_mass in self.dp_cache:
            # first, if value is in table, return cached value
            return self.dp_cache[query_mass]

        # lets define base cases and what needs to be returned
        if -error <= query_mass <= error:
            # here the query mass is within the target value, return an empty initial solution
            return [{"path": [], "nap": 1.0, "val": 2, "formula": [], "mass": 0, "counts": {}}]
        if step == 0:
            # here we have traversed the table but no solution was found
            return []

        # now, we didn't have it cached and it wasn't zero and we still have table left to traverse
        self.dp_cache[query_mass] = []
        for comp_i, component in enumerate(self.components):
            used = set()

            # solution is too massive
            if query_mass < component['mz_delta'] - error:
                continue
            
            # here is the dynamic programming.
            # okay, the dp_cache is a table indexed by mass, we are going to remove component mass and 
            # see if the resulting mass difference can be explained. It will be explained if there are a chain of 
            # elements that reach the base case or we run out of table. 
            MIN_NAP_eff = MIN_NAP / component['nap']
            for previous_path in self.backpropagate(query_mass - component['mz_delta'], step - 1, error, ROUND, MIN_NAP):                     
            
                pp = previous_path
                # skip if we've seen this path before
                # solution reuses table entry 
                if frozenset(sorted(previous_path["path"])) in used:
                    continue

                # skip if elements are reused
                # solution reuses element
                if component['name'] in [self.components[i]['name'] for i in pp['path']]:
                    continue
                
                # skip if NAP not high enough
                # nap is not sufficient

                #if pp and np.prod([self.components[i]['nap'] for i in pp]) < (MIN_NAP / component['nap']):
                if np.prod([self.components[i]['nap'] for i in pp['path']]) < MIN_NAP_eff:
                #if previous_path['nap'] < MIN_NAP / component['nap']:
                    continue

                # skip if valence not positive
                # too many hydrogen or equivalents
                if np.sum([self.components[i]['valence'] for i in pp['path']]) <= -1 * component['valence']:
                    continue

                #used.add(frozenset(sorted(previous_path["path"])))
                new_path = sorted(previous_path["path"] + [comp_i])
                if frozenset(new_path) not in used:
                    
                    new_counts = previous_path["counts"].copy()
                    for k in component['tags']:
                        new_counts[k] = new_counts.get(k, 0) + 1

                    self.dp_cache[query_mass].append({
                        "path": new_path,
                        #"nap": previous_path["nap"] * component['nap'],
                        #"val": previous_path["val"] + component['valence'],
                        #"formula": previous_path["formula"] + [component['name']],
                        #"mass": component['mz_delta'] + previous_path["mass"],
                        "counts": new_counts
                    })
        # self.dp_cache[query_mass] would have been populated during the traversal
        return self.dp_cache[query_mass]

    def query(self, query_mass, max_steps=None):
        """
        Query paths leading to a specific mass with lazy evaluation.

        Parameters:
        - query_mass: The target mass to compute paths for.
        - max_steps: Maximum steps allowed (default: self.max_solutions).

        Returns:
        - List of paths leading to the query_mass.
        """
        max_steps = max_steps if max_steps is None else self.max_solutions
        error = (query_mass / 1e6 * 5) + (10**-ROUND)
        return self.backpropagate(query_mass, max_steps, error)
    
class MassExplainer:
    def __init__(self, csv_file, max_solutions=10):
        """
        Initialize the MassExplainer with isotopologue data.

        Args:
            csv_file (str): Path to the CSV file containing isotopologue data.
            max_solutions (int): Maximum number of components in a solution.
        """
        self.max_solutions = max_solutions
        self.components = self._parse_csv(csv_file)
        self.lazy_dp = BackpropagateDP(self.components, max_solutions=5)

    def _parse_csv(self, file_path):
        """
        Parse the CSV file and calculate isotope mass deltas relative to the most abundant isotopologue.

        Args:
            file_path (str): Path to the CSV file.

        Returns:
            list: List of components with relative mass differences.
        """
        isotopes = {}
        with open(file_path, 'r') as f:
            next(f)  # Skip header
            for line in f:
                name, mz_delta, prob, formula = line.strip().split(',')
                mz_delta = float(mz_delta)
                prob = float(prob)
                if formula not in isotopes:
                    isotopes[formula] = []
                isotopes[formula].append({'name': name, 'mz_delta': mz_delta, 'prob': prob})

        # Normalize isotopes to the most abundant (highest probability)
        components = []
        for formula, isotope_list in isotopes.items():
            isotope_list.sort(key=lambda x: -x['prob'])  # Sort by descending probability
            for i in range(1, ELEMENT_RANGE.get(formula, 3)):
                isotope = isotope_list[0]
                if isotope['prob'] ** i > MIN_NAP:
                    components.append({
                        'name': isotope['name'],
                        'name2': isotope['name'] + "_" + str(i),
                        'mz_delta': float(isotope['mz_delta']) * i,
                        'nap': isotope['prob'] ** i,
                        'formula': {formula: i},
                        'tags': TAGS.get(formula, []) * i,
                        'valence': VALENCE_TABLE.get(formula, -1) * i
                    })
        return components

    def explains(self, mz_value, ppm_tolerance):
        """
        Determine if the mz_value can be explained by a linear combination of components.

        Args:
            mz_value (float): Observed m/z value.
            ppm_tolerance (float): Error tolerance in parts-per-million.

        Returns:
            list: Solutions that explain the mz_value within the given tolerance.
        """
        # Initialize lazy DP
        
        print("Len Cache: ", len(self.lazy_dp.dp_cache))
        return self.lazy_dp.query(mz_value, ppm_tolerance)

# Example Usage
if __name__ == "__main__":
    X, Y = [], []
    # Initialize the MassExplainer with the provided CSV file
    explainer = MassExplainer("./isotopes.csv", max_solutions=5)
    import time
    i = 0
    while True:
        print(i)
        i += 1
        t1 = time.time()
        c = explainer.explains(i*12, 5)
        for z in c:
            print("\t", i*12, z)
        t2 = time.time()
        print("Time Elapsed: ", str(t2-t1), " mz: ", str(i * 12)) 
        #print("Len Cache: ", str(len(explainer.dp_cache)))
        X.append(i)
        Y.append(t2-t1)
        if (t2-t1) > 10.0:
            break

    import matplotlib.pyplot as plt
    plt.scatter(X, Y)
    plt.show()
