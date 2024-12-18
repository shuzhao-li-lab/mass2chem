import numpy as np
from mass2chem.formula import dict_to_hill_formula
import pandas as pd

class MassExplainer:
    """
    A class to explain a given target mass by identifying possible combinations of isotopic elements.

    Attributes:
        dp (dict): Memoized solutions for dynamic programming.
        bad (set): Known invalid states to avoid redundant calculations.
        masses (np.array): Vector of isotopic masses (in milliDaltons).
        naps (np.array): Vector of isotopic abundances (NAP values).
        valence (np.array): Vector of valence electron counts for each isotope.
        max_counts (np.array): Maximum allowed counts for each isotope.
        elements (list): List of element symbols (e.g., ["C", "H", "O"]).
        tolerance (float): Mass tolerance for solutions.
        allow_negative (bool): Allow negative counts for isotopes (neutral loss explanation).
        skip_valence (bool): Disable valence check on combinations.
        __configured (bool): Indicates if the class is configured and ready to use.
    """
    def __init__(self) -> None:
        """Initializes the MassExplainer with default values."""
        self.dp = {}  # Memoization dictionary for solutions.
        self.bad = set()  # Set to track invalid states.

        # Isotopic and element-specific properties
        self.masses = None
        self.naps = None
        self.valence = None
        self.max_counts = None
        self.elements = None

        # Configurable parameters
        self.tolerance = 0.01  # Mass tolerance for solutions.
        self.allow_negative = False  # Allow negative counts for isotopes.
        self.skip_valence = False  # Skip valence check.

        self.__configured = False  # Configuration status.

    def explain(self, target):
        """
        Finds combinations of isotopic elements that explain the target mass.

        Parameters:
            target (float): Target mass to explain.

        Returns:
            list: List of valid combinations or an empty list if none exist.
        """
        if not self.__configured:
            raise Exception("Please call configure first")
        key = (0, self.tolerance, int(target * 1000))

        def backprop(current_sum, consider=0):
            """
            Recursive helper function to find valid combinations of isotopes.

            Parameters:
                current_sum (int): Current mass to explain (scaled by 1000).
                consider (int): Index of the current isotope under consideration.

            Returns:
                np.array: Array of valid combinations for the given state.
            """
            key = (consider, self.tolerance, current_sum)
            if key in self.dp:
                return self.dp[key]
            elif abs(current_sum) < 1000 * self.tolerance:
                return np.array([[0 for _ in self.masses]], dtype=np.int8)
            elif key in self.bad or consider >= len(self.masses) or (not self.allow_negative and current_sum < 1000 * self.tolerance):
                return []
            else:
                mass = self.masses[consider]
                max_coin = int(min(np.ceil((current_sum + 1)/mass), self.max_counts[consider], 127))
                for i in range(-max_coin if self.allow_negative else 0, max_coin):
                    for comb in backprop(current_sum - mass * i if i != 0 else current_sum, consider + 1):
                        new_comb = np.array(comb, copy=True, dtype=np.int8)
                        new_comb[consider] = i
                        if consider == 0 and np.sum(np.abs(new_comb[1:])) == 0 and i == 2:
                            if key not in self.dp:
                                self.dp[key] = []
                            self.dp[key].append(new_comb)
                        elif (np.prod(np.power(self.naps, np.abs(new_comb))) > 0.4 and consider != 0) or (self.skip_valence or (np.dot(self.valence, np.abs(new_comb)) - (np.sum(new_comb[1:]) * np.sum(new_comb[1:]) > 1) >= 0)):
                            if key not in self.dp:
                                self.dp[key] = []
                            self.dp[key].append(new_comb)
                if key in self.dp:
                    self.dp[key] = np.array(self.dp[key], dtype=np.int8)
                    return self.dp[key]
                else:
                    self.bad.add(key)
                    return []

        backprop(int(target * 1000))
        return self.dp.get(key, [])
    
    def load_isotopes(self, csv_path="/Users/mitchjo/mass2chem/mass2chem/utils/isotopes.csv"):
        """
        Loads isotopic data from a CSV file.

        Parameters:
            csv_path (str): Path to the isotopes CSV file.

        Returns:
            tuple: Tuple containing masses, valence, naps, max_counts, and elements arrays.
        """
        iso_data = pd.read_csv(csv_path)
        masses = np.array(iso_data['mz_delta'] * 1000, dtype=np.int32)
        valence = np.array(iso_data['valence'], dtype=np.int8)
        naps = np.array(iso_data['prob'], dtype=np.float32)
        maxes = np.array(iso_data['max_val'], dtype=np.int8)
        elements = np.array(iso_data['formula'])

        W = [z for z in maxes]
        for i in range(len(masses)):
            for j in range(maxes[i]):
                if naps[i] ** j < 0.4:
                    break
            W[i] = j
        maxes = W

        return masses, valence, naps, maxes, elements
    
    def configure(self, tolerance=0.01, allow_negative=False):
        """
        Configures the MassExplainer for use.

        Parameters:
            tolerance (float): Mass tolerance for solutions.
            allow_negative (bool): Allow negative counts for isotopes.
        """
        self.dp = {} 
        self.bad = set()
        masses, valence, naps, maxes, elements = self.load_isotopes()
        self.masses = masses
        self.valence = valence
        self.naps = naps
        self.max_counts = maxes
        self.elements = elements
        self.tolerance = tolerance
        self.allow_negative = allow_negative
        self.__configured = True

    def explain_neutral_mass(self, qmass):
        """
        Once configured, explain a neutral mass using the provided components

        Complexity is the sum of the added and removed elements divided by the gcd of the non-zero elements. 

        Parameters:
            qmass (float): Query neutral mass

        Returns:
            to_return (list): list of solutions in form [qmass, component_list, solution_mass, solution_mass_error, complexity, nap]
        """
        to_return = []
        seen = set()

        r1 = self.explain(qmass) # search the qmass
        r2 = self.explain(qmass - self.tolerance) # add the tolerance and search to account for rounding
        r3 = self.explain(qmass + self.tolerance) # sub the tolerance and serach to account for rounding

        if len(r1) or len(r2) or len(r3):
            result = np.unique(np.vstack([x for x in [r1, r2, r3] if len(x)]), axis=0)
            if result.shape[1]:
                for c in sorted(result, key=lambda x: (abs((np.dot(self.masses, x)/1000)-qmass), np.sum(abs(x)) / np.gcd.reduce(np.abs(x)))):
                    c_formula = dict_to_hill_formula({k: v for k, v in zip(self.elements, c) if v != 0})
                    c_mass = np.dot(self.masses, c) / 1000
                    if c_formula not in seen and abs(c_mass - qmass) < self.tolerance:
                        to_return.append([qmass, c.tolist(), c_formula, float(c_mass), float(round(c_mass - qmass, 6)), np.sum(abs(c)) / np.gcd.reduce(np.abs(c)), np.prod(np.power(self.naps, np.abs(c)))])
                        seen.add(c_formula)
        return to_return
    
    @staticmethod
    def _explain_figure1C(path):
        """
        Given a csv with neutral loss patterns (such as figure 1C from 12/18/2024), generate explanations

        Args:
            path (str): path to csv
        """
        rdf = []
        explainer = MassExplainer()
        explainer.configure(tolerance=0.01, allow_negative=True)
        for qm in sorted(set(pd.read_csv(path)["query_mass"])):
            results = explainer.explain_neutral_mass(qm)
            if results:
                for r in results:
                    qmass, clist, formula, mass, error, complexity, nap = r
                    rdf.append({
                        "query_mass": qmass,
                        "element_list": clist,
                        "formula": formula,
                        "solution_mass": mass,
                        "solution_error": error,
                        "complexity": complexity,
                        "nap": nap
                    })
            else:
                rdf.append({
                    "query_mass": qmass
                })
        return pd.DataFrame(rdf).to_csv(path.replace(".csv", "_v3.csv"))

if __name__ == '__main__':
    import sys
    MassExplainer._explain_figure1C(sys.argv[1])