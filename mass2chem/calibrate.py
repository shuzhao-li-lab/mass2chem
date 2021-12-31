'''
Computational calibration of meassured masses,
assuming that a common set of metabolites exist in the samples.
The chance of matching the accurate mass of large number of ions is unambiguous.

Multiple methods to calibrate:

    1. A regression model is built for (theoretical m/z) = a * (measured m/z) + b + e,
    then the model is applied to all features to calculate a calibrated mz.
    Because the model is too sensitive to `a` and mass shift in this slope is often negligible,
    method 2 is preferred.

    2. (measured m/z)/(theoretical m/z) - 1 =  b + e. Systematic shift can be reflected in `b`.
    b*10e6 is conventional ppm.
    Assuming a normal distribution of mass errors, b and std are fitted.

    3. (theoretical m/z) = b + e*f(measured m/z). The accuracy is a function of m/z, as 
    larger m/z has larger errors. The convention is to use a dalton window (0 order) or ppm (1st order),
    but Orbitrap data suggest that we need a 2nd order model. 

We implement method 2 for now. Method 3 will be on to-do list.

Suggested use:
from .calibrate import normal_error_distribution as mass_calibrate


These have to be based on high-selectivity ions, otherwise ground truth is hard to establish,
and the statistics is not valid.

This should be done after correspondence of high-selectivity features 
(after peak detection and initial annotation).
If the mass calibration is off, correct shift (mu) and use empirical stdev to recalculate selectivity.



'''

from scipy.stats import norm

from .io import index_features

def mass_calibrate(mz_List_1, mz_List_2):
    '''
    mz_List_1, mz_List_2 as numpy arrays.
    return
    ------
    The mean (List2 - List1) and standard deviation of mass differences. 
    The latter can be used as a guide of precision in an experiment.

    Example
    -------
    p = norm.pdf(xx, mu, std)
    plt.plot(xx, p, 'kx', linewidth=2)

    '''
    return norm.fit(mz_List_2 - mz_List_1)



# -----------------------------------------------------------------------------
#
# Prototype code with short-handed data structure
#
# -----------------------------------------------------------------------------

def old_mass_calibrate(query_features, anchor_features, limit_ppm=10):
    '''
    Use a list of common metabolites to calibrate m/z values in query_features 
    (e.g. all m/z values in a smaple).
    Assuming a large number of compounds in ref_mass_dict exist in the list_of_features,
    which should be the case in common biological matrices.
    This function does not require complete use of anchor_features.
    The mass errors are assumed to fit a normal distribution, 
    and it's parameters (e.g. 2*stdev) are used for subsequent confidence in matching.
    
    Input
    -----
    query_features: List of m/z values from the same Experiment.
                
    anchor_features: list of m/z values for metabolites known to be in the sample.
                    They need not to be 100% present, but in a sufficient number 
                    to estimate error distribution.

    limit_ppm: default 10. Expecting a high-res instrument performs within 25 ppm, 

    Return
    ------
    coefficients: (mean, std) in the fitted model of normal distribution of mass errors.
    updated_query_features: List of updated m/z values for the query features in same order.
    '''
    matched = []
    for m in calibration_mass_dict.values():
        best_pair = _find_closest(m, indexed_features, limit_ppm)
        if best_pair:
            matched.append(best_pair)

    def _find_closest(query_mz, indexed_features, limit_ppm):
        # to find closest match of a theoretical m/z in indexed_features, under limit_ppm
        # return (delta, feature_mz, query_mz) or None
        _delta = query_mz * limit_ppm * 0.000001
        _low_lim, _high_lim = query_mz - _delta, query_mz + _delta
        result = []
        for ii in range(int(_low_lim), int(_high_lim)+1):
            for F in indexed_features[ii]:
                result.append( (abs(query_mz-F['mz']), F['mz'], query_mz) )

        result.sort()
        if result and result[0][0] < _delta:
            return result[0]
        else:
            return None

    def _objective(x, a, b): return a * x + b

    indexed_features = index_features(list_features)
    matched = []
    for m in calibration_mass_dict.values():
        best_pair = _find_closest(m, indexed_features, limit_ppm)
        if best_pair:
            matched.append(best_pair)

    num_matched = len(matched)
    if num_matched < 3:
        print("Aborted. Found %d matched mass features." %num_matched)

    else:
        if num_matched < 10:
            print("Warning. Using %d matched mass features, too few for regression." %num_matched)
        X, Y = [x[1] for x in matched], [x[2] for x in matched]
        popt, _ = curve_fit(_objective, X, Y)
        a, b = popt
        if abs(b) > 0.01 or abs(a-1) > 0.000025:
            print(" -- WARNING -- extreme parameters are usually bad, ", a, b)

        for F in list_features:
            F['calibrated_mz'] = a * F['mz'] + b

        return (a,b), list_features


def fake_mass_calibrate(list_features):
    for F in list_features:
        F['calibrated_mz'] = F['mz']
    return list_features
