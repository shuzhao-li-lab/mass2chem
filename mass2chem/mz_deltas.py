'''
This module provides functions to retrieve the data from source data (immutable) and lib (ongoing lists).
1. Source Data:
  1. Top N mz_deltas for a given instrument and ionization mode.mz deltas can be fragments, isotopes, or modifications. Various filters are given to retrieve the desired mz_deltas. The functions are:
      1. top_N_modication_mz_deltas_for_instrument_and_mode()
      2. top_N_isotope_mz_deltas_for_instrument_and_mode()
  2. known xenobiotic modifications as described in Zhao, Haoqi Nina, et al with function: zhao2024_drug_exposure()
  3. known biotransformations as described in Xing, Shipei, et al with function: xing2020_hypothetical_neutral_losses()
2. lib:
  1. lib_mzdiff_bioreaction: List of common mass differences in bioreactions from difference sources. Function: lib_mzdiff_bioreaction()
  2. lib_mzdiff_in_source: List of common mass differences in in-source fragmentation if different instruments/ion modes. Function: lib_mzdiff_in_source()
'''

import os
import json
import pkg_resources
import pandas as pd

def harmonize_instrument_mode(instrument, mode):
    # todo - this should be moved to a common asari utils repository
    """
    Users may not call the instruments the exact string we do. 
    This function helps to harmonize instrument type and ionization 
    modes to our internal canonical names. Certainly, this function 
    needs to be somewhere else. 

    Canonical instrument types:
        "orbi", "tof"

    Canonical ionization modes:
        "pos", "neg"
    """

    instrument, mode = instrument.lower(), mode.lower()
    harmonize_instrument = {
        "orbitrap": "orbi",
        "orbital": "orbi",
        "tof": "tof"
    }
    harmonize_mode = {
        "positive": "pos",
        "negative": "neg",
        "+": "pos",
        "-": "neg"
    }
    if instrument not in harmonize_instrument.values():
        if instrument in harmonize_instrument:
            instrument = harmonize_instrument[instrument]
        else:
            raise ValueError(f"The instrument {instrument} is not supported and could not be mapped to a supported instrument type!\n see: {__doc__}")
    if mode not in harmonize_mode.values():
        if mode in harmonize_mode:
            mode = harmonize_mode[mode]
        else:
            raise ValueError(f"The ionization mode of {mode} is unknown, please provide 'pos' or 'neg'\n see: {__doc__}")
    return instrument, mode

def retrieve_frequent_deltas(instrument, mode):
    """
    Given an instrument type and mode, both strings (see __harmonize_instrument_mode), return the list of frequent mz deltas
    as described in ___ref_to_preprint___ and saved in the source_data folder.
    """
    instrument, mode = harmonize_instrument_mode(instrument, mode)
    frag_lists_by_instrument_by_mode = {
        "orbi": {
            "pos": pkg_resources.resource_filename('mass2chem', 'source_data/chi2025_isf_orbi_neg.tsv'),
            "neg": pkg_resources.resource_filename('mass2chem', 'source_data/chi2025_isf_orbi_pos.tsv')
        },
        "tof": {
            "pos": pkg_resources.resource_filename('mass2chem', 'source_data/chi2025_isf_tof_pos.tsv'),
            "neg": pkg_resources.resource_filename('mass2chem', 'source_data/chi2025_isf_tof_neg.tsv')
        }
    }
    if frag_lists_by_instrument_by_mode.get(instrument, {}).get(mode, None) is None:
        raise ValueError(f"The combination of {instrument} and {mode} is not supported")
    if not os.path.exists(frag_lists_by_instrument_by_mode[instrument][mode]):
        raise FileNotFoundError(f"The file {frag_lists_by_instrument_by_mode[instrument][mode]} does not exist")
    deltas_for_intrument_mode = pd.read_csv(frag_lists_by_instrument_by_mode[instrument][mode], sep='\t', encoding='unicode_escape')
    # these may be removed later if we change format, the above exceptions will not change.
    assert "delta_mz" in deltas_for_intrument_mode.columns, "The column name 'delta_mz' is not found in the source data"
    assert "type" in deltas_for_intrument_mode.columns, "The column name 'type' is not found in the source data"
    return deltas_for_intrument_mode

def top_N_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N=20, filter_type=None, sort_by='count_estimate'):
    """
    For a given instrument type and ionization mode, this function retrieves the top N mz_deltas from the source data
    Optionally, you can filter the results by the type of mz_deltas (e.g., 'modification', 'isotope')
    Instrument should be either 'orbitrap' or 'tof'
    Ionization mode should be either 'pos' or 'neg'
    """
    mz_deltas = retrieve_frequent_deltas(instrument, ionization_mode)
    if sort_by:
        if sort_by not in mz_deltas.columns:
            raise ValueError(f"The column name {sort_by} is not found in the source data")
        mz_deltas = mz_deltas.sort_values(by=sort_by, ascending=False)
    if filter_type:
        mz_deltas = mz_deltas[mz_deltas['type'] == filter_type]
    return mz_deltas.head(N)

def top_N_modication_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N=20):
    """
    This is a utility function for retrieving the top N modification mz_deltas from the source data
    """
    return top_N_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N=N, filter_type="modification")

def top_N_isotope_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N=20):
    """
    This is a utility function for retrieving the top N isotope mz_deltas from the source data
    """
    return top_N_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N=N, filter_type="isotope")

def xing2020_hypothetical_neutral_losses():
    """
    This function returns a dataframe of known biotransformations as described in Xing, Shipei, et al. "Retrieving and utilizing hypothetical neutral losses from tandem mass 
    spectra for spectral similarity analysis and unknown metabolite annotation." Analytical Chemistry 92.21 (2020): 14476-14483. The data is saved in the source_data folder.
    """
    assert os.path.exists(pkg_resources.resource_filename('mass2chem', 'source_data/xing2020_hypothetical_neutral_losses.tsv')), "The file xing2020_hypothetical_neutral_losses.tsv does not exist"
    return pd.read_csv(pkg_resources.resource_filename('mass2chem', 'source_data/xing2020_hypothetical_neutral_losses.tsv'), sep = "\t")

def zhao2024_drug_exposure(skip_isotopologues=True):
    """
    This function returns a dataframe of known xenobiotic modifications as described in Zhao, Haoqi Nina, et al. 
    "Empirically establishing drug exposure records directly from untargeted metabolomics data." bioRxiv (2024).

    The input list contains isotopologue differences, which are not modifications. 
    If you want to include them, set skip_isotopologues=False. 
    """
    assert os.path.exists(pkg_resources.resource_filename('mass2chem', 'source_data/zhao2024_drug_exposure.tsv')), "The file zhao2024_drug_exposure.tsv does not exist"
    xeno_mods = pd.read_csv(pkg_resources.resource_filename('mass2chem', 'source_data/zhao2024_drug_exposure.tsv'), sep = "\t")
    
    if skip_isotopologues:
        # fix this more elegantly
        filtered = []
        for xm in xeno_mods.to_dict(orient='records'):
            if isinstance(xm['Rationale'], str):
                if 'isotope' in xm['Rationale'] or '15N(1)' in xm['Rationale']:
                    pass
                else:
                    filtered.append(xm)
            else:
                filtered.append(xm)
        xeno_mods = pd.DataFrame(filtered)
    return xeno_mods

def lib_mzdiff_bioreaction():
    """Function used for retrieving mzdiff bioreactions under lib folder. 

    Returns
    -------
    dict
        A python dict containing mzdiff bioreactions deltas from different sources. 
        E.g. {'zhao2024_drug_exposure': [[-212.0086, 'denucleotide (monophosphate)', "{'C': -5, 'H': -9, 'O': -7, 'P': -1}"], 
        [-207.068414, 'de-3-phenoxyphenylacetonitrile', "{'C': -14, 'H': -9, 'N': -1, 'O': -1}"]]}
        Three elements in the list are mass delta signatures, description and formula dict respectively. 

    Examples
    --------
    >>> import mass2chem.mz_deltas
    >>> mass2chem.mz_deltas.lib_mzdiff_bioreaction()
    """
    try:
        with open(pkg_resources.resource_filename('mass2chem', 'lib/mzdiff_bioreaction.json'), 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print("Error: mzdiff_bioreaction.json file not found.")
        return {}  
    except json.JSONDecodeError:
        print("Error: Failed to parse mzdiff_bioreaction.json (Invalid JSON format).")
        return {}  
    except Exception as e:
        print(f"Unexpected error: {e}")
        return {}
    
def lib_mzdiff_in_source():
    """Function used for retrieving mzdiff of in-source fragments under lib folder. 

    Returns
    -------
    dict
        A python dict containing mzdiff of in-source fragments of different instruments and ion modes. 
        E.g. {'orbi_pos': [[1.0033, 4155, 1.00335, 'C isotope', "{'(C12)': -1, '(C13)': 1}", 'isotope'], 
        [14.0155, 900, 14.015649, 'addition of acetic acid and loss of CO2. Reaction: (+C2H2O2) and (-CO2)', "{'C': 1, 'H': 2}", 'modification']]}
        Five elements in the list are experimental deltas, count estimate, mass delta signatures in database, description and formula dict respectively. 

    Examples
    --------
    >>> import mass2chem.mz_deltas
    >>> mass2chem.mz_deltas.lib_mzdiff_in_source()
    """
    try:
        with open(pkg_resources.resource_filename('mass2chem', 'lib/mzdiff_in_source.json'), 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print("Error: mzdiff_in_source.json file not found.")
        return {}  
    except json.JSONDecodeError:
        print("Error: Failed to parse mzdiff_in_source.json (Invalid JSON format).")
        return {}  
    except Exception as e:
        print(f"Unexpected error: {e}")
        return {}