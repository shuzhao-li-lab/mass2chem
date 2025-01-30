#
# This module provides functions to retrieve the top N mz_deltas for a given instrument and ionization mode
# mz deltas can be fragments, isotopes, or modifications. Various filters are given to retrieve the desired mz_deltas
#

import pkg_resources
import os
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
            "pos": pkg_resources.resource_filename('mass2chem', 'source_data/S2_top_frequent_delta_mz_orbi_pos.tsv'),
            "neg": pkg_resources.resource_filename('mass2chem', 'source_data/S3_top_frequent_delta_mz_orbi_neg.tsv')
        },
        "tof": {
            "pos": pkg_resources.resource_filename('mass2chem', 'source_data/S4_top_frequent_delta_mz_tof_pos.tsv'),
            "neg": pkg_resources.resource_filename('mass2chem', 'source_data/S5_top_frequent_delta_mz_tof_neg.tsv')
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
    return top_N_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N, "modification")

def top_N_isotope_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N=20):
    """
    This is a utility function for retrieving the top N isotope mz_deltas from the source data
    """
    return top_N_mz_deltas_for_instrument_and_mode(instrument, ionization_mode, N, "isotope")

def known_biological_modifications():
    """
    This function returns a dataframe of known biotransformations as described in ___ref_to_preprint___, sourced
    from Tau Huan's paper <reference to paper>. The data is saved in the source_data folder.
    """
    assert os.path.exists(pkg_resources.resource_filename('mass2chem', 'source_data/known_biological_modifications.tsv')), "The file S1_known_biotransformations.tsv does not exist"
    return pd.read_csv(pkg_resources.resource_filename('mass2chem', 'source_data/known_biological_modifications.tsv'), sep = "\t")

def known_xenobiotic_modifications(skip_isotopologues=True):
    """
    This function returns a dataframe of known xenobiotic modifications as described in ___ref_to_preprint___, sourced
    from _________.

    The input list contains isotopologue differences, which are not modifications. 
    If you want to include them, set skip_isotopologues=False. 
    """
    assert os.path.exists(pkg_resources.resource_filename('mass2chem', 'source_data/known_xenobiotic_modifications.tsv')), "The file S1_known_biotransformations.tsv does not exist"
    xeno_mods = pd.read_csv(pkg_resources.resource_filename('mass2chem', 'source_data/known_xenobiotic_modifications.tsv'), sep = "\t")
    
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