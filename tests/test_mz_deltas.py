import unittest
import pandas as pd

from mass2chem.mz_deltas import (
    harmonize_instrument_mode,
    retrieve_frequent_deltas,
    top_N_mz_deltas_for_instrument_and_mode,
    top_N_modication_mz_deltas_for_instrument_and_mode,
    top_N_isotope_mz_deltas_for_instrument_and_mode,
    xing2020_hypothetical_neutral_losses,
    zhao2024_drug_exposure,
    lib_mzdiff_bioreaction,
    lib_mzdiff_in_source
)

class TestMzDeltas(unittest.TestCase):

    def test_harmonize_instrument_mode(self):
        self.assertEqual(harmonize_instrument_mode("orbitrap", "positive"), ("orbi", "pos"))
        self.assertEqual(harmonize_instrument_mode("tof", "negative"), ("tof", "neg"))
        self.assertEqual(harmonize_instrument_mode("orbital", "+"), ("orbi", "pos"))
        with self.assertRaises(ValueError):
            harmonize_instrument_mode("unknown", "positive")
        with self.assertRaises(ValueError):
            harmonize_instrument_mode("orbitrap", "unknown")

    def test_retrieve_frequent_deltas(self):
        df = retrieve_frequent_deltas("orbitrap", "positive")
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("delta_mz", df.columns)
        self.assertIn("type", df.columns)

    def test_top_N_mz_deltas_for_instrument_and_mode(self):
        df = top_N_mz_deltas_for_instrument_and_mode("orbitrap", "positive", N=10)
        self.assertEqual(len(df), 10)
        df = top_N_mz_deltas_for_instrument_and_mode("orbitrap", "positive", N=5, filter_type="modification")
        self.assertEqual(len(df), 5)
        self.assertTrue((df['type'] == "modification").all())

    def test_top_N_modication_mz_deltas_for_instrument_and_mode(self):
        df = top_N_modication_mz_deltas_for_instrument_and_mode("orbitrap", "positive", N=10)
        self.assertEqual(len(df), 10)
        self.assertTrue((df['type'] == "modification").all())

    def test_top_N_isotope_mz_deltas_for_instrument_and_mode(self):
        df = top_N_isotope_mz_deltas_for_instrument_and_mode("orbitrap", "positive", N=10)
        self.assertEqual(len(df), 7)
        self.assertTrue((df['type'] == "isotope").all())

    def test_xing2020_hypothetical_neutral_losses(self):
        df = xing2020_hypothetical_neutral_losses()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("Mass difference", df.columns)
        self.assertIn("Description", df.columns)
        self.assertIn("Formula Change", df.columns)

    def test_zhao2024_drug_exposure(self):
        df = zhao2024_drug_exposure()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("Mass difference", df.columns)
        self.assertIn("Rationale", df.columns)
        self.assertIn("Formula Change", df.columns)
        df_with_isotopologues = zhao2024_drug_exposure(skip_isotopologues=False)
        self.assertGreater(len(df_with_isotopologues), len(df))
    
    def test_lib_mzdiff_bioreaction(self):
        d = lib_mzdiff_bioreaction()
        sources = d.keys()
        self.assertIsInstance(d, dict)
        self.assertIn('zhao2024_drug_exposure', sources)
        self.assertIn('xing2020_hypothetical_neutral_losses', sources)
        
    def test_lib_mzdiff_in_source(self):
        d = lib_mzdiff_in_source()
        methods = d.keys()
        self.assertIsInstance(d, dict)
        self.assertIn('orbi_pos', methods)
        self.assertIn('orbi_neg', methods)
        self.assertIn('tof_pos', methods)
        self.assertIn('tof_neg', methods)

if __name__ == '__main__':
    unittest.main()