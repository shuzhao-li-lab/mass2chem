import unittest
import pandas as pd

from mass2chem.mz_deltas import (
    harmonize_instrument_mode,
    retrieve_frequent_deltas,
    top_N_mz_deltas_for_instrument_and_mode,
    top_N_modication_mz_deltas_for_instrument_and_mode,
    top_N_isotope_mz_deltas_for_instrument_and_mode,
    known_biological_modifications,
    known_xenobiotic_modifications
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

    def test_known_biological_modifications(self):
        df = known_biological_modifications()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("Mass difference", df.columns)
        self.assertIn("Description", df.columns)
        self.assertIn("Formula Change", df.columns)

    def test_known_xenobiotic_modifications(self):
        df = known_xenobiotic_modifications()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("Mass difference", df.columns)
        self.assertIn("Rationale", df.columns)
        self.assertIn("Formula Change", df.columns)
        df_with_isotopologues = known_xenobiotic_modifications(skip_isotopologues=False)
        self.assertGreater(len(df_with_isotopologues), len(df))

if __name__ == '__main__':
    unittest.main()