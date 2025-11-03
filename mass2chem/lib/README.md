# lib
> This folder contains curated data from different sources. The files will keep changing as time goes by.

The mzdiff_bioreaction and mzdiff_in_source JSON files are not cleanly separated. 

The isotope, adduct and ISF lists used in Chi et al are in `isf_default_mzdiff.py`.


| File    | Note | Usage |
| :------ | :----: | :--- |
| common_mass.py |   List of most common mass differences (mostly but not limited to adducts) | All common scenarios needs mass difference comparison |
| formula_coordinate.py |  Common formulas and their masses | Mass calibration |
| LCMS_contaminants.py  |  Common contaminants in LCMS.   | Contaminants detection |
| mzdiff_bioreaction.json |  List of common mass differences in bioreactions from difference sources  | Intended for Bioreactions exploration. Can be retrieved by mz_deltas.py |
| mzdiff_in_source.json |  List of common mass differences by instruments/ion modes  | From In-source fragmentation study, including isotopes, adducts and ISFs. Can be retrieved by mz_deltas.py |
| reaction_rules.py |  Rules of chemical reactions based on mass differences  |  |