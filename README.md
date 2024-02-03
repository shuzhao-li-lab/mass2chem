# mass2chem - low level utilities in interpreting mass spectrometry data

This package provides 
- functions on handling chemical formulas
- formula based adduct calculation 
- indexing and search functions on mass spec data
- libraries of common metabolites, contaminants, mass differences
- [to-do] functions of chemical similary, dataset similarity

## Related tools
- Generalized computing of isotopes and adducts: khipu (https://github.com/shuzhao-li-lab/khipu, https://pubs.acs.org/doi/10.1021/acs.analchem.2c05810)

- High-level metabolite functions and metabolic models: Json's Metabolite Services (JMS, https://github.com/shuzhao-li-lab/JMS)

- Metabolomics data processing: asari (https://github.com/shuzhao-li-lab/asari, https://www.nature.com/articles/s41467-023-39889-1)

- Python-Centric Pipeline for Metabolomics (https://github.com/shuzhao-li-lab/PythonCentricPipelineForMetabolomics)

- Common data models for metabolomics: metDataModel (https://github.com/shuzhao-li/metDataModel)

## Third party references:

https://github.com/opencobra/cobrapy/blob/devel/cobra/core/formula.py (using average molecular weight at the time of retrieval, not mass spec oriented)

https://github.com/domdfcoding/chemistry_tools

Pychemy (https://github.com/ginkgobioworks/pychemy). 
Pychemy at this time isn't good fit for high-resolution metabolomics because its mass calculation is not of enough precision. E.g. in pychemy.adducts, it's wrong to use ('M+3H', 0.33,  1.0073),
because the computing/rounding error in 0.33 (correct is 1/3) is far too large for mass precision.
For high-resolution measurements, electrons should be considered too.

------------------------
Please do not hesitate to contact us via the GitHub issues.

Citation to come.
