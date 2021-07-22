# mass2chem - common utilities in interpreting mass spectrometry data

Annotation and Inferrence

## to include

* Handling chemical formula via pychemy (which needs update to high res data)

* A list of common mass values, including contaminants

* A list of common adducts, rules, while they are more directly ready in future Azimuth

* Chemical similary computing

* Reaction inference, including mass diff corresponding to common reactions

* Annotation via in-house libraries 

* hook/adaptor to other tools


## Added basic formula based calculations

Note: RE based formula parsing is still limited.

Pychemy isn't good fit, as 
1) high-resolution calculation needs update
2) Open babel binding is not worthy the trouble

E.g. in pychemy.adducts, it's wrong to use ('M+3H', 0.33,  1.0073),
because the computing error in 0.33 (correct is 1/3) is far too large for mass precision.

Pychemy is included as stripped version in "mass2chem.chem" for now, but only used for formula handling.

For high-resolution measurements, electrons should be considered too.


## Related

https://github.com/shuzhao-li/pychemy

https://github.com/opencobra/cobrapy/blob/devel/cobra/core/formula.py (they are using average molecular weight, not mass spec oriented)

https://github.com/shuzhao-li/Azimuth 

