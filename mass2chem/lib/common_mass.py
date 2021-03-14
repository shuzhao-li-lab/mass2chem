# Common contaminants, interfering ions, adducts etc. will need to be sync to Azimuth.
#

PROTON = 1.00727646677

# atomic masses adopted from the National Institute of Standarts and Technology (NIST)
#   (http://physics.nist.gov/PhysRefData/Compositions/index.html)
atoms = {
    'H': 
}


Carbon-12	C-12	12.000000
Carbon-13	C-13	13.003355
Hydrogen	H	1.007825
Deuterium	D	2.014101
Oxygen	O	15.994915
Oxygen-18	O-18	17.999160
Nitrogen-14	N-14	14.003074
Nitrogen-15	N-15	15.000109
Sulfur-32	S-32	31.972071
Sulfur-34	S-34	33.967867
Phosphorus	P	30.973762
Silicium	Si	27.976927
Sodium	Na	22.989769
Potassium	K	38.963706
Chlorine-35	Cl-35	34.968853
Chlorine-37	Cl-37	36.965903
Bromine-79	Br-79	78.918337
Bromine 81	Br-81	80.916291
Iodine	Iodine	126.904473
Fluorine	F	18.998403
Iron-54	Fe-54	53.939611
Iron-56	Fe-56	55.934938
Lithium-6	Li-6	6.015123
Lithium-7	Li-7	7.016005
Boron-10	B-10	10.012937
Boron-11	B-11	11.009305
Copper-63	Cu-63	62.929601
Copper-65	Cu-65	64.927794
Silver-107	Ag-107	106.905094
Silver-109	Ag-109	108.904756
Tin-120	Sn-120	119.902199
Cesium	Cs-133	132.905450
		
Proton	H+	1.007276
Electron	e-	0.000549





isotope_shifts = {
    'C13': 1.0034,
    '2C13': 2.0068,
    'S34': 1.9958,
    'Cl37': 1.9972,

}

adduct_shifts_pos = {
    'M+Na[1+]': 22.9893,

}

adduct_shifts_neg = {
    'M+Cl[-]': 34.9689,

}






'''
    addList = [(mw, 'M[1+]', ''), 
         (mw + PROTON, 'M+H[1+]', ''),
         (mw/2 + PROTON, 'M+2H[2+]', ''),
         (mw/3 + PROTON, 'M+3H[3+]', ''),
         (mw +1.0034 + PROTON, 'M(C13)+H[1+]', 'C'),
         (mw/2 + 0.5017 + PROTON, 'M(C13)+2H[2+]', 'C'),
         (mw/3 + 0.3344 + PROTON, 'M(C13)+3H[3+]', 'C'),
         (mw +1.9958 + PROTON, 'M(S34)+H[1+]', 'S'),
         (mw +1.9972 + PROTON, 'M(Cl37)+H[1+]', 'Cl'),
         (mw + 21.9820 + PROTON, 'M+Na[1+]', ''),        # Na = 21.9820 + PROTON = 22.9893
         (mw/2 + 10.991 + PROTON, 'M+H+Na[2+]', ''),
         (mw + 37.9555 + PROTON, 'M+K[1+]', ''),         # K = 37.9555 + PROTON = 38.9628
         (mw + 18.0106 + PROTON, 'M+H2O+H[1+]', ''), 
         (mw - 18.0106 + PROTON, 'M-H2O+H[1+]', 'H2O'), 
         (mw - 36.0212 + PROTON, 'M-H4O2+H[1+]', 'H4O2'),
         (mw - 17.0265 + PROTON, 'M-NH3+H[1+]', 'NH3'),
         (mw - 27.9950 + PROTON, 'M-CO+H[1+]', 'CO'),
         (mw - 43.9898 + PROTON, 'M-CO2+H[1+]', 'CO2'),
         (mw - 46.0054 + PROTON, 'M-HCOOH+H[1+]', 'H2CO2'),
         (mw + 67.9874 + PROTON, 'M+HCOONa[1+]', ''),
         (mw - 67.9874 + PROTON, 'M-HCOONa+H[1+]', 'HCO2Na'),
         (mw + 57.9586 + PROTON, 'M+NaCl[1+]', ''), 
         (mw - 72.0211 + PROTON, 'M-C3H4O2+H[1+]', 'C3H4O2'),
         (mw + 83.9613 + PROTON, 'M+HCOOK[1+]', ''),
         (mw - 83.9613 + PROTON, 'M-HCOOK+H[1+]', 'HCO2K'),
         ] + [
             (mw - PROTON, 'M-H[-]', ''),
        (mw/2 - PROTON, 'M-2H[2-]', ''),
        (mw + 1.0034 - PROTON, 'M(C13)-H[-]', 'C'),
        (mw + 1.9958 - PROTON, 'M(S34)-H[-]', 'S'),
        (mw + 1.9972 - PROTON, 'M(Cl37)-H[-]', 'Cl'),
        (mw + 22.9893 - 2*PROTON, 'M+Na-2H[-]', ''),
        (mw + 38.9628 - 2*PROTON, 'M+K-2H[-]', ''),
        (mw - 18.0106 - PROTON, 'M-H2O-H[-]', 'H2O'),
        (mw + 34.9689, 'M+Cl[-]', ''),
        (mw + 36.9659, 'M+Cl37[-]', ''),
        (mw + 78.9183, 'M+Br[-]', ''),
        (mw + 80.9163, 'M+Br81[-]', ''),
        (mw + 2*12 + 3*1.007825 + 14.00307 - PROTON, 'M+ACN-H[-]', ''),
        (mw + 1.007825 + 12 + 2*15.99491, 'M+HCOO[-]', ''),
        (mw + 3*1.007825 + 2*12 + 2*15.99491, 'M+CH3COO[-]', ''),
        (mw - PROTON + 15.99491, 'M-H+O[-]', ''),
        ]

'''
