'''
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


'''
# Common contaminants, interfering ions, adducts etc. will need to be sync to Azimuth.
#

PROTON = 1.00727646677


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



# common mz diff in LC-MS
# 
mass_signatures = [
(-82.95402353323, 'M-HCOOK+H[1+]'),
(-71.01382353323001, 'M-C3H4O2+H[1+]'),
(-66.98012353323, 'M-HCOONa+H[1+]'),
(-46.005479, 'H2O+CO'),
(-44.99812353323, 'M-HCOOH+H[1+]'),
(-42.98252353323, 'M-CO2+H[1+]'),
(-35.037114, 'NH3+H2O'),
(-35.01392353323, 'M-H4O2+H[1+]'),
(-26.98772353323, 'M-CO+H[1+]'),
(-19.01787646677, 'M-H2O-H[-]'),
(-18.010565, 'H2O'),
(-17.026549, 'NH3'),
(-17.00332353323, 'M-H2O+H[1+]'),
(-16.019223533229997, 'M-NH3+H[1+]'),
(-1.00727646677, '-H[-]'),
(-0.0038764667699999755, 'M(C13)-H[-]'),
(0.5017, 'double charged with C13'),
(0.98402, 'OH <-> NH2, e.g. de-amidiation, CHNO compounds'),
(0.98852353323, 'M(S34)-H[-]'),
(0.98992353323, 'M(Cl37)-H[-]'),
(1.0034, 'M(C13)'),
(1.00727646677, '+H[1+]'),
(1.34167646677, 'M(C13)+3H[3+]'),
(1.50897646677, 'M(C13)+2H[2+]'),
(1.97975, 'K+ <-> Cl-+2H2+, salt adduct'),
(1.99566, 'F <-> OH, halogen exchange with hydroxy group (typically -F + OH)'),
(1.9958, 'M(S34)'),
(1.9972, 'M(Cl37)'),
(2.01067646677, 'M(C13)+H[1+]'),
(2.014552, '2H'),
(2.01565, '± 2H, opening or forming of double bond'),
(3.00307646677, 'M(S34)+H[1+]'),
(3.00447646677, 'M(Cl37)+H[1+]'),
(3.021828, '3H'),
(4.9554, 'Na+<-> NH4+, salt adduct'),
(7.00467, 'F <-> CN, halogen exchange with cyano group'),
(8.96578, 'Cl <-> CN, halogen exchange with cyano group'),
(11.99827646677, 'M+H+Na[2+]'),
(13.97927, 'O <-> 2H, e.g. Oxidation follwed by H2O elimination'),
(13.99419, 'Cl-+2H2+ <-> Na+, salt adduct'),
(14.01565, "± CH2, alkane chains, waxes, fatty acids, methylation; or '-[C3H6ON] <-> -[C2H4ON], acrylamide versus iodoacetamide in cysteine alkylation (gels)"),
(14.987633533230001, 'M-H+O[-]'),
(15.97394, 'Na+<-> K+, salt adduct'),
(15.97716, 'S <-> O, sulfur compounds'),
(15.99492, '± O, e.g. oxidation/reduction'),
(17.02655, '± NH3, ammonium adduct/neutral ammonium loss; or NH4+<-> H+, salt adduct'),
(17.96611, 'Cl <-> OH, halogen exchange with hydroxy group (typically -Cl + OH)'),
(17.99058, 'F <-> H, halogen exchange'),
(18.01057, '± H2O, water addition/loss'),
(18.033823, 'NH4+'),
(18.99458, 'CN <-> COOH, nitrile compounds'),
(19.01787646677, 'M+H2O+H[1+]'),
(20.92933, 'K+<-> NH4+, salt adduct'),
(20.97474706646, 'M+Na-2H[-]'),
(21.98194, 'Na+<-> H+, salt adduct'),
(22.989218, 'Na'),
(23.996494, 'H+Na'),
(24.99525, 'CN <-> H, nitrile compounds'),
(25.00377, '2H+Na'),
(27.0109, '± HCN, nitrile compounds'),
(27.99492, '± CO'),
(28.00615, '- 2N, nitrogen loss, e.g. azido compounds (N2)'),
(28.0313, '± C2H4, natural alkane chains such as fatty acids'),
(29.97418, 'NO2 <-> NH2, nitro compounds'),
(29.99799, '-NO, nitroso compounds'),
(31.97207, '± S, sulfur compounds'),
(31.98983, '± 2O, oxygen loss'),
(32.026215, 'MeOH'),
(33.02146, '-NH2OH, loss from hydroxamic acids'),
(33.96103, 'Cl <-> H, halogen exchange'),
(33.98772, '± H2S, sulfur compounds'),
(34.969402, 'Cl-'),
(35.9767, 'Cl-+2H2+ <-> H+, salt adduct'),
(36.94824706646, 'M+K-2H[-]'),
(36.9664, 'M+Cl37[-]'),
(37.955882, 'K-H'),
(37.98916, '2CN <-> 2COOH, nitrile compounds'),
(38.963158, 'K'),
(39.970434, 'H+K'),
(39.99248, 'Na+ <-> -H2O+H'),
(40.01926853323, 'M+ACN-H[-]'),
(40.0313, '"(C3H6O - H2O), acetone condensation after dehydration"'),
(40.97771, '2H+K'),
(41.026549, 'Acetonitrile'),
(42.01057, '± COCH2'),
(42.04695, '± C3H6, propylation'),
(43.00581, '± CONH (wrong calc. in ref.)'),
(43.94948, '79Br <-> Cl, halogen exchange (typically -Br +Cl)'),
(43.96389, '2Na-2H'),
(43.98983, '± CO2'),
(44.998201, 'COOH-'),
(45.94744, '81Br <-> Cl, halogen exchange (typically -Br +Cl)'),
(45.978436, '2Na'),
(46.00548, '± CO+H2O (carboxylic acid)'),
(46.985712, 'H+2Na'),
(47.96699, '± SO, sulfur compounds'),
(52.91311, 'Fe3+ - 4H- <-> H-, salt adduct (probably only for polymers)'),
(52.91526, '79Br <-> CN, halogen exchange with cyano group'),
(54.91322, '81Br <-> CN, halogen exchange with cyano group'),
(56.0626, '± C4H8, butylation'),
(57.958622, 'NaCl'),
(58.00548, '± CO2CH2'),
(58.04187, '"+C3H6O, acetone condensation"'),
(58.96587646677, 'M+NaCl[1+]'),
(59.013295, 'M+CH3COO[-]'),
(59.03711, '± C2H5NO, N-acetyl or cut from side chain of Asn in peptides'),
(59.93783, 'Na+ + K+ - H+ <-> H+, salt adduct'),
(61.9156, '79Br <-> OH, halogen exchange with hydroxy group (typically -Br + OH)'),
(61.952376, 'Na+K'),
(62.959652, 'H+Na+K'),
(63.91355, '81Br <-> OH, halogen exchange with hydroxy group (typically -Br + OH)'),
(63.9619, '± SO2, sulfur compounds'),
(63.99829, '-CH3SOH, specific loss from singly oxidized methionine residues'),
(65.945835, '3Na-3H'),
(67.987424, 'NaCOOH'),
(68.967654, '3Na'),
(68.99467646676999, 'M+HCOONa[1+]'),
(77.91051, '79Br <-> H, halogen exchange'),
(77.926316, '2K'),
(78.9183, 'M+Br[-]'),
(78.933592, 'H+2K'),
(79.90847, '81Br <-> H, halogen exchange'),
(79.95682, '± SO3, sulfur compounds'),
(80.9163, 'M+Br81[-]'),
(83.961361, 'KCOOH'),
(84.941594, '2Na+K'),
(84.96857646676999, 'M+HCOOK[1+]'),
(89.969369, 'NaCOOH +Na-H'),
(91.93562, 'I <-> Cl, halogen exchange (typically -I +Cl)'),
(97.96738, '± H2SO4, sulfur compounds'),
(97.9769, '± H3PO4, phosphorous compounds'),
(100.9014, 'I <-> CN, halogen exchange with cyano group'),
(100.915534, 'Na+2K'),
(109.90173, 'I <-> OH, halogen exchange with hydroxy group (typically -I + OH)'),
(113.992903, 'HCOOH + NaCOOH'),
(115.917245, '2X NaCl'),
(116.889474, '3K'),
(125.89665, 'I <-> H, halogen exchange'),
(125.946046, 'NaCl + NaCOOH'),
(135.974848, '2X NaCOOH'),
(146.05791, '± [Deoxy-Hexose-H2O, C6O4H10], e.g. Fucose'),
(162.05283, '± [Hexose-H2O, C6O5H10], e.g. Glucose, Galactose, Mannose, Fructose'),
(164.06848, '± [Deoxy-Hexose-H2O, C6O5H12], e.g. Fucose'),
(176.03209, '± [Glucuronic acid-H2O, C6O6H8]'),
(180.06339, '± [Hexose, C6O6H12], e.g. Glucose, Galactose, Mannose, Fructose'),
(194.04266, '± [Glucuronic acid, C6O7H10]'),
(203.07937, '± [HexNAc-H2O, C8O5NH13], N-acetylhexoseamine'),
(206.05791, '+ sinapic acid - H2O, (C11H12O5-H2O), MALDI-matrix for proteins, possible photo- or condensation adduct'),
(221.08994, '± [HexNAc, C8O6NH15], N-acetylhexoseamine'),
(289.07324, '± [Glutathione-H2O, C10O5N3SH15]'),
(291.09542, '± [Neu5Ac-H2O, C11O8NH17], N-acetylneuraminic acid, sialic acid, NANA'),
(305.06816, '± [Glutathione+O-H2O, C10O6N3SH15]'),
(307.08381, '± [Glutathione, C10O6N3SH17]'),
(309.10598, '± [Neu5Ac, C11O9NH19], N-acetylneuraminic acid, sialic acid, NANA'),
(324.10565, '± [Sucrose-H2O, C12O10H20]'),
(342.11622, '± [Sucrose, C12O11H22]'),
(484.17101, '+Tetraphenyl-tetramethyl-trisiloxane, (C28H32O2Si3), ( from silicon-based diffusion pump oil (DC7040, Dow corning)'),
(546.18666, '+Pentaphenyl-trimethyl-trisiloxane, (C33H34O2Si3), ( from silicon-based diffusion pump oil (DC7050, Dow corning)')
]

