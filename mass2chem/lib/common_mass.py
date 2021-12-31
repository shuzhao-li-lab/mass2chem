'''
# atomic masses adopted from the National Institute of Standarts and Technology (NIST)
#   (http://physics.nist.gov/PhysRefData/Compositions/index.html)

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
'''

PROTON = 1.00727646677
electron = 0.000549

# common mz diff in LC-MS
# Electron is not considered in formula calculation, but may affect mass calculation.
# The dictionary in each entry is used for formula calculation, not mass.

# To-dos:
# to add 15N, 18O etc;
# to separate pos and neg ionization, and biological reactions if possible

mass_signatures = [
(-82.95402353323, 'M-HCOOK+H[1+]' , {'C':-1,'O':-2, 'K':-1, }),
(-71.01382353323001, 'M-C3H4O2+H[1+]', {'C':-3,'H':-3, 'O':-2,}),
(-66.98012353323, 'M-HCOONa+H[1+]', {'C':-1 , 'O':-2,'Na':-1 ,}),
(-46.005479, 'H2O+CO',{'C':-1 , 'H':-2, 'O':-2 }),
(-44.99812353323, 'M-HCOOH+H[1+]',{'C': -1, 'H':-1, 'O':-2, }),
(-42.98252353323, 'M-CO2+H[1+]',{'C':-1 , 'H':1, 'O':-2,}),
(-35.037114, 'NH3+H2O', {'H':-5, 'O':-1, 'N':-1} ),
(-35.01392353323, 'M-H4O2+H[1+]',{'H':-3, 'O':-2,}),
(-26.98772353323, 'M-CO+H[1+]',{'C':-1 , 'H':1, 'O':-1 ,}),
(-19.01787646677, 'M-H2O-H[-]', {'H':-3, 'O':-1 ,}),
(-18.010565, 'H2O',{'H':-2, 'O':-1}),
(-17.026549, 'NH3',{'N':-1 , 'H':-3}),
(-17.00332353323, 'M-H2O+H[1+]',{'H':-1,'O':-1 ,}),
(-16.019223533229997, 'M-NH3+H[1+]', {'N': -1, 'H':-2 ,}),
(-1.00727646677, '-H[-]',{'H':-1,}),
# 
(-0.0038764667699999755, 'M(C13)-H[-]',{ 'C':-1, '(C13)':1, 'H':-1}),
(0.5017, 'double charged with C13',{'C':-1 , '(C13)':1 }),
(0.98402, 'OH <-> NH2, e.g. de-amidiation, CHNO compounds', {'O':1 , 'H':-1, 'N':-1}),
(0.98852353323, 'M(S34)-H[-]', {'S':-1, '(S34)':1, 'H':-1}),
(0.98992353323, 'M(Cl37)-H[-]',{ 'Cl':-1,'(Cl37)':1, 'H':-1 }),
(1.0034, 'M(C13)',{'C':-1 , '(C13)':1 }),
(1.00727646677, '+H[1+]',{'H':1 }),
(1.34167646677, 'M(C13)+3H[3+]',{'C':-1 , '(C13)':1 , 'H':3}),
(1.50897646677, 'M(C13)+2H[2+]',{'C':-1 , '(C13)':1 , 'H':2}),
(1.97975, 'K+ <-> Cl-+2H2+, salt adduct', {'K':-1, 'Cl':1, 'H':2}),
(1.99566, 'F <-> OH, halogen exchange with hydroxy group (typically -F + OH)', {'F':1 , 'O':-1, 'H':-1 }),
(1.9958, 'M(S34)', {'S':-1, '(S34)':1,}),
(1.9972, 'M(Cl37)',{'Cl':-1,'(Cl37)':1}),
(2.01067646677, 'M(C13)+H[1+]',{'C':-1 , '(C13)':1 , 'H':1}),
(2.014552, '2H', {'H':2 }),
(2.01565, '± 2H, opening or forming of double bond', {'H':2}),
(3.00307646677, 'M(S34)+H[1+]',{'S':-1,'(S34)':1 , 'H':1 }),
(3.00447646677, 'M(Cl37)+H[1+]', {'Cl':-1, '(Cl37)':1, 'H':1}),
(3.021828, '3H',{'H':3}),
(4.9554, 'Na+<-> NH4+, salt adduct',{'Na':1 , 'N':-1 , 'H':-4}),
(7.00467, 'F <-> CN, halogen exchange with cyano group',{'F':1 , 'C':-1, 'N':-1 }),
(8.96578, 'Cl <-> CN, halogen exchange with cyano group',{'Cl':1 , 'C':-1, 'N':-1 }),
(11.99827646677, 'M+H+Na[2+]',{'H':1 , 'Na':1}),
#
(13.97927, 'O <-> 2H, e.g. Oxidation follwed by H2O elimination', {'H':-2, 'O':1 }),
(13.99419, 'Cl-+2H2+ <-> Na+, salt adduct',{'Na':-1, 'Cl':1, 'H':2}),
(14.01565, "± CH2, alkane chains, waxes, fatty acids, methylation; or '-[C3H6ON] <-> -[C2H4ON], acrylamide versus iodoacetamide in cysteine alkylation (gels)",{'C':1, 'H':2 }),
(14.987633533230001, 'M-H+O[-]',{'H':-1, 'O':1 }),
#
(15.97394, 'Na+<-> K+, salt adduct', {'Na':-1, 'K':1 }),
(15.97716, 'S <-> O, sulfur compounds',{'S':1, 'O':-1 }),
(15.99492, '± O, e.g. oxidation/reduction', {'O':1 }),
(17.02655, '± NH3, ammonium adduct/neutral ammonium loss; or NH4+<-> H+, salt adduct',{'N':1, 'H':3 }),
(17.96611, 'Cl <-> OH, halogen exchange with hydroxy group (typically -Cl + OH)',{'Cl':1, 'O':-1 , 'H':-1}),
(17.99058, 'F <-> H, halogen exchange',{'F':1, 'H':-1 }),
(18.01057, '± H2O, water addition/loss',{'H':2, 'O':1 } ),
(18.033823, 'NH4+',{'N':1, 'H':4 }),
(18.99458, 'CN <-> COOH, nitrile compounds', { 'N':-1 , 'H':1, 'O':2}),
(19.01787646677, 'M+H2O+H[1+]',{ 'H': 3, 'O':1}),
(20.92933, 'K+<-> NH4+, salt adduct', { 'K':1 , 'N':-1, 'H':-4}),
(20.97474706646, 'M+Na-2H[-]', { 'Na':1 , 'H':-2 }),
(21.98194, 'Na+<-> H+, salt adduct', { 'Na':1 , 'H':-1}),
(22.989218, 'Na',{ 'Na':1 }),
(23.996494, 'H+Na',{ 'Na':1 , 'H':1}),
(24.99525, 'CN <-> H, nitrile compounds', { 'C':1,'N':1 , 'H':-1}),
(25.00377, '2H+Na',{ 'Na':1 , 'H':2}),
(27.0109, '± HCN, nitrile compounds',{ 'N':1 , 'H':1, 'C':1}),
(27.99492, '± CO', { 'C':1 , 'O':1}),
#
(28.00615, '- 2N, nitrogen loss, e.g. azido compounds (N2)',{ 'N':2}),
(28.0313, '± C2H4, natural alkane chains such as fatty acids',{ 'C':2 , 'H':4}),
(29.97418, 'NO2 <-> NH2, nitro compounds',{ 'H':-2, 'O':2}),
(29.99799, '-NO, nitroso compounds',{ 'N': -1, 'O':-1}),
(31.97207, '± S, sulfur compounds',{ 'S': 1}),
(31.98983, '± 2O, oxygen loss',{ 'O':2}),
(32.026215, 'MeOH',{ 'C':1 , 'H':4, 'O':1}),
(33.02146, '-NH2OH, loss from hydroxamic acids', { 'N':1 , 'H':3, 'O':1}),
(33.96103, 'Cl <-> H, halogen exchange',{ 'Cl': 1, 'H':-1}),
(33.98772, '± H2S, sulfur compounds',{ 'H': 2, 'S':1}),
(34.969402, 'Cl-',{ 'Cl': -1}),
# 
(35.9767, 'Cl-+2H2+ <-> H+, salt adduct',{'Cl': 1, 'H':1}),
(36.94824706646, 'M+K-2H[-]',{ 'K':1, 'H':-2}),
(36.9664, 'M+Cl37[-]', { '(Cl37)':1}),
(37.955882, 'K-H',{ 'K':1 , 'H':-1}),
(37.98916, '2CN <-> 2COOH, nitrile compounds',{'N':-2,'O':4 , 'H':2}),
(38.963158, 'K',{ 'K':1} ),
(39.970434, 'H+K',{ 'K':1 , 'H':1}),
(39.99248, 'Na+ <-> -H2O+H',{ 'Na':1, 'O':1, 'H':1}),
# Acetonitrile is C2H3N, 24 + 3 * 1.007825 + 14.003074 - 1.00727 = 40.019279
(40.01926853323, 'M+ACN-H[-]', { 'C':2, 'N':1,'H':2}),
(40.0313, '"(C3H6O - H2O), acetone condensation after dehydration"', { 'H':4 ,'C':3}),
(40.97771, '2H+K',{ 'H':2 ,'K':1}),
(41.026549, 'Acetonitrile', { 'C':2 ,'H':3, 'N':1}),
(42.01057, '± COCH2',{ 'C':2 ,'O':1, 'H':2} ),
(42.04695, '± C3H6, propylation', { 'C':3 ,'H':6}),
(43.00581, '± CONH (wrong calc. in ref.)',{ 'C':1 ,'O':1, 'N':1, 'H':1}),
(43.94948, '79Br <-> Cl, halogen exchange (typically -Br +Cl)', { 'Br':1 ,'Cl':-1}),
(43.96389, '2Na-2H', { 'Na':2 ,'H':-2}),
(43.98983, '± CO2', { 'C':1 ,'O':2}),
(44.998201, 'COOH-', { 'C':1 ,'O':2, 'H':1}),
(45.94744, '81Br <-> Cl, halogen exchange (typically -Br +Cl)',{'(Br81)':1 ,'Cl':-1}),
(45.978436, '2Na', { 'Na':2}),
(46.00548, '± CO+H2O (carboxylic acid)', { 'C':1 ,'O':2, 'H':2}),
(46.985712, 'H+2Na',{ 'H':1 ,'Na':2}),
(47.96699, '± SO, sulfur compounds',{ 'S':1 ,'O':1}),
# 55.934938 - 1.00727*3 = 52.913128
(52.91311, 'Fe3+ - 4H- <-> H-, salt adduct (probably only for polymers)',{'Fe':1 ,'H':-3 }),
(52.91526, '79Br <-> CN, halogen exchange with cyano group',{'Br':1 ,'C':-1, 'N':-1} ),
(54.91322, '81Br <-> CN, halogen exchange with cyano group', {'(Br81)':1 ,'C':-1, 'N':-1}),
(56.0626, '± C4H8, butylation', {'C':4 ,'H':8}),
(57.958622, 'NaCl', {'Na':1 ,'Cl':1}),
(58.00548, '± CO2CH2', {'C':2 ,'O':2,'H':2}),
(58.04187, '"+C3H6O, acetone condensation"',{'C':3 ,'H':6, 'O':1}),
(58.96587646677, 'M+NaCl[1+]', {'Na':1 ,'Cl':1}),
(59.013295, 'M+CH3COO[-]', {'C':2 ,'H':3, 'O':2}),
(59.03711, '± C2H5NO, N-acetyl or cut from side chain of Asn in peptides', {'C':2 ,'H':5, 'N':1, 'O':1}),
#?Na+K=61.95347 ; 61.95347 - 1.007825*2 = 59.93782
(59.93783, 'Na+ + K+ - H+ <-> H+, salt adduct',  {'K':1, 'Na':1, 'H':-2}),
(61.9156, '79Br <-> OH, halogen exchange with hydroxy group (typically -Br + OH)', {'Br':1,'O':-1, 'H':-1}),
(61.952376, 'Na+K', {'K':1, 'Na':1}),
(62.959652, 'H+Na+K', {'H':1,'Na':1, 'K':1}),
(63.91355, '81Br <-> OH, halogen exchange with hydroxy group (typically -Br + OH)', {'(Br81)':1,'O':-1, 'H':-1}),
(63.9619, '± SO2, sulfur compounds',{'S':1,'O':2}),
(-63.99829, '-CH3SOH, specific loss from singly oxidized methionine residues',{'C':-1,'H':-4, 'S':-1, 'O':-1}),
(65.945835, '3Na-3H',{'Na':3,'H':-3}),
(67.987424, 'NaCOOH', {'C':1,'O':2, 'Na':1, 'H':1}),
(68.967654, '3Na', {'Na':3}),
(68.99467646676999, 'M+HCOONa[1+]',{'C':1,'O':2, 'Na':2, 'H':1}),
(77.91051, '79Br <-> H, halogen exchange', {'Br':1,'H':-1}),
(77.926316, '2K', {'K':2}),
(78.9183, 'M+Br[-]', {'Br':1}),
(78.933592, 'H+2K', {'H':1,'K':2}),
(79.90847, '81Br <-> H, halogen exchange', {'(Br81)':1,'H':-1}),
(79.95682, '± SO3, sulfur compounds', {'O':3, 'S':1}),
(80.9163, 'M+Br81[-]', {'(Br81)':1}),
(83.961361, 'KCOOH', {'O':2,'C':1, 'K':1, 'H':1}),
(84.941594, '2Na+K', {'Na':2,'K':1}),
(84.96857646676999, 'M+HCOOK[1+]', {'C':1,'H':1, 'O':2, 'K':1}),
#
(89.969369, 'NaCOOH +Na-H',{'C':1,'O':2, 'Na':2}),
(91.93562, 'I <-> Cl, halogen exchange (typically -I +Cl)', {'I':1 ,  'Cl':-1}),
(97.96738, '± H2SO4, sulfur compounds', {'O':4,'H':2, 'S':1}),
(97.9769, '± H3PO4, phosphorous compounds', {'O':4,'H':3, 'P':1}),
(100.9014, 'I <-> CN, halogen exchange with cyano group', {'C':-1 ,  'N':-1, 'I':1}),
(100.915534, 'Na+2K', {'K':2, 'Na':1}),
(109.90173, 'I <-> OH, halogen exchange with hydroxy group (typically -I + OH)', {'I':1 ,  'H':-1, 'O':-1}),
(113.992903, 'HCOOH + NaCOOH', {'C':2 ,  'O':4,'H':3, 'Na':1}),
(115.917245, '2X NaCl', {'Cl':2, 'Na':2}),
(116.889474, '3K', {'K':3 }),
(125.89665, 'I <-> H, halogen exchange', {'I':1 ,  'H':-1}),
(125.946046, 'NaCl + NaCOOH',{'C':1 ,  'O':2,'H':1, 'Na':2, 'Cl':1}),
(135.974848, '2X NaCOOH', {'C':2 ,  'O':4,'H':2, 'Na':2}),
#
(146.05791, '± [Deoxy-Hexose-H2O, C6O4H10], e.g. Fucose', {'C':6 ,'O':4,'H':10}),
(162.05283, '± [Hexose-H2O, C6O5H10], e.g. Glucose, Galactose, Mannose, Fructose', {'C':6 ,'O':5,'H':10}),
(164.06848, '± [Deoxy-Hexose-H2O, C6O5H12], e.g. Fucose', {'C':6 , 'O':5,'H':12}),
(176.03209, '± [Glucuronic acid-H2O, C6O6H8]',{'C':6,  'O':6,'H':8}),
(180.06339, '± [Hexose, C6O6H12], e.g. Glucose, Galactose, Mannose, Fructose', {'C':6 ,  'O':6,'H':12}),
(194.04266, '± [Glucuronic acid, C6O7H10]', {'C':6 ,  'O':7,'H':10}),
(203.07937, '± [HexNAc-H2O, C8O5NH13], N-acetylhexoseamine', {'C':8 ,  'O':5,'N':1,'H':13}),
(206.05791, '+ sinapic acid - H2O, (C11H12O5-H2O), MALDI-matrix for proteins, possible photo- or condensation adduct', {'C':11 ,  'O':4,'H':10}),
(221.08994, '± [HexNAc, C8O6NH15], N-acetylhexoseamine', {'C':8,  'O':6, 'N':1, 'H':15}),
(289.07324, '± [Glutathione-H2O, C10O5N3SH15]',{'C':10 ,  'O':5, 'N':3, 'S':1,'H':15}),
(291.09542, '± [Neu5Ac-H2O, C11O8NH17], N-acetylneuraminic acid, sialic acid, NANA',{'C':11,  'O':8, 'N':1, 'H':17}),
(305.06816, '± [Glutathione+O-H2O, C10O6N3SH15]', {'C':10 ,  'O':6, 'N':3, 'S':1,'H':15}),
(307.08381, '± [Glutathione, C10O6N3SH17]',{'C':10 , 'O':6,'N':3, 'S':1, 'H':17}),
(309.10598, '± [Neu5Ac, C11O9NH19], N-acetylneuraminic acid, sialic acid, NANA',{'C':11 , 'H':19, 'O':9, 'N':1}),
(324.10565, '± [Sucrose-H2O, C12O10H20]',{'C':12 , 'H':20, 'O':10}),
(342.11622, '± [Sucrose, C12O11H22]',{'C':12 , 'H':22, 'O':11}),
(484.17101, '+Tetraphenyl-tetramethyl-trisiloxane, (C28H32O2Si3), ( from silicon-based diffusion pump oil (DC7040, Dow corning)',{'C':28 , 'H':32, 'O':2,'Si':3}),
(546.18666, '+Pentaphenyl-trimethyl-trisiloxane, (C33H34O2Si3), ( from silicon-based diffusion pump oil (DC7050, Dow corning)', {'C':33 , 'H':34, 'O':2, 'Si':3}),
]
