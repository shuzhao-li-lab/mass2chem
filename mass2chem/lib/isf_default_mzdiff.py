# The isotope, adduct and ISF lists used in Chi et al

# pos ionization
isotope_search_patterns_pos = [ (1.003355, '13C/12C', (0, 0.8)),
                            (2.00671, '13C/12C*2', (0, 0.8)),
                            (3.010065, '13C/12C*3', (0, 0.8)),
                            # (3.9948, '44Ca/40Ca', (0, 0.1)), # 2%
                            (1.9970, '37Cl/35Cl', (0.1, 0.8)), # 24.24%
                            ]

isotope_search_patterns_neg = [ (1.003355, '13C/12C', (0, 0.8)),
                            (2.00671, '13C/12C*2', (0, 0.8)),
                            (3.010065, '13C/12C*3', (0, 0.8)),
                            (1.9970, '37Cl/35Cl', (0.1, 0.8)), # 24.24%
                            (1.9958, '32S/34S', (0, 0.1)), # 4%
                            ]

adduct_search_patterns_pos = [  # initial patterns are relative to M+H+
                            (21.98194, 'Na/H'),
                            (41.026549, 'ACN'),     # Acetonitrile
                            (67.987424, 'NaCOOH'),
                            (37.955882, 'K/H'),
                            (32.026215, 'CH3OH')
                            ]
adduct_search_patterns_neg = [  
                            (21.98194, 'Na/H'), 
                            (67.987424, 'NaCOOH'),
                            (82.0030, 'NaCH2COOH'),
                            # (1.99566, 'F <-> OH'), 
                            (41.026549, 'ACN'),
                            (37.955882, 'K/H'),
                            ]
extended_adducts = [  # excluding neutral loss here; include as a step after khipu
                            (1.0078, 'H'),
                            (17.02655, 'NH3'),
                            (18.0106, 'H2O'),      # easy to confuse with bio reactions
                            (18.033823, 'NH4'),
                            (27.01089904, 'HCN'),
                            (27.99492, 'CO'),
                            (32.026215, 'CH3OH'),
                            (35.9767, 'HCl'),
                            (37.94694, 'Ca/H2'),
                            (43.96389, 'Na2/H2'),
                            (46.00548, 'CO2H2'),
                            (67.987424, 'NaCOOH'),
                            (83.961361, 'KCOOH'),
                            (97.96737927, 'H2SO4'),
                            (97.97689507, 'H3PO4'),
]

neutral_loss_patterns_pos = [(14.015649,
                'addition of acetic acid and loss of CO2. Reaction: (+C2H2O2) and (-CO2)'),
                (18.010565, 'water'),
                (2.01565, '± 2H, opening or forming of double bond'),
                (44.0262, 'hydroxyethylation'),
                (28.0313, '± C2H4, natural alkane chains such as fatty acids'),
                (15.9949, 'oxidation'),
                (17.0265, 'addition of ammonia. Reaction: (+NH3)'),
                (26.01565, 'acetylation and loss of oxygen. Reaction: (+C2H2O) and (-O)'),
                (27.9949, 'addition of CO. Reaction: (+CO)'),
                (12.0, 'methylation and reduction'),
                (42.010564, 'malonylation and loss of CO2. Reaction: (+C3H2O3) and (-CO2)'),
                (67.987424, 'NaCOOH'),
                (13.979264,
                'nitrification and loss of oxygen. Reaction: (NH2 -> NO2) and (-O)'),
                (24.0, 'acetylation and loss of water. Reaction: (+C2H2O) and (-H2O)'),
                (16.0313, 'Methylation + reduction'),
                (42.04695, '± C3H6, propylation'),
                (46.005305, 'formic acid adduct'),
                (88.052429, 'butanoic acid'),
                (41.026549, 'Acetonitrile'),
                (30.04695, 'addition of C2H4 and hydrogenation. Reaction: (+C2H4) and (+H2)'),
                (35.037114,
                'addition of water and addition of ammonia. Reaction: (+H2O) and (+NH3)')]

neutral_loss_patterns_neg = [
                (67.9873, '\'NaCOOH\', "{\'C\': 1, \'O\': 2, \'Na\': 1, \'H\': 1}"'),
                (14.0155,
                '"± CH2, alkane chains, waxes, fatty acids, methylation; or \'-C3H6ON <-> -C2H4ON, acrylamide versus iodoacetamide in cysteine alkylation (gels)", "{\'C\': 1, \'H\': 2}"'),
                (2.0156, '\'± 2H, opening or forming of double bond\', "{\'H\': 2}"'),
                (82.0029, ''),
                (15.9948, '\'± O, e.g. oxidation/reduction\', "{\'O\': 1}"'),
                (43.9897, '\'± CO2\', "{\'C\': 1, \'O\': 2}"'),
                (18.0104, '\'H2O\', "{\'H\': -2, \'O\': -1}"'),
                (26.0155, 'C2H2'),
                (46.0053, '\'H2O+CO\', "{\'C\': -1, \'H\': -2, \'O\': -2}"'),
                (28.0312,
                '\'± C2H4, natural alkane chains such as fatty acids\', "{\'C\': 2, \'H\': 4}"'),
                (27.9948, '\'± CO\', "{\'C\': 1, \'O\': 1}"'),
                (1.9957,
                '\'F <-> OH, halogen exchange with hydroxy group (typically -F + OH)\', "{\'F\': 1, \'O\': -1, \'H\': -1}"'),
                (42.0104, '\'± COCH2\', "{\'C\': 2, \'O\': 1, \'H\': 2}"')
]