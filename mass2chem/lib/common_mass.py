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
