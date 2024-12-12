from matchms.importing import load_from_msp
from matchms.similarity import CosineHungarian
from matchms.exporting import save_as_msp
from matchms.filtering import default_filters, reduce_to_number_of_peaks

import os
import numpy as np
import tqdm
import sys

path = sys.argv[1]
min_peaks = int(sys.argv[2])
max_peaks = int(sys.argv[3])
spectral_registry = {}
total = 0
for x in tqdm.tqdm(load_from_msp(path)):
    try:
        inchikey = x.metadata_dict()['inchikey']
        x = reduce_to_number_of_peaks(default_filters(x), min_peaks, max_peaks)
        if inchikey:
            if inchikey not in spectral_registry:
                spectral_registry[inchikey] = []
            spectral_registry[inchikey].append(x)
            total += 1
    except:
        pass

print("Unfiltered MS2 Size: ", str(len(spectral_registry)), " Compounds with ", str(total), " MS2 Spectra")

selected_spectra = {}
comparator = CosineHungarian()
for inchikey, spectrum_set in spectral_registry.items():
    spec_matrix = np.array([[None for _ in spectrum_set] for _ in spectrum_set])
    for i, spectrum_1 in enumerate(spectrum_set):
        for j, spectrum_2 in enumerate(spectrum_set):
            if spec_matrix[i][j] is None:
                if i == j:
                    spec_matrix[i][j] = spec_matrix[j][i] = 0
                else:
                    cosine_result = comparator.pair(spectrum_1, spectrum_2)
                    cosine_meta_score = cosine_result['matches'] * cosine_result['score']
                    spec_matrix[i][j] = spec_matrix[j][i] = cosine_meta_score
    totals = np.sum(spec_matrix, axis=1)
    selected_spectra[inchikey] = spectrum_set[np.argmax(totals)]

print("Filtered MS2 Size: ", str(len(spectral_registry)), " Compounds with ", str(len(list(selected_spectra.values()))), " MS2 Spectra")

save_as_msp(list(selected_spectra.values()), "filtered_" + os.path.basename(path))



