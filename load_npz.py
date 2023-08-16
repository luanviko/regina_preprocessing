import numpy as np, sys
import ROOT

loaded_arrays = np.load(sys.argv[1])

print(loaded_arrays.files)
for array in loaded_arrays.files:
    print(loaded_arrays[array].shape)