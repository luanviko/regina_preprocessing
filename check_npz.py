import numpy as np, sys

# Load up numpy arrays in file
print(f"Loading npz archive: {sys.argv[1]}")
loaded_arrays = np.load(sys.argv[1])

# Print keys/files
print("\nFiles (keys) found in the archive: ")
print(loaded_arrays.files)

# Print dimensions
print("\nShape of the numpy arrays of each file (key):")
for i in range(0, len(loaded_arrays.files)):
    print(f"{loaded_arrays.files[i]}: ", loaded_arrays[loaded_arrays.files[i]].shape)