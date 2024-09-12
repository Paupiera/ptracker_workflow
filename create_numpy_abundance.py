import numpy as np
import sys

#file_to_read = "./abundances_sample_4.tsv.bck"
#out_file = "abundance.npz"
file_to_read = sys.argv[1:-1]
out_file = sys.argv[-1]
print(file_to_read, out_file)

if len(file_to_read) < 1:
    sys.exit()

# Read in the abundances to a list
array = []
sample_names = []

# Ssample_8CNODE_3078_length_2004_cov_3.886095    3.792415 

for i, file in enumerate(file_to_read):
    array.append([]) 
    sample_names.append(file)

    with open(file, "r") as f:
        for line in f.readlines():
            line = line.split()
            name = line[0]
            length = name.split("_")[4]
            try:
                print(line)
                int(length)
            except ValueError:
                raise ValueError(length)
            
            if int(length) >= 2000:
                abundance = float(line[1])
                array[i].append(abundance)
            
np_array = np.asanyarray(array, dtype="float32").T
print(np_array, sample_names)
#print(type(np_array[0]))
#print(np_array.dtype)
print(len(np_array.T[0]))

# Write to file as npz
np.savez(out_file, matrix=np_array, samplenames=sample_names, minid=np.array([0.0]), refhash=np.array([None]))

#a = np.load("abundance.npz")
#print(a.files)




