from itertools import product
import numpy as np
import pandas as pd
import os
import sys

in_file_name = sys.argv[1]
fastq_file = open(in_file_name,'r')
seq_file = open('trial_seq.fq', 'w+')

lines = 100
# creating sequence file out of fastq file 

count = 1

for line in fastq_file:
  if (count%4 == 2):
    count += 1
    seq_file.write(line)
  else:
    count += 1

seq_file.close()

k = 5
nucleotides = 'ATGC'
k_mers = [''.join(c) for c in product(nucleotides, repeat=k)]

occurence = { }
for i in range(len(k_mers)):

	k_mer = k_mers[i]
	file_name_1 = str(i+1) + ".txt"

	cmd = "grep -n %s 'trial_seq.fq' | cut -d : -f 1 > %s"%(k_mer, file_name_1)
	os.system(cmd)

	file = open(file_name_1,'r')

	line_numbers = file.readlines()
	total_lines = len(line_numbers)

	k = 0
	occurence[k_mer] = []
	for j in range(lines):
		if (k < total_lines):
			if (str(j+1) == line_numbers[k].rstrip()):
				occurence[k_mer].append(1)
				k += 1 
			else:
				occurence[k_mer].append(0)
		else:
			occurence[k_mer].append(0)

	file.close()

matrix_file = "matrix.txt"

data_frame = pd.DataFrame(occurence)

data_frame.to_csv(matrix_file, sep='\t',index=False)

os.system("rm [0-9]*.txt")