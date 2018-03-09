from itertools import product
import numpy as np
import sys

in_file_name = sys.argv[1]

fastq_file = open(in_file_name,'r')
seq_file = open('trial_seq.fq', 'w+')

# creating sequence file out of fastq file 

count = 1

for line in fastq_file:
  if (count%4 == 2):
    count += 1
    seq_file.write(line)
  else:
    count += 1

# creating all possible 4-mers

k = 5
nucleotides = 'ATGC'
k_mers = [''.join(c) for c in product(nucleotides, repeat=k)]

# Knuth Morris Pratt (KMP) Algorithm for exact string matching

# pre processig the k-mer

def prefix_array(pattern):
  j = 0
  k = 1
  skips = [0]
  while k < len(pattern):
    if j == 0 and pattern[j] != pattern[k]:
      skips.append(j)
      k += 1
      
    elif pattern[j] == pattern[k]:
      skips.append(j+1)
      j += 1
      k += 1
      
    elif pattern[j] != pattern[k]:
      j = skips[j-1]

  return skips

# KMP String matching function

def KMP(text, pattern):
    i = 0
    j = 0
    result = 0
    skips = prefix_array(pattern)
    while i < len(text) and j < len(pattern):
        
        if text[i] == pattern[j] and j == len(pattern)-1:
            result = 1
            if j != 0 and i != 0:
                j = skips[j-1]
            elif j==0:
                i += 1
            break
            
        elif text[i] == pattern[j]:
            i += 1
            j += 1

        elif text[i] != pattern[j] and j != 0:
            j = skips[j-1]

        elif text[i] != pattern[j] and j == 0:
            i += 1

    return result

# creating bitwise matrix


seq_file = open('trial_seq.fq', 'r')

seq_number = 0

output_file = open('matrix_2.txt', 'w+')

seq_file = open('trial_seq.fq', 'r')

for line in seq_file:
  for i in range(len(k_mers)):
    k_mer = k_mers[i]
    bit = KMP(line.rstrip(), k_mer)
    output_file.write("%i\t" % bit)
  output_file.write("\n")

