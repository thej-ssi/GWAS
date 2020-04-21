#!/usr/bin/env python3

# usage: python3 batch_transcribe.py pan_genome_sequences pan_genome_protein_sequences

input_dir = sys.argv[1]
output_dir = sys.argv[2]
script_path = sys.arv[3]

import os
import sys

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

files = os.listdir(input_dir)

for file in files:
	src = os.path.join(input_dir,file)
	dst = os.path.join(output_dir,file)
	cmd = 'python3 ' + script_path + ' ' + src + ' > ' + dst
	print(cmd)
	os.system(cmd)