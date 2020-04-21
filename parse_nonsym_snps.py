#!/usr/bin/env python3

import os
import sys

#nucl_dir = sys.argv[1]
prot_dir = sys.argv[1]
output_dir = sys.argv[2]

if not os.path.exists(output_dir):
	os.makedirs(output_dir)


def read_fasta(input_file):
	headers = []
	seqs = []
	with open(input_file) as f:
		for line in f:
			line = line.rstrip('\n')
			if line[0] == '>':
				headers.append(line[1:])
				seqs.append('')
			else:
				seqs[-1] += line
	return(headers,seqs)


def identify_snps(headers,seqs,file_name):
	seq1 = seqs[0]
	ref_header = headers[0]
	return_line = ''
	if seq1.startswith('M') and seq1.endswith('*'):
		for i in range(len(headers)):
			for n in range(min(len(seqs[i]),len(seq1))):
				if seq1[n] != seqs[i][n]:
					prot_change = seq1[n]+str(n+1)+seqs[i][n]
					return_line += file_name+'\t'+headers[i]+'\t'+prot_change+'\t'+str(n+1)+'\t'+seq1[n]+'\t'+seqs[i][n] + '\n'
	return(return_line)


files = os.listdir(prot_dir)
for file in files:
	src = os.path.join(prot_dir,file)
	headers,seqs = read_fasta(src)
	snp_table = identify_snps(headers,seqs,file)
	dst = os.path.join(output_dir,file+'.txt')
	o = open(dst,'w')
	o.write(snp_table)
	o.close()
