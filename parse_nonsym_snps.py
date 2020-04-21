#!/usr/bin/env python3

import os
import sys

nucl_dir = sys.argv[1]
prot_dir = sys.argv[2]
output_dir = sys.argv[3]

if not os.path.exists(output_dir):
	os.makedirs(output_dir)


def read_fasta(input_file):
	headers = []
	seqs = []
	with open(input_file) as f:
		for line in f:
			line = line.rstrip('\n')
			if line[0] == '>':
				headers.append[line[1:]]
				seqs.append['']
			else:
				seqs[-1] += line
	return(headers,seqs)


def identify_snps(headers,seqs):
	seq1 = seqs[0]
	ref_header = headers[0]
	prot_change_list = [0]*length(seqs)
	unique_prot_changes = []
	if seq1.startswith('M') and seq1.endswith('*'):
		for i in range(length(headers)):
			for n in range(length(seq1)):
				if seq1[n] != seqs[i][n]:
					prot_change = seq1[n]+str(n+1)+seqs[i][n]
					prot_change[i].append(prot_change)
					if prot_change not in uniq_prot_changes:
						uniq_prot_changes.append(prot_changes)
	return(uniq_prot_changes,prot_changes)

headers,seqs = read_fasta(sys.argv[1])
uniq,prot = identify_snps(headers,seqs)
