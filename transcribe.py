#!/usr/bin/env python3

import sys

d = {}
with open(sys.argv[1]) as f:
	for line in f:
		line = line.rstrip('\n')
		if line[0] == '>':
			header = line
			d[header] = ''
		else:
			d[header] += line
f.close()

for header in d:
	count = 0
	gapcount = 0
	peptide = ''
	seq = ''
	original_seq = d[header]
	seq_without_start_gaps = original_seq.lstrip('-')
	start_gap_length = len(original_seq)-len(seq_without_start_gaps)
	peptide += '-'*int(start_gap_length/3)
	modulo = start_gap_length%3
	if modulo == 1:
		start_remove = 2
		peptide += '-'
	elif modulo == 2:
		start_remove = 1
		peptide += '-'
	elif modulo == 0:
		start_remove=0
	seq_without_start_gaps = seq_without_start_gaps[start_remove:]
	for nt in seq_without_start_gaps:
		if nt != '-':
			seq += nt
			count += 1
		else:
			gapcount +=1
			if gapcount >= 3:
				peptide+='-'
				gapcount = 0
		if count >= 3:
			seq = seq.upper()
			count = 0
			if gapcount != 0:
				gapcount = 0
			if seq == 'TTT' or seq == 'TTC': peptide += 'F'
			elif seq == 'TTA' or seq == 'TTG': peptide += 'L'
			elif seq[:2] == 'CT': peptide += 'L'
			elif seq == 'ATT' or seq == 'ATC' or seq == 'ATA': peptide += 'I'
			elif seq == 'ATG': peptide += 'M'
			elif seq[:2] == 'GT': peptide += 'V'
			elif seq[:2] == 'TC': peptide += 'S'
			elif seq[:2] == 'CC': peptide += 'P'
			elif seq[:2] == 'AC': peptide += 'T'
			elif seq[:2] == 'GC': peptide += 'A'
			elif seq == 'TAT' or seq == 'TAC': peptide += 'Y'
			elif seq == 'TAA' or seq == 'TAG': peptide += '*'
			elif seq == 'CAT' or seq == 'CAC': peptide += 'H'
			elif seq == 'CAA' or seq == 'CAG': peptide += 'Q'
			elif seq == 'AAT' or seq == 'AAC': peptide += 'N'
			elif seq == 'AAA' or seq == 'AAG': peptide += 'K'
			elif seq == 'GAT' or seq == 'GAC': peptide += 'D'
			elif seq == 'GAA' or seq == 'GAG': peptide += 'E'
			elif seq == 'TGT' or seq == 'TGC': peptide += 'C'
			elif seq == 'TGA': peptide += '*'
			elif seq == 'TGG': peptide += 'W'
			elif seq[:2] == 'CG': peptide += 'R'
			elif seq == 'AGT' or seq == 'AGC': peptide += 'S'
			elif seq == 'AGA' or seq == 'AGG': peptide += 'N'
			elif seq[:2] == 'GG': peptide += 'G'
			else: peptide += 'X'

			seq = ''

	print(header)
	n=0
	while n < len(peptide):
		print(peptide[n:n+60])
		n+=60
