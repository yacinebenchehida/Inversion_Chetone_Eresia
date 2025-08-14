#!/usr/bin/python
from Bio import SeqIO
import sys
import os

#usage: python -W ignore *fasta

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")

for record in input_seq_iterator:
	print(record.id + "\t" + str(len(record.seq)))
