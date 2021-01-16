#!/usr/bin/env python3

# BASE FILE: provided by professor

import sys
import logging

if len(sys.argv) != 1+2:
	print(f"Syntax: {sys.argv[0]} FILE.sam GENOME_LENGTH")
	print()
	print("Description:")
	print("  Returns a track with the sequence coverage in .wig format.")
	print()
	print("Example:")
	print(f"  $ {sys.argv[0]} lact.sam 3079196")
	print(f"  $ {sys.argv[0]} lact.sam 3079196 > lact_physicalcov.wig")
	exit(1)

# Get args
input_file = sys.argv[1]
genome_length = int(sys.argv[2])

# Initialize genome_change variable as a list constituted by 0 
# with length = genome_length
genome_change = [0]*genome_length

# Open sam file
sam_file = open(input_file)

# Compute physical coverage: the physical coverage of a genomic 
# position is given by the number of fragments laying on that 
# position rather than the number of read
for line in sam_file:
	if line[0] == '@':
		continue

	fields = line.split("\t")   # makes a list of individual tab-separated fields
	if ((int(fields[1]) & 4) == 0):   # flag unset indicates that segment maps 
		starting_mate_position = int(fields[3])
		mate_length = 100

		# increment start position by one
		genome_change[starting_mate_position] += 1
		# decrement end position by one
		genome_change[starting_mate_position + mate_length] -= 1

# Close the sam file
sam_file.close()
	
# Print genomic profile as a wiggle file
print("fixedStep chrom=genome start=1 step=1 span=1")

current_coverage = 0

# Cicle over all positions of the genome
for position in range(genome_length):
	current_coverage += genome_change[position]
	print(genome_change[position])
