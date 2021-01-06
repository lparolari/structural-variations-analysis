#!/usr/bin/env python3

import sys
import logging

if len(sys.argv) != 1+2:
	print(f"Syntax: {sys.argv[0]} FILE.sam GENOME_LENGTH")
	print()
	print("Description:")
	print("  Returns a track with the physical coverage in .wig format.")
	print()
	print("  The physical coverage of a genomic position is given by the")
	print("  number of fragments laying on that position rather than the")
	print("  number of read.")
	print()
	print("Example:")
	print(f"  $ {sys.argv[0]} lact.sam 3079196")
	print(f"  $ {sys.argv[0]} lact.sam 3079196 > lact.wig")
	exit(1)

# Get args
input_file = sys.argv[1]
genome_length = int(sys.argv[2])

# Initialize genome_change variable as a list constituted by 0 
# with length = genome_length
genome_change = [0]*genome_length

# Open sam file
sam_file = open(input_file)

# Compute sequence coverage
for line in sam_file:
	if line[0] == '@':
		continue

	fields = line.split("\t")   # makes a list of individual tab-separated fields

	flag =  int(fields[1])  # bitwise flag
	pos =   int(fields[3])  # 1-based leftmost mapping positionition 
	pnext = int(fields[7])  # position of the mate/next read 
	tlen =  int(fields[8])  # observed template length 
	seq =   fields[9]       # segment sequence
 
	# Both reads align correctly and len is greater than zero
	if ((flag & 4) == 0 and tlen > 0):   # flag unset indicates that segment maps 
		# Increment start position by one
		genome_change[pos] += 1
		# Decrement fragment end position by one
		genome_change[pnext + len(seq)] -= 1

# Close the sam file
sam_file.close()
	
# Print genomic profile as a wiggle file
print("fixedStep chrom=genome start=1 step=1 span=1")

current_coverage = 0

# Cicle over all positions of the genome
for position in range(genome_length):
	current_coverage += genome_change[position]
	print(current_coverage)
