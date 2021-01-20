import argparse
import logging
import numpy as np
import pathlib

from sam_utils import MATE_LENGTH, FLAG_SEGMENT_UNMAPPED, FLAG_UNSET
from sam_utils import read_mates, to_wig


def get_hard_soft_clippings_change(mates, genome_length):
	hard_soft_clippings_change = [0] * genome_length

	is_read_mapped = lambda mate: (mate["flag"] & FLAG_SEGMENT_UNMAPPED) == FLAG_UNSET
	is_hard_clipping = lambda mate: "H" in mate["cigar"]
	is_soft_clipping = lambda mate: "S" in mate["cigar"]
	is_hard_or_soft_clipping = lambda mate: is_hard_clipping(mate) or is_soft_clipping(mate)

	mates = list(filter(is_read_mapped, mates))  # keep mapped reads
	mates = list(filter(is_hard_or_soft_clipping, mates))  # keep reads with H or S in cigar field

	for mate in mates:
		pos = mate["pos"]

		hard_soft_clippings_change[pos] += 1
		hard_soft_clippings_change[pos + MATE_LENGTH] -= 1

	return hard_soft_clippings_change


def get_hard_soft_clippings(mates, genome_length):
	hard_soft_clippings = [0] * genome_length
	hard_soft_clippings_change = get_hard_soft_clippings_change(mates, genome_length)

	current_change = 0

	for i in range(genome_length):
		current_change += hard_soft_clippings_change[i]
		hard_soft_clippings[i] = current_change
	
	return hard_soft_clippings


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="""
		Computes a track where, for every genomic position, we count
		the number of reads with H or S in CIGAR field.
	""")

	# Main args
	parser.add_argument("file", type=pathlib.Path, help="A genome file in .sam format")
	parser.add_argument("genome_length", type=int, help="The genome length")
	
	# Options
    parser.add_argument("--verbose", default=False, dest="verbose", action="store_true", help="Verbose output")

	# Get args
	args = parser.parse_args()
	
	input_file = args.file
	genome_length = args.genome_length
	verbose = args.verbose

	# Set loggin level
	logging.basicConfig(level=(logging.DEBUG if verbose else None))
	
	logging.debug(f"input_file = {input_file}")
	logging.debug(f"genome_length = {genome_length}")

	# Read mates
	mates = read_mates(input_file, keep_fields=["pos", "flag", "cigar"])

	# Get multiple aligment track
	hard_soft_clippings = get_hard_soft_clippings(mates, genome_length)

	# Print track to wig file
	to_wig(hard_soft_clippings)
