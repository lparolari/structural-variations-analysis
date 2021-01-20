import argparse
import logging
import numpy as np
import pathlib

from sam_utils import read_mates, to_wig
from sam_utils import MATE_LENGTH


def get_multiple_alignments_change(mates, genome_length):
	multiple_alignments_change = [0] * genome_length

	for mate in mates:
		pos = mate["pos"]
		aligments = mate["ma"]  # other alignments where read maps, can also be empty!

		if (len(aligments) == 0):
			continue

		# Mark this read as mutiple alignment
		multiple_alignments_change[pos] += 1
		multiple_alignments_change[pos + MATE_LENGTH] -= 1

		# Mark each read identified by current as multiple
		for aligment in aligments:
			# We need to "abs" bacause of + or - preceding aligment position
			aligment_pos = abs(int(aligment[1]))

			multiple_alignments_change[aligment_pos] += 1
			multiple_alignments_change[aligment_pos + MATE_LENGTH] -= 1

	return multiple_alignments_change


def get_multiple_alignments(mates, genome_length):
	multiple_alignments = [0] * genome_length
	multiple_alignments_change = get_multiple_alignments_change(mates, genome_length)

	current_change = 0

	for i in range(genome_length):
		current_change += multiple_alignments_change[i]
		multiple_alignments[i] = current_change
	
	return multiple_alignments


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="""
		Computes a track where, for every genomic position, we count
		the number reads mapping there.
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
	mates = read_mates(input_file, keep_fields=["pos", "ma"])

	# Get multiple aligment track
	multiple_alignments = get_multiple_alignments(mates, genome_length)

	# Print track to wig file
	to_wig(multiple_alignments)
