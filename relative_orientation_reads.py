import sys
import logging
import numpy as np
from typing import List

from sam_utils import FLAG_UNSET, FLAG_SEGMENT_UNMAPPED, FLAG_NEXT_SEGMENT_UNMAPPED
from sam_utils import MATE_LENGTH
from sam_utils import read_mates, to_wig
from sam_utils import filter_out_invalid_mates
from sam_utils import is_mapping_fr


def print_usage(argv):
    print(f"Syntax: {argv[0]} FILE.sam GENOME_LENGTH")
    print()
    print("Description:")
    print("  Returns a track with relative orientation for every genomic position.")
    print()
    print("Example:")
    print(f"  $ {argv[0]} lact.sam 3079196")
    print(f"  $ {argv[0]} lact.sam 3079196 > lact_orientation.wig")
    exit(1)


def get_relative_orientations_change(mates, genome_length):
    relative_orientations = [0] * genome_length

    # With standard sequencing technologies reads should align 
    # in FR mode, if they align FF, RR or RF thay are probably
    # wrong.

    mates = list(filter(is_mapping_fr, mates))

    for mate in mates:
        pos = mate["pos"]
        pnext = mate["pos"]

        relative_orientations[pos] += 1
        relative_orientations[pnext + MATE_LENGTH] -= 1

    return relative_orientations


def get_relative_orientations(mates, genome_length):
    relative_orientations = [0] * genome_length
    
    changes = get_relative_orientations_change(mates, genome_length)
    current_change = 0
    
    for i in range(genome_length):
        current_change += changes[i]
        relative_orientations[i] = current_change
    
    return relative_orientations


if __name__ == "__main__":
    # Set verbose logging
    logging.basicConfig(level=logging.DEBUG)

    # Check args
    if len(sys.argv) != 1 + 2:
        print_usage(sys.argv)

    # Get args
    input_file = sys.argv[1]
    genome_length = int(sys.argv[2])

    logging.debug(f"input_file = {input_file}")
    logging.debug(f"genome_length = {genome_length}")

    # Read mates
    mates = read_mates(input_file, keep_fields=["pos", "pnext", "tlen", "flag"])
    mates = filter_out_invalid_mates(mates)

    # Compute single mates percentage
    relative_orientations = get_relative_orientations(mates, genome_length)
    assert len(relative_orientations) == genome_length, f"Track length must be equal to genome lenth, but {len(relative_orientations)} != {genome_length}"

    # Print single_mates_percentage in wig format
    to_wig(relative_orientations)
