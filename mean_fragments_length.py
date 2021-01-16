import sys
import logging
import numpy as np
from typing import List

from sam_utils import FLAG_UNSET, FLAG_SEGMENT_UNMAPPED, FLAG_NEXT_SEGMENT_UNMAPPED
from sam_utils import MATE_LENGTH
from sam_utils import read_mates, to_wig
from sam_utils import filter_out_invalid_mates
from sam_utils import is_first_and_second_read_mapped


def print_usage(argv):
    print(f"Syntax: {argv[0]} FILE.sam GENOME_LENGTH")
    print()
    print("Description:")
    print("  Returns a track with mean fragments length.")
    print()
    print("Example:")
    print(f"  $ {argv[0]} lact.sam 3079196")
    print(f"  $ {argv[0]} lact.sam 3079196 > lact_avgfraglen.wig")
    exit(1)


def get_framents_length(mates, genome_length):
    """
    Returns an array where for every genomic position we have
    computed the fragment length.
    """

    fragments_length = [0] * genome_length

    for mate in mates:
        pos = mate["pos"]
        pnext = mate["pnext"]
        tlen = mate["tlen"]
        
        if is_first_and_second_read_mapped(mate):

            first_mate_starting_position = pos
            second_mate_finish_position = pnext + MATE_LENGTH

            # For each genomic position we sum the its current tlen
            # from start mate pair to end mate pair
            fragments_length[first_mate_starting_position] += tlen
            fragments_length[second_mate_finish_position] -= tlen
   
    return fragments_length


def get_physical_coverage(mates, genome_length):
    """
    Return the physical coverage for mates without percentage.
    """    
    physical_coverage = [0] * genome_length

    for mate in mates:
        pos = mate["pos"]
        pnext = mate["pnext"]
        
        if is_first_and_second_read_mapped(mate):

            first_mate_starting_position = pos
            second_mate_finish_position = pnext + MATE_LENGTH

            physical_coverage[first_mate_starting_position] += 1
            physical_coverage[second_mate_finish_position] -= 1

    return physical_coverage


def get_avg_fragments_length(mates, genome_length):
    """
    Return an array where every genomic position represents the average
    length of fragments at that position.
    """
    avg_fragments_length = [0] * genome_length

    physical_coverage = get_physical_coverage(mates, genome_length)
    fragments_length = get_framents_length(mates, genome_length)

    current_phisical_coverage = 0
    current_fragments_length = 0
    
    for i in range(genome_length):
        current_phisical_coverage += physical_coverage[i]
        current_fragments_length += fragments_length[i]

        if (current_phisical_coverage == 0):
            avg_fragments_length[i] = 0
            continue
    
        result = current_fragments_length / current_phisical_coverage
        avg_fragments_length[i] = int(result)

    return avg_fragments_length


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

    # Compute avg frag length
    avg_fragments_length = get_avg_fragments_length(mates, genome_length)

    # Print single_mates_percentage in wig format
    to_wig(avg_fragments_length)
