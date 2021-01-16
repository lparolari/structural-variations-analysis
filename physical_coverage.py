#!/usr/bin/env python3

import sys
import logging

from sam_utils import FLAG_SEGMENT_UNMAPPED, FLAG_NEXT_SEGMENT_UNMAPPED, FLAG_UNSET
from sam_utils import MATE_LENGTH
from sam_utils import read_mates, to_wig
from sam_utils import filter_out_invalid_mates
from sam_utils import is_first_and_second_read_mapped


def print_usage():
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
    print(f"  $ {sys.argv[0]} lact.sam 3079196 > lact_physicalcov.wig")


def get_physical_change(mates, genome_length):
    """
    Return an array where every genomic position represents the 
    number of mates starting at that genomic position.
    """
    genome_physical_change = [0] * genome_length

    for mate in mates:
        pos = mate["pos"]
        pnext = mate["pnext"]

        if is_first_and_second_read_mapped(mate):
            
            first_mate_starting_position = pos
            second_mate_finish_position = pnext + MATE_LENGTH

            genome_physical_change[first_mate_starting_position] += 1
            genome_physical_change[second_mate_finish_position] -= 1

    return genome_physical_change


def get_sum_fragment_change(mates, genome_length):
    """
    Return an array where every genomic position represents the sum
    of the length of fragments with the first mate starting at that
    genomic position.
    """
    genome_sum_fragment_change = [0] * genome_length

    for mate in mates:
        pos = mate["pos"]
        pnext = mate["pnext"]

        if is_first_and_second_read_mapped(mate):
            
            first_mate_starting_position = pos
            second_mate_finish_position = pnext + MATE_LENGTH

            fragment_length = second_mate_finish_position - first_mate_starting_position

            genome_sum_fragment_change[first_mate_starting_position] += fragment_length
            genome_sum_fragment_change[second_mate_finish_position] -= fragment_length

    return genome_sum_fragment_change


def get_physical_coverage_percentage(mates, genome_length):
    """
    Return an array with physical coverage percentage for every
    genomic position.

    Parameters
    ----------
    * mates: A list of plausible mates.
    * genome_length: Length of genome.
    """
    genome_physical_coverage = []

    genome_physical_change = get_physical_change(mates, genome_length)
    genome_sum_fragment_change = get_sum_fragment_change(mates, genome_length)

    current_physical_coverage = 0
    current_sum_fragment_length = 0

    for i in range(genome_length):
        current_physical_coverage += genome_physical_change[i]
        current_sum_fragment_length += genome_sum_fragment_change[i]

        if(current_physical_coverage > 0):
            coverage = current_sum_fragment_length / current_physical_coverage
        else:
            coverage = 0
        
        genome_physical_coverage.append(coverage)
    
    return genome_physical_coverage


if __name__ == "__main__":
    
    # Set verbose logging
    logging.basicConfig(level=logging.DEBUG)

    # Check args
    if len(sys.argv) != 1 + 2:
        print_usage()
        exit(1)

    # Get args
    input_file = sys.argv[1]
    genome_length = int(sys.argv[2])

    logging.debug(f"input_file = {input_file}")
    logging.debug(f"genome_length = {genome_length}")

    # Get mates
    mates = read_mates(input_file, keep_fields=["pos", "pnext", "flag", "tlen"])
    mates = filter_out_invalid_mates(mates)

    # Compute physical coverage
    physical_coverage_percentage = get_physical_coverage_percentage(mates, genome_length)
    assert len(physical_coverage_percentage) == genome_length, f"Track length must be equal to genome lenth, but {len(physical_coverage_percentage)} != {genome_length}"

    # Print to wig file
    to_wig(physical_coverage_percentage)
