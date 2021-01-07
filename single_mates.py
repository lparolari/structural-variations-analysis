import sys
import logging
import numpy as np
from typing import List

from sam_utils import FLAG_UNSET, FLAG_SEGMENT_UNMAPPED, FLAG_NEXT_SEGMENT_UNMAPPED
from sam_utils import read_mates, to_wig


MATE_LENGTH = 100


def print_usage(argv):
    print(f"Syntax: {argv[0]} FILE.sam GENOME_LENGTH")
    print()
    print("Description:")
    print("  Returns a track with percentage of single mates.")
    print("  Single mates are those where only one of the two reads is mapping.")
    print()
    print("Example:")
    print(f"  $ {argv[0]} lact.sam 3079196")
    print(f"  $ {argv[0]} lact.sam 3079196 > lact_singlemates.wig")
    exit(1)


def is_plausible_tlen(mate):
    """
    Returns true whether the tlen of the mate is geq 0 and leq 20000, false otherwise.
    """
    return (mate["tlen"] >= 0) and (mate["tlen"] <= 20000)

def is_first_read_exlusively_mapped(mate):
    """
    Returns true whether the first read maps exclusively, false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_SEGMENT_UNMAPPED


def is_second_read_exlusively_mapped(mate): 
    """
    Returns true whether the second read maps exclusively, false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_NEXT_SEGMENT_UNMAPPED
    

def is_first_and_second_read_mapped(mate): 
    """
    Returns true whether first and second read are mapped, false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_UNSET


def get_tlen_distribution_params(mates):
    """
    Returns the average and standard deviation of mates with template length geq 0
    and leq 2000. Two mates must map in order to use their template length.
    """

    add_second_read_length = lambda mate: mate["tlen"] + MATE_LENGTH 

    tlens = list(filter(is_plausible_tlen, mates))
    tlens = list(filter(is_first_and_second_read_mapped, tlens))
    tlens = list(map(add_second_read_length, tlens))

    logging.debug(f"number_of_mapping_reads = {len(tlens)}")

    # Get mean and std
    tlen_avg = np.mean(np.array(tlens))
    tlen_std = np.std(np.array(tlens))

    return tlen_avg, tlen_std


def get_single_mates_count(mates):
    """
    Returns the number of single mates, i.e., mates where only one of the 
    two maps.
    """

    # Get the number of mates where, respectively, the first read maps exclusively and
    # the second read maps exclusively
    singles_first_no = len(list(filter(is_first_read_exlusively_mapped, mates)))
    singles_second_no = len(list(filter(is_second_read_exlusively_mapped, mates)))

    # Sum the two counts
    single_mates_no = singles_first_no + singles_second_no

    return single_mates_no


def get_single_mates(mates, genome_length):
    """
    Returns an array with 
    *  1, at the start of one unmapped mate
    * -1, at the end of one unmapped mate
    *  0, elsewhere 
    """

    # Get average and standard deviation for distance mate pairs
    tlen_avg, tlen_std = get_tlen_distribution_params(mates)

    logging.info(f"tlen_avg = {tlen_avg}")
    logging.info(f"tlen_avg = {tlen_std}")
    
    # Initialize array for storing all single mate pairs percentage
    single_mates = [0] * genome_length

    for mate in mates:
        pos = mate["pos"]
        
        if is_first_read_exlusively_mapped(mate):
            # First mate maps, while the second does not
            mapped_mate = pos
            unmapped_mate = mapped_mate + int(tlen_avg)

            single_mates[unmapped_mate] += 1
            single_mates[unmapped_mate + MATE_LENGTH] -= 1
        
        if is_second_read_exlusively_mapped(mate):
            # First mate does not map, while the second does
            mapped_mate = pos
            unmapped_mate = mapped_mate - int(tlen_avg)

            single_mates[unmapped_mate] += 1
            single_mates[unmapped_mate + MATE_LENGTH] -= 1
    
    return single_mates


def get_single_mates_percentage(mates, genome_length):
    """
    Returns single mates percentage by computing the sum of single mates on each
    genomic position, dividing by single mates count and multiplying by 100.
    """

    single_mates_percentage = []

    single_mates = get_single_mates(mates, genome_length)
    single_mates_count = get_single_mates_count(mates)

    logging.info(f"single_mates_count = {single_mates_count}")

    current_single_mates = 0

    for i in range(genome_length):
        current_single_mates += single_mates[i]
        single_mates_percentage.append(current_single_mates * 100 / single_mates_count)

    return single_mates_percentage


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

    # Compute single mates percentage
    single_mates_percentage = get_single_mates_percentage(mates, genome_length)

    # Print single_mates_percentage in wig format
    to_wig(single_mates_percentage)
