import sys
import logging
import numpy as np
from typing import List

FLAG_SEGMENT_UNMAPPED = 4
FLAG_NEXT_SEGMENT_UNMAPPED = 8


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


def read_mates(input_file: str, raw_fields: bool = False, keep_comments: bool = False):
    """
    Return a list of dict where every entry is representing a field from the
    sam file.

    :param raw_fields: Allow to do not unpack fields and get the whole text line
    :param keep_comments: Allow to keep comments in sam file, i.e, line with `@` at the beginning
    """
    with open(input_file) as f:
        lines = f.readlines()

        if not keep_comments:
            lines = list(filter(lambda line: line[0] != "@", lines))
        
        if not raw_fields:
            def split_fields(line):
                values = line.split("\t")

                fields = {
                    "qname": values[0],
                    "flag": int(values[1]),
                    "rname": values[2],
                    "pos": int(values[3]),
                    "mapq": int(values[4]),
                    "cigar": values[5],
                    "rnext": values[6],
                    "pnext": int(values[7]),
                    "tlen": int(values[8]),
                    "seq": values[9],
                    "qual": values[10]
                }

                return fields

            lines = list(map(split_fields, lines))
        
        return lines


def is_first_read_unmapped(mate):
    """
    Returns true whether the first read maps exclusively, false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_SEGMENT_UNMAPPED


def is_second_read_unmapped(mate): 
    """
    Returns true whether the second read maps exclusively, false otherwise.
    """
    return (mate["flag"] & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_NEXT_SEGMENT_UNMAPPED
    

def is_first_or_second_unmapped(mate): 
    """
    Returns true whether first read maps exclusively or second read maps esclusively,
    false otherwise.
    """
    return is_first_read_unmapped(mate) or is_second_read_unmapped(mate)


def get_tlen_distribution_params(mates):
    """
    Returns the average and standard deviation of template length only for 
    mates with tlen
    * >= 0
    * <= 20000
    * only one of the two mates maps
    
    Please note that avg and std are calculated on unmapped tlen.
    """

    tlens = list(filter(lambda mate: (mate["tlen"] >= 0) and (mate["tlen"] <= 20000), mates))
    tlens = list(filter(lambda mate: is_first_or_second_unmapped, tlens))
    tlens = list(map(lambda mate: mate["tlen"] + len(mate["seq"]), tlens))

    logging.debug(f"len(tlens) = {len(tlens)}")

    tlen_avg = np.mean(np.array(tlens))
    tlen_std = np.std(np.array(tlens))

    return tlen_avg, tlen_std


def get_single_mates_count(mates):
    """
    Returns the number of single mates, i.e., mates where only one of the 
    two maps.
    """

    return len(list(filter(is_first_or_second_unmapped, mates)))


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
        flag = mate["flag"]
        pos = mate["pos"]
        seq = mate["seq"]
        
        if ((flag & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_NEXT_SEGMENT_UNMAPPED):
            # First mate maps, while the second does not
            mapped_mate = pos
            unmapped_mate = mapped_mate + int(tlen_avg)

            single_mates[unmapped_mate] += 1
            single_mates[unmapped_mate + len(seq)] -= 1
        
        if ((flag & (FLAG_SEGMENT_UNMAPPED | FLAG_NEXT_SEGMENT_UNMAPPED)) == FLAG_NEXT_SEGMENT_UNMAPPED):
            # First mate does not map, while the second does
            mapped_mate = pos
            unmapped_mate = mapped_mate - int(tlen_avg)

            single_mates[unmapped_mate] += 1
            single_mates[unmapped_mate + len(seq)] -= 1
    
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


def to_wig(ls, step_type: str = "fixedStep", chrom: str = "genome", start: int = 1, step: int = 1, span: int = 1):
    """
    Prints to the standard output a list in wig format.
    """
    print(f"{step_type} chrom={chrom} start={start} step={step} span={span}")

    for item in ls:
        print(item)


def main():

    logging.basicConfig(level=logging.DEBUG)

    if len(sys.argv) != 1 + 2:
        print_usage(sys.argv)

    # Get args
    input_file = sys.argv[1]
    genome_length = int(sys.argv[2])

    logging.debug(f"input_file = {input_file}")
    logging.debug(f"genome_length = {genome_length}")

    # Read mates
    mates = read_mates(input_file)

    # Compute single mates percentage
    single_mates_percentage = get_single_mates_percentage(mates, genome_length)

    # Print single_mates_percentage in wig format
    to_wig(single_mates_percentage)


if __name__ == "__main__":
    main()
