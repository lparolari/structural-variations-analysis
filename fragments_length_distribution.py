import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import pathlib
import scipy
import sys
from typing import List

from sam_utils import FLAG_UNSET, FLAG_SEGMENT_UNMAPPED, FLAG_NEXT_SEGMENT_UNMAPPED
from sam_utils import MATE_LENGTH
from sam_utils import read_mates, to_wig
from sam_utils import filter_out_invalid_mates
from sam_utils import is_first_and_second_read_mapped
from sam_utils import is_first_read_exlusively_mapped, is_second_read_exlusively_mapped


def get_tlen_distribution_params(mates):
    """
    Returns the average and standard deviation for mates tlen.
    """

    add_second_read_length = lambda mate: mate["tlen"] + MATE_LENGTH 

    tlens = mates
    tlens = list(filter(is_first_and_second_read_mapped, tlens))
    tlens = list(map(add_second_read_length, tlens))

    logging.debug(f"number_of_mapping_reads = {len(tlens)}")

    # Get mean and std
    tlen_avg = np.mean(np.array(tlens))
    tlen_std = np.std(np.array(tlens))

    logging.debug(f"mu = {tlen_avg}")
    logging.debug(f"std = {tlen_std}")

    return tlen_avg, tlen_std


def plot_distribution(mu, std):
    """
    Plo
    """
    from scipy import stats

    # We want to plot a finite part of the distribution and
    # we keep only points between `± 3 * std` from `mu`.
    # This values are the 99.7% of all values.
    x = np.linspace(mu - 3 * std, mu + 3 * std, 100)

    # We generate the pdf
    y = stats.norm.pdf(x, mu, std)

    assert np.abs(1 - np.trapz(y, x=x)) < 0.01, "Integral of distrib not equal to 1"

    plt.plot(x, y)
    plt.savefig("fragments-length-distribution.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Saves to fragment-length-distribution.png a plot with fragments length distribution.")

    parser.add_argument("file", type=pathlib.Path, help="A .sam file")
    parser.add_argument("genome_length", type=int, help="Genome length")
    parser.add_argument("--verbose", type=bool, default=False, help="Verbose output")

    # Get args
    args = parser.parse_args()
    
    input_file = args.file
    genome_length = args.genome_length
    verbose = args.verbose
    
    logging.basicConfig(level=(logging.DEBUG if verbose else None))

    logging.debug(f"input_file = {input_file}")
    logging.debug(f"genome_length = {genome_length}")

    # Read mates
    mates = read_mates(input_file, keep_fields=["pos", "pnext", "tlen", "flag"])
    mates = filter_out_invalid_mates(mates)

    # Get tlen distribution params
    mu, std = get_tlen_distribution_params(mates)

    # Plot distribution
    plot_distribution(mu, std)