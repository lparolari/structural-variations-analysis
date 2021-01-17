import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from scipy import stats
import sys
from typing import List

from sam_utils import FLAG_UNSET, FLAG_SEGMENT_UNMAPPED, FLAG_NEXT_SEGMENT_UNMAPPED
from sam_utils import MATE_LENGTH
from sam_utils import read_mates, to_wig
from sam_utils import filter_out_invalid_mates
from sam_utils import is_first_and_second_read_mapped
from sam_utils import is_first_read_exlusively_mapped, is_second_read_exlusively_mapped

from mean_fragments_length import get_avg_fragments_length


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


def get_distribution_space(mu, std):
    """
    Return a mock space for a basic pdf centered in `mu` with samples 
    between `± 3 * std`.
    """

    # We want to keep a finite part of the distribution and
    # we keep only values between `± 3 * std` from `mu`.
    # This values are in 99.7% pdf.

    return np.linspace(mu - 3 * std, mu + 3 * std, 100)


def get_pdf(x, mu, std):
    """
    Return the `y` for pdf with mean `mu` and standard deviation `std`
    for values on `x`.
    """
    return stats.norm.pdf(x, mu, std)


def get_probability_for_insertion(x, mu, std):
    """
    Return the probability of `x` to be an insertion based on the pdf
    obtained by `mu` and `std`. Probability of insertion is higher if
    average fragments length is less than `mu`.

    Please note that `x` can be an array.
    """
    # stats.norm.sf is the survival function
    
    return stats.norm.sf(x, mu, std)


def get_probability_for_deletion(x, mu, std):
    """
    Return the probability of `x` to be a deletion based on the pdf
    obtained by `mu` and `std`. Probability of deletion is higher if
    average fragments length is more than `mu`.

    Please note that `x` can be an array.
    """
    # stats.norm.cdf is the cumulative distrib function
    return stats.norm.cdf(x, mu, std)


def plot_pdf(x, y, mu, std):
    name = f"probability density function"

    plt.title(name)
    plt.plot(x, y, label=f"pdf")
    plt.legend()
    plt.xlabel("fragments avg length")
    plt.ylabel("density of probability")
    plt.savefig(f"{name}.png")
    plt.close()


def plot_insertion_probability(x, y, name="plot"):
    name = "Probability to be insertion"

    plt.title(name)
    plt.plot(x, y, label="probability of insertion")
    plt.legend()
    plt.xlabel("fragments avg length")
    plt.ylabel("probability")
    plt.savefig(f"{name}.png")
    plt.close()


def plot_deletion_probability(x, y):
    name = "Probability to be deletion"

    plt.title(name)
    plt.plot(x, y, label="probability of deletion")
    plt.legend()
    plt.xlabel("fragments avg length")
    plt.ylabel("probability")
    plt.savefig(f"{name}.png")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
        Computes the probability distribution for fragments length and plots 
        two charts respectively with probability of insertion and deletion or
        produces two tracks.
    """)

    # Main args
    parser.add_argument("file", type=pathlib.Path, help="A .sam file")
    parser.add_argument("genome_length", type=int, help="Genome length")
    
    # Options
    parser.add_argument("--verbose", type=bool, default=False, help="Verbose output")
    parser.add_argument("--plot", default=False, dest="plot", action="store_true",
        help="Plot and save chart with probability distribution for insertions and deletions")
    parser.add_argument("--track", choices=["insertion", "deletion"],
        help="""Compute the probability to be an insertion/deletion for every genomic position and prints it as a wig track""")

    # Get args
    args = parser.parse_args()
    
    input_file = args.file
    genome_length = args.genome_length
    verbose = args.verbose
    plot = args.plot
    track = args.track

    # Set loggin level
    logging.basicConfig(level=(logging.DEBUG if verbose else None))
    
    logging.debug(f"input_file = {input_file}")
    logging.debug(f"genome_length = {genome_length}")

    # Check if we have a task to do
    if not plot and not track:
        logging.warning("Nothing to do, maybe you should specify --plot or --track")
        exit(0)

    # Read mates
    mates = read_mates(input_file, keep_fields=["pos", "pnext", "tlen", "flag"])
    mates = filter_out_invalid_mates(mates)

    # Get tlen distribution params
    mu, std = get_tlen_distribution_params(mates)

    # Plot distribution
    if plot:
        x = get_distribution_space(mu, std)
        prob_insertion = get_probability_for_insertion(x, mu, std)
        prob_deletion = get_probability_for_deletion(x, mu, std)
        pdf = get_pdf(x, mu, std)

        plot_pdf(x, pdf, mu, std)
        plot_insertion_probability(x, prob_insertion)
        plot_deletion_probability(x, prob_deletion)

    if track == "insertion":
        logging.info("Computing insertion probabilities...")

        ls = get_avg_fragments_length(mates, genome_length)
        ls = get_probability_for_insertion(ls, mu, std)

        to_wig(ls)

    if track == "deletion":
        logging.info("Computing deletion probabilities...")
        
        ls = get_avg_fragments_length(mates, genome_length)
        ls = get_probability_for_deletion(ls, mu, std)

        to_wig(ls)
