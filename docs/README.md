# Report

The main goal of this document is to show how to implement a simple
analysis pipeline to detect structural variations on a reference
genome by aligning unknown bacterium reads.

The reference genome is from a bacterium called _Lactobacillus casei_
which is _3079196 bp_ long.

## Prerequisities

We need some tools in order to let the pipeline works properly:

- _samtools_, the tool for alignments in the SAM format;
- _bwa_, the Burrows-Wheeler Aligner is a software package for mapping
  low-divergent sequences against a large reference genome;
- _igw_, the Integrative Genomics Viewer is a high-performance,
  easy-to-use, interactive tool for the visual exploration of genomic
  data;

and we also need some data:

- *Lactobacillus_casei_genome.fasta*, the fasta file for the reference genome;
* *lact_sp.read1.fastq* and *lact_sp.read2.fastq*, the fastq file with reads for donor genome;

## Project

### Makefile

TODO

### Utils

In order to simplify some tasks we implemented the `sam_utils` module. In this module we grouped some function for multiple purposes.

Here follows the list of functions with their descriptions.

* `read_mates`, a helper function for reading pair mates from `.sam` files. 

* `to_wig`, a helper function for printing a sequence to `.wig` format.

* `is_plausible_tlen`, a predicate that checks if the `tlen` field from a mate pair is a plausible template length, i.e., it is greater than *0* and less than a threshold usually set to *20000*.

* `is_first_read_exlusively_mapped`, a predicate that checks using the `flag` field from a mate pair if only the first read is mapping, while the second does not.

* `is_second_read_exlusively_mapped`, the dual predicate of `is_first_read_exlusively_mapped`.

* `is_first_and_second_read_mapped`, a predicate that checks using the `flag` field from mate pair if both reads map.

* `is_mapping_fr`, a predicate that checks with `flag` field from mate pair whether reads map in FR schema.

* `filter_out_invalid_mates`, an helper function that given a list of mate pairs, filters out invalid mate paris, i.e. mate pairs with implausible template length.

### 1. Sequence Coverage

> Implemented in `sequence_coverage.py`

Sequence coverage is a number that, for every genomic position, counts number of reads mapping there. The count is done only for the read.

More precisely, sequence coverage is the average number of reads that align to, or "cover", known reference bases.

The sequencing coverage level often determines whether variant discovery can be made with a certain degree of confidence at particular base positions, in fact, higher levels of coverage means each base is covered by a greater number of aligned sequence read.

### 2. Physical Coverage

> Implemented in `physical_coverage.py`

Physical coverage is a number that, for every genomic position, counts number of fragments mapping there. The main difference between sequence coverage and physical coverage is that in the first the count is done on reads while the second is done on fragments.

Physical coverage in implemented in module `physical_coverage` through three main functions:

* `get_physical_change`, compute the physical coverage change on the genome setting *+= 1* where the fragment starts and *-= 1* where the fragment ends.

* `get_sum_fragment_change`, compute the real physical coverage values by summing changes for every genomic position.

* `get_physical_coverage_percentage`, compute the physical coverage percentage by scaling the physical coverage values by fragments length.

### 3. Single Mates

> Implemented in `single_mates.py`

Single mates are mate pairs where only one of the two reads maps on the genome.



## Reference

1. https://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
2. https://towardsdatascience.com/pdf-is-not-a-probability-5a4b8a5d9531
