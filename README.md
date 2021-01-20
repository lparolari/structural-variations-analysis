# Structural Variation Analysis

> Finding structural variations in a genome after comparison with a
> reference genome.

_Note_. This is an academic project for the course of Bioinformatics
(a.y. 2020-2021) at University of Padua.

## Goal

The main goal of this project is to analyze a genome through mate
pairs of an unkndown bacterium and find structural variation wrt
_Lactobacillus casei_ bacterium genome.

On the Bioinformatics 2020-21
[course page](https://elearning.unipd.it/math/course/view.php?id=647)
you can find the assignment and the genome files, but you need UNIPD
account in order to view it. In case you can't access send an email to
`luca.parolari@studenti.unipd.it`.

## [Documentation](docs/README.md)

You can find a report for the exam containing the documentation for
the project in [docs/README.md](docs/README.md). See [Docs](#docs)
section to build documentation in `pdf` format.

## Prerequisities

### Tools

You need some tools in order to work with this pipeline:

- [samtools](https://www.htslib.org/), the tool for alignments in the
  SAM format;

- [bwa](https://github.com/lh3/bwa), the Burrows-Wheeler Aligner is a
  software package for mapping low-divergent sequences against a large
  reference genome;

- [igv](http://software.broadinstitute.org/software/igv/), the
  Integrative Genomics Viewer is a high-performance, easy-to-use,
  interactive tool for the visual exploration of genomic data;

This repository includes a compiled version of _samtools_ and _bwa_
for Linux/Ubuntu distribution, so you may need only to download _igv_
from the
[download page](http://software.broadinstitute.org/software/igv/download)
and unpack it in the `tools` folder to start playing with this
project.

However is more secure to dowload and install (or make from sources)
the tools for your operating system.

### Data

You need also some data to process, in particular

- _Lactobacillus_casei_genome.fasta_, the fasta file for the reference
  genome;

* _lact_sp.read1.fastq_ and _lact_sp.read2.fastq_, the fastq file with
  reads for donor genome;

## Install

You just need to clone the project to your local machine.

## Usage

### Build

You can align the genome on reference genome and analyze it by running

```
make all
```

This command generates the aligned file with _bwa_ (with some index
files) and many analysis tracks in _wig_ format that can be easily
loaded on _igv_.

Take a look at all available recipes in [`Makefile`](Makefile).

Other useful "bulk" recipes are

- `make align`, to generate only the aligned file;
- `make analysis`, to generate all the analysis track;

You may also want to take a look at

- `make lact_seqcov.wig`, to build sequence coverage analysis;
- `make lact_physicalcov.wig`, to build physical coverage analysis;
- `make lact_singlemates.wig`, to build single mates count analysis;
- `make lact_avgfraglen.wig`, to build avg fragments length analysis;
- `make lact_orientation.wig`, to build relative reads orientation
  analysis;
- `make lact_probins.wig`, to build probability of insertion analysis;
- `make lact_probdel.wig`, to build probability of deletion analysis;
- `make lact_multiplealignments.wig`, to build multiple alignments
  analysis;
- `make lact_hsclipping.wig`, to build hard/soft clippings analysis;

### View

After building aligned genome file and analysis track you can load
with `igv` the file `igv_session.xml` and browse the genome to find
anomalies.

### Docs

You can also build documentation as a pdf with

```bash
cd docs
make
```

## Author

Luca Parolari

- Email: [luca.parolari23@gmail.com](mailto:luca.parolari23@gmail.com)
- GitHub: [@lparolari](https://github.com/lparolari)
- Telegram: [@lparolari](https://t.me/lparolari)

## License

This project is MIT licensed.
