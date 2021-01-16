SHELL=/bin/bash
GENOME=Lactobacillus_casei_genome
READ_1=lact_sp.read1
READ_2=lact_sp.read2
ALIGNED=lact
ALIGNED_SAM=${ALIGNED}.sam
ALIGNED_BAM=${ALIGNED}.bam
ALIGNED_SORTED=${ALIGNED}_sorted.bam
ALIGNED_SORTED_INDEXED=${ALIGNED}_sorted.bam.bai
SEQ_COVERAGE=lact_seqcov
PHYSICAL_COVERAGE=lact_physicalcov
SINGLE_MATES=lact_singlemates
AVG_FRAGMENTS_LENGTH=lact_avgfraglen
RELATIVE_ORIENTATIONS=lact_orientation
PROB_INS=lact_probins
PROB_DEL=lact_probdel

### **********************************************************
### MAIN TARGETS

all: aligned
.PHONY: all

aligned: ${ALIGNED_SORTED} ${ALIGNED_SORTED_INDEXED}
.PHONY: aligned

### **********************************************************
### RECIPES

# Genome indexing
${GENOME}.fasta.amb: index_genome
${GENOME}.fasta.ann: index_genome
${GENOME}.fasta.bwt: index_genome
${GENOME}.fasta.fai: index_genome
${GENOME}.fasta.pac: index_genome
${GENOME}.fasta.sa: index_genome

index_genome:
	./bwa index ${GENOME}.fasta
.PHONY: index_genome

# Genome alignment with reads
${ALIGNED_SAM}: ${GENOME}.fasta.amb ${GENOME}.fasta.ann ${GENOME}.fasta.bwt ${GENOME}.fasta.fai ${GENOME}.fasta.pac ${GENOME}.fasta.sa
	./bwa mem ${GENOME}.fasta ${READ_1}.fastq ${READ_2}.fastq > ${ALIGNED_SAM}

# From .sam to .bam
${ALIGNED_BAM}: ${ALIGNED_SAM}
	./samtools view -bS ${ALIGNED_SAM} > ${ALIGNED_BAM}

# Sort alignmnets by genomic positions
${ALIGNED_SORTED}: ${ALIGNED_BAM}
	./samtools sort ${ALIGNED_BAM} > ${ALIGNED_SORTED}

# Index sorted file
${ALIGNED_SORTED_INDEXED}: ${ALIGNED_SORTED}
	./samtools index ${ALIGNED_SORTED} > ${ALIGNED_SORTED_INDEXED}

# Create sequence coverage track
${SEQ_COVERAGE}.wig:
	source ./.venv/bin/activate; python3 sequence_coverage.py ${ALIGNED_SAM} 3079196 > ${SEQ_COVERAGE}.wig

# Create physical coverage track
${PHYSICAL_COVERAGE}.wig:
	source ./.venv/bin/activate; python3 physical_coverage.py ${ALIGNED_SAM} 3079196 > ${PHYSICAL_COVERAGE}.wig

# Create single mates track
${SINGLE_MATES}.wig:
	source ./.venv/bin/activate; python3 single_mates.py ${ALIGNED_SAM} 3079196 > ${SINGLE_MATES}.wig

# Create avg fragments length tracks
${AVG_FRAGMENTS_LENGTH}.wig:
	source ./.venv/bin/activate; python3 mean_fragments_length.py ${ALIGNED_SAM} 3079196 > ${AVG_FRAGMENTS_LENGTH}.wig

# Create relative orientation track
${RELATIVE_ORIENTATIONS}.wig:
	source ./.venv/bin/activate; python3 relative_orientation_reads.py ${ALIGNED_SAM} 3079196 > ${RELATIVE_ORIENTATIONS}.wig

# Create prob ins track
${PROB_INS}.wig:
	source ./.venv/bin/activate; python3 fragments_length_distribution.py ${ALIGNED_SAM} 3079196 --track insertion > ${PROB_INS}.wig

# Create prob del track
${PROB_DEL}.wig:
	source ./.venv/bin/activate; python3 fragments_length_distribution.py ${ALIGNED_SAM} 3079196 --track deletion > ${PROB_DEL}.wig

### **********************************************************
### UTILS

clean: clear
clear:
	# Remove index utils
	rm -f ${GENOME}.fasta.amb ${GENOME}.fasta.ann ${GENOME}.fasta.bwt ${GENOME}.fasta.fai ${GENOME}.fasta.pac ${GENOME}.fasta.sa
	# Remove aligmnet files
	rm -f ${ALIGNED_SAM} 
	rm -f ${ALIGNED_BAM} 
	rm -f ${ALIGNED_SORTED} 
	rm -f ${ALIGNED_SORTED_INDEXED}
	# Remove tmp.*.bam files
	rm -f *.tmp.*.bam
	# Remove .wig files
	rm -f ${SEQ_COVERAGE}.wig
	rm -f ${PHYSICAL_COVERAGE}.wig
	rm -f ${SINGLE_MATES}.wig
	rm -f ${AVG_FRAGMENTS_LENGTH}.wig
	rm -f ${RELATIVE_ORIENTATIONS}.wig
	rm -f ${PROB_INS}.wig
	rm -f ${PROB_DEL}.wig

.PHONY: clear clean