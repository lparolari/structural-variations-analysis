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

${SEQ_COVERAGE}.wig:
	source ./.venv/bin/activate; python3 sequence_coverage.py ${ALIGNED_SAM} 3079196 > ${SEQ_COVERAGE}.wig

${PHYSICAL_COVERAGE}.wig:
	source ./.venv/bin/activate; python3 physical_coverage.py ${ALIGNED_SAM} 3079196 > ${PHYSICAL_COVERAGE}.wig

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

.PHONY: clear clean