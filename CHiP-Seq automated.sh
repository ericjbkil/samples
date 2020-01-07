#!/bin/sh

FASTQ_DIR=/data00/datafiles/concatenate-hif-chip
TAG_DIR=/data00/datafiles/bowtie-chip
TRIM_DIR=/data00/ekil/ChIP

trim_galore --fastqc -j 4 -o $TRIM_DIR $FASTQ_DIR/cIP1.fastq $FASTQ_DIR/input.fastq

makeTagDirectory cIP1_tag $TAG_DIR/cIP1bowtie2-mm10.sam
makeTagDirectory input_tag $TAG_DIR/Input-mm10.sam

makeUCSCfile cIP1_tag -o auto -fsize 2e7

findPeaks cIP1_tag -style factor -F 2 -P .05 -o auto -i Input_tag

findMotifsGenome.pl peaks.txt $TRIM_DIR/mm10 $TRIM_DIR/motifs/

annotatePeaks.pl peaks.txt $TRIM_DIR/mm10/ > annotations.txt

