#!/bin/sh

#Specify your directory paths below

FASTQ_DIR=/data00/datafiles/concatenate-nkx
TRIM_DIR=/data00/ekil/hw1pt2

#Step through the analysis steps below add comments so you 
#know what you are doing and can look back and make sense of everything

#run fastqc analysis on your seq files (works for fastq format)
#mkdir fastqcruns
#mv fastqcruns
#fastqc -t 2 -o $TRIM_DIR/ $FASTQ_DIR/*.fastq
#cd ..

# trim the rest of the paired end samples below

#mkdir 10_KO
#mkdir 5_CTRL
#mkdir 5_KO
#mkdir 2_CTRL
#mkdir 8_KO
#mkdir 9_CTRL

trim_galore --fastqc --paired -j 4 -o $TRIM_DIR/10_KO/ $FASTQ_DIR/10*_R1.fastq $FASTQ_DIR/10*_R2.fastq

trim_galore --fastqc --paired -j 4 -o $TRIM_DIR/5_CTRL/ $FASTQ_DIR/5_CTRL_R1.fastq $FASTQ_DIR/5_CTRL_R2.fastq

trim_galore --fastqc --paired -j 4 -o $TRIM_DIR/5_KO $FASTQ_DIR/5_KO_R1.fastq $FASTQ_DIR/5_KO_R2.fastq

trim_galore --fastqc --paired -j 4 -o $TRIM_DIR/2_CTRL/ $FASTQ_DIR/2*_R1.fastq $FASTQ_DIR/2*_R2.fastq

trim_galore --fastqc --paired -j 4 -o $TRIM_DIR/9_CTRL/ $FASTQ_DIR/9*_R1.fastq $FASTQ_DIR/9*_R2.fastq

trim_galore --fastqc --paired -j 4 -o $TRIM_DIR/8_KO/ $FASTQ_DIR/8*_R1.fastq $FASTQ_DIR/8*_R2.fastq

# Create an index with Kallisto. You have already done this so 
# you can comment this step out but keep it in here for future 
# reference

#kallisto index -i Mus_musculus.idx Mus_musculus.GRCm38.cdna.all.fa.gz

# Map reads to kallisto index for the rest of your samples 
kallisto quant -i Mus_musculus.idx -t 2 -o $TRIM_DIR/10_KO/ $TRIM_DIR/10_KO/10*_R1_val_1.fq $TRIM_DIR/10_KO/10*_R2_val_2.fq

kallisto quant -i Mus_musculus.idx -t 2 -o $TRIM_DIR/5_CTRL/ $TRIM_DIR/5_CTRL/5_CTRL_R1_val_1.fq $TRIM_DIR/5_CTRL/5_CTRL_R2_val_2.fq

kallisto quant -i Mus_musculus.idx -t 2 -o $TRIM_DIR/5_KO/ $TRIM_DIR/5_KO/5_KO_R1_val_1.fq $TRIM_DIR/5_KO/5_KO_R2_val_2.fq

kallisto quant -i Mus_musculus.idx -t 2 -o $TRIM_DIR/2_CTRL/ $TRIM_DIR/2_CTRL/2*_R1_val_1.fq $TRIM_DIR/2_CTRL/2*_R2_val_2.fq

kallisto quant -i Mus_musculus.idx -t 2 -o $TRIM_DIR/9_CTRL/ $TRIM_DIR/9_CTRL/9*_R1_val_1.fq $TRIM_DIR/9_CTRL/9*_R2_val_2.fq

kallisto quant -i Mus_musculus.idx -t 2 -o $TRIM_DIR/8_KO/ $TRIM_DIR/8_KO/8*_R1_val_1.fq $TRIM_DIR/8_KO/8*_R2_val_2.fq

