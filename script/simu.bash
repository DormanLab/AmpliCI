#!/bin/bash

DATA_FOLDER=~/Documents/Research/data/hiv/src/simulation/new
UMI_FASTA=test_barcode_nalign_lq.fa
HAPLOTYPE_FASTA=ASV.collect.l2.fa

## changing parameters
FILENAME=sim8.5
CYCLES=10
EFFICIENCY=0.5

mkdir "$DATA_FOLDER/$FILENAME"

python3 simulation.py -b $UMI_FASTA -i $HAPLOTYPE_FASTA -o $DATA_FOLDER/$FILENAME/$FILENAME -c $CYCLES -e $EFFICIENCY

/Users/xiyupeng/Documents/Research/analysis/art_simulator/art_src_MountRainier_MacOS/art_illumina -ss MSv1 -amp -na -i $DATA_FOLDER/$FILENAME/$FILENAME.fasta -l 250 -f 1 -o $DATA_FOLDER/$FILENAME/$FILENAME -ir 0.00002 -dr 0.00002

##preparing UMI and sequence fastq files.
seqkit subseq -r 1:9 $DATA_FOLDER/$FILENAME/$FILENAME.fq > $DATA_FOLDER/$FILENAME/$FILENAME.bc.fq
seqkit subseq -r 10:250 $DATA_FOLDER/$FILENAME/$FILENAME.fq > $DATA_FOLDER/$FILENAME/$FILENAME.trim.fq

##preparing true haplotype set
head -n 55 $FILENAME.para.txt | tail -50 > $FILENAME.true.fa
