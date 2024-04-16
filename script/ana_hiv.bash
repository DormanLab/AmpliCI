#!/bin/bash

## changing pho
PHO=20.01
FILENAME=SRR2241783_2

## barcode data
./run_AmpliCI cluster --umi --fastq ../test/$FILENAME.trim.noN.bc.fastq --outfile ../test/$FILENAME.bc

## Estimate error profile with partitions
./run_AmpliCI error --fastq ../test/$FILENAME.trim.noN.fastq --outfile ../test/error.out --partition ../test/$FILENAME.bc.out  --exclude --abundance 2 

## merged data
./run_AmpliCI cluster --fastq ../test/$FILENAME.merge.noN.fastq --outfile ../test/$FILENAME.merge --profile ../test/error.out -trim 9 --nJC69
cat ../test/$FILENAME.merge.fa |  seqkit subseq -r 10:250 | seqkit rmdup -s | seqkit seq -w 0 > ../test/$FILENAME.merge.trim_dedup.fa

## run main algorithm
./run_AmpliCI daumi --fastq ../test/$FILENAME.merge.noN.fastq --umifile ../test/$FILENAME.bc.fa --haplotype ../test/$FILENAME.merge.trim_dedup.fa --umilen 9 --outfile ../test/test_MPLE_${PHO}_band0_amplici_merge_epar_modified_collision_exc2 --profile ../test/error.out -rho $PHO
