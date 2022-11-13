#!/bin/bash

## Step 1: Prepare candidate UMIs
./run_AmpliCI cluster --umi --fastq ../test/sim2.umi.fq --outfile ../test/sim2.umi

## Step 2(i): Estimate error profile with UMI partitions
./run_AmpliCI error --fastq ../test/sim2.trim.fq  --outfile ../test/sim2.trim.err --partition ../test/sim2.umi.out -exclude

## Step 2(ii): Prepare candidate haplotypes
./run_AmpliCI cluster --fastq ../test/sim2.fq --outfile ../test/sim2 --profile ../test/sim2.trim.err --trim 9 
cat ../test/sim2.fa | seqkit subseq -r 10:250 | seqkit rmdup -s | seqkit seq -w 0 > ../test/sim2.trim.fa

## Step 3: DAUMI! Estimate true haplotypes and their abundance.
./run_AmpliCI daumi --fastq ../test/sim2.fq --umifile ../test/sim2.umi.fa --haplotype ../test/sim2.trim.fa --umilen 9 --outfile ../test/sim2 --profile ../test/sim2.trim.err -rho 46
