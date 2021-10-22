#!/bin/bash

## changing pho
PHO=20.01
FILENAME=SRR2241783_2

## barcode data
./run_AmpliCI --umi --fastq ../test/$FILENAME.trim.noN.bc.fastq --outfile ../test/$FILENAME.bc

## Estimate error profile with partitions
grep "assign" ../test/$FILENAME.bc.out | awk '{$1="";print}' > ../test/partition.txt
./run_AmpliCI --fastq ../test/$FILENAME.trim.noN.fastq --outfile ../test/error.out --partition ../test/partition.txt  --exclude --abundance 2 --error

## merged data
./run_AmpliCI --fastq ../test/SRR2241783_2.merge.noN.fastq --outfile ../test/$FILENAME.merge --profile ../test/error.out -trim 9
cat ../test/$FILENAME.merge.fa |  seqkit subseq -r 10:250 | seqkit rmdup -s | seqkit seq -w 0 > ../test/$FILENAME.merge.trim_dedup.fa

## run main algorithm
./run_AmpliCI --fastq ../test/SRR2241783_2.merge.noN.fastq --umifile ../test/$FILENAME.bc.fa --haplotype ../test/$FILENAME.merge.trim_dedup.fa --umilen 9 --outfile ../test/test_MPLE_${PHO}_band0_amplici_merge_epar_modified_collision_exc2 --profile ../test/error.out -rho $PHO
