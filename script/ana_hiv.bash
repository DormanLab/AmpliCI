#!/bin/bash

## changing pho
PHO=20.01
FILENAME=SRR2241783_2

## barcode data
./run_AmpliCI --umi -f ../test/$FILENAME.trim.noN.bc.fastq -o ../test/$FILENAME.bc

## Estimate error profile with partitions
grep "assign" ../test/$FILENAME.bc.out | awk '{$1="";print}' > ../test/partition.txt
./run_AmpliCI -f ../test/$FILENAME.trim.noN.fastq -o ../test/error.out --partition ../test/partition.txt  --exclude --abundance 2 --error

## merged data
./run_AmpliCI -f ../test/SRR2241783_2.merge.noN.fastq -o ../test/$FILENAME.merge -p ../test/error.out -trim 9
cat ../test/$FILENAME.merge.fa |  seqkit subseq -r 10:250 | seqkit rmdup -s | seqkit seq -w 0 > ../test/$FILENAME.merge.trim_dedup.fa

## run main algorithm
./run_AmpliCI -f ../test/SRR2241783_2.merge.noN.fastq -u ../test/$FILENAME.bc.fa -i ../test/$FILENAME.merge.trim_dedup.fa -x 9 -o ../test/test_MPLE_${PHO}_band0_amplici_merge_epar_modified_collision_exc2 -p ../test/error.out -rho $PHO
