#!/bin/bash


## bc data
./run_AmpliCI --umi --fastq ../test/sim2.bc.fq --outfile ../test/sim2.bc

## Estimate error profile with partitions
grep "assign" ../test/sim2.bc.out | awk '{$1="";print}' > ../test/partition.txt
./run_AmpliCI --fastq ../test/sim2.trim.fq  --outfile ../test/error.out --partition ../test/partition.txt -exclude --abundance 2 --error

## merged data
./run_AmpliCI --fastq ../test/sim2.fq --outfile ../test/sim2.merge --profile ../test/error.out -trim 9 
cat ../test/sim2.merge.fa |  seqkit subseq -r 10:250 | seqkit rmdup -s | seqkit seq -w 0 > ../test/sim2.merge.trim_dedup.fa

## run main algorithm
./run_AmpliCI --fastq ../test/sim2.fq -u ../test/sim2.bc.fa -i ../test/sim2.merge.trim_dedup.fa -x 9 --outfile ../test/test --profile ../test/error.out -rho 46
