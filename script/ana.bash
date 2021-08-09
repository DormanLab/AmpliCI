#!/bin/bash

## changing file name
FILENAME=sim8.2
PHO=40.01

## bc data
/Users/xiyupeng/Documents/Research/github/AmpliCI-UMI/src/run_AmpliCI --umi -f $FILENAME.bc.fq  -o $FILENAME.bc

## Estimate error profile with partitions
grep "assign" $FILENAME.bc.out > tmp.txt
awk '{$1="";print}' tmp.txt >partition.txt
rm tmp.txt
/Users/xiyupeng/Documents/Research/github/AmpliCI-UMI/src/run_AmpliCI -f $FILENAME.trim.fq  -o error.out --partition partition.txt -exclude --abundance 2 -e 

## merged data
/Users/xiyupeng/Documents/Research/github/AmpliCI-UMI/src/run_AmpliCI -f $FILENAME.fq -o $FILENAME.merge -p error.out -trim 9 
cat $FILENAME.merge.fa |  seqkit subseq -r 10:250 | seqkit rmdup -s | seqkit seq -w 0 > $FILENAME.merge.trim_dedup.fa


## run main algorithm
/Users/xiyupeng/Documents/Research/github/AmpliCI-UMI/src/run_AmpliCI -f $FILENAME.fq -u $FILENAME.bc.fa -i $FILENAME.merge.trim_dedup.fa -x 9 -o test_MPLE_${PHO}_band0_amplici_merge_epar_modified_collision_exc2 -p error.out -rho $PHO
