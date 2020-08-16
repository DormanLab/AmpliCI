AmpliCI
=======

AmpliCI, Amplicon Clustering Inference, denoises Illumina amplicon data by approximate model-based clustering.

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Installation](#installation)
1. [Preparing input](#input)
	1. [Demultiplexing](#demultiplexing)
	1. [Quality control and read processing](#quality)
	1. [Input files](#inputfiles)
1. [Usage](#usage)
1. [Output](#output)
1. [Downstream analysis](#downstream)
1. [Troubleshooting](#troubleshooting)
1. [Detailed options](#options)
1. [C library](#library)
1. [Acknowledgements](#acknowledgements)
1. [Citation](#citation)
1. [Contact](#contact)

# Prerequisites <a name = "preresuisites" />

- AmpliCI requires [cmake](https://cmake.org) (3.5.0 or higher version) and [gcc](https://gcc.gnu.org) (5.4.0 or higher version).
- AmpliCI requires some C and FORTRAN libraries provided by R.  You can download and install R from [https://www.r-project.org](https://www.r-project.org).
- AmpliCI requires Rmath, the [R Standalone Math Library](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).  Often, the Rmath library (libRmath.a or libRmath.so for Linux or libRmath.dylib for MacOS) will be installed with R, but not always.  Here are some other locations for the library.
	- r-mathlib on [Ubuntu](https://ubuntu.com/) and [Debian](https://www.debian.org/)
	- libRmath on [Fedora](https://ubuntu.com/), [CentOS](https://centos.org/), [Mageia](https://www.mageia.org/en/), and [Mandriva](https://www.openmandriva.org/)
	- Or if all else fails, you can install the Rmath standalone library from the repository [https://github.com/statslabs/rmath](https://github.com/statslabs/rmath)

# Installation <a name = "installation" />

AmpliCI has been tested under Linux and MacOS.

1. Clone the repository.

    ```sh
    git clone https://github.com/DormanLab/AmpliCI.git
    ```

2. Configure the project.

   ```sh
   cd AmpliCI/src
   cmake .
   ```

3. Compile AmpliCI.  The executable is called ```run_AmpliCI```.  It will appear in the ```src``` directory you are currently in.

   ```sh
   make
   ```

# Preparing input <a name="input" />

The input of AmpliCI is a FASTQ file, but there is some necessary preprocessing.

## **Demultiplexing** <a name="demultiplexing" />

Like all other denoising methods, the starting point of the analysis is FASTQ sequence data after demultiplexing.  If you start with separate barcode and read FASTQ files, you can use the qiime script [split\_libraries\_fastq.py](http://qiime.org/scripts/split_libraries_fastq.html) for demultiplexing.  Use the option ```--store_demultiplexed_fastq``` to keep demultiplexed fastq files.

## **Quality control and read preprocessing** <a name="quality" />

AmpliCI requires all input reads have the **same** length, with no **ambiguous** nucleotides (only A, C, G, T base calls allowed).  (One way to truncate or filter reads with ambiguous nucleotides is via the R package [ShortRead](https://rdrr.io/bioc/ShortRead/).)

## **Input files** <a name="inputfiles" />

AmpliCI takes a single demultiplexed FASTQ file (one per sample) generated from the Illumina sequencing platform, with reads trimmed to the same length and containing no ambiguous nucleotides (see above steps).  If you have paired end data, AmpliCI can analyze the forward reads, the reverse reads, or the merged reads, but not both forward and reverse reads simultaneously.

You can find example input FASTQ files in the [test](https://github.com/DormanLab/AmpliCI/tree/master/test) directory.

One read in the input FASTQ file should fit in exactly **four** lines, as in the format below.

```
@SRR2990088.351
TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTTTTAAGTCAGCGGTGAAAGTCTGTGGCTCAACCATAGAATTGCCGTTGAAACTGGGAGGCTTGAGTATGTTTGAGGCAGGCGGAATGCGTGGTGTAGCGGTGAAATGCGTAGATATCAAGCAGAACACCGATTGCGAAGGCAGCTTGCTAAGCCATGACTGACGCTGATGCACGAAAGCGTGGGGATGAAACA
+
CCCCCGGGG8CFCFGGEGGGGGGGGGB@FFEEFFGFCFFFGGGGGGGEFGGG9@@F@FF9EFFG<EEGD@EFFGGGG,ECBCEFGCAFEFEEF<E?FEFFG<F@FFFGGG9FG@FGGG8DEGGGD,A=4,AEDF+F3BCCEEE7DFCGEEDEFEGFEGEGE<@@F>*:?BB7@;,>,5,*;CC:,4C957*:AB5<=DF6:>/*5*121/(/*500.<<;52(444+164-83::>:B021;91-(.<6).
```

First line: @sequence name

Second line: DNA sequence (A, T, C, G)

Third line: +any content on a single line

Fourth line: quality score sequence ([ASCII](https://en.wikipedia.org/wiki/FASTQ_format#Encoding) [!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ])

If your read or quality scores are split over multiple lines, AmpliCI will not work.  One possible script for fixing your FASTQ-formatted files is given by [Damian Kao on BioStars](https://www.biostars.org/p/14828/).

# Usage

AmpliCI runs in two major steps:

1. Use AmpliCI to estimate the error profile directly from the data (the executable is called run_AmpliCI):

```sh
./run_AmpliCI --fastq <input_fastq_file> --outfile <output_error_profile_file> --error
```

An example (from the ```src``` directory):

```sh
./run_AmpliCI --fastq ../test/sim3.8.1.fastq --outfile ../test/error.out  --error
```

2.  Use Amplici to estimate the haplotypes and their abundance using the estimated error profile:

```sh
./run_AmpliCI --fastq <input_fastq_file> --outfile <output_base_filename> --abundance 2 --profile <input_error_profile_file>
```

An example (from the ```src``` directory):

```sh
./run_AmpliCI --fastq ../test/sim3.8.1.fastq --outfile ../test/test --abundance 2 --profile ../test/error.out
```

If you provide no input error profile with the `--profile` option, AmpliCI will assume the error rates are the error rates dictated by [Phred quality scores](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf).
Assuming Phred quality scores is not a good idea.
Using Phred quality scores tends to generate high numbers of false positives and runs very slowly.

- You can also use AmpliCI to reassign reads given input haplotypes. You can provide your haplotype set with '--haplotypes' option:

```sh
./run_AmpliCI --fastq <input_fastq_file> --outfile <output_assignment_filename> --profile <input_error_profile_file> --haplotypes <input_haplotypes_fasta_file>
```

An example (from the ```src``` directory):

```sh
./run_AmpliCI --fastq ../test/sim3.8.1.fastq --outfile ../test/test.id --profile ../test/error.out --haplotypes ../test/test.fa
```


- Detailed help can be obtained with:

```sh
./run_AmpliCI --help
```

- If you apply AmpliCI on longer reads with length > 300 (like merged reads), you may want to decrease the default Lower bound for screening reads during cluster assignment with `--log_likelihood` [DEFAULT: -100.000000]. For example, you can set the lower bound at -200.

```sh
./run_AmpliCI --fastq ../test/sim3.8.1.fastq --outfile ../test/test.id --profile ../test/error.out --haplotypes ../test/test.fa --log_likelihood -200
```


# Output Files <a name = "output" />

## Estimating error profile.

When run to estimate the error profile, AmpliCI will output an error profile `<output_error_profile_file>` in text format.  This is simply a list of comma-separated probabilities (times 1000) of the probability that haplotype nucleotide `n` is misread as read nucleotide `m` with quality score `q`.  They are ordered as `(n,m,q)`, with the last index varying the fastest. Both haplotype nucleotide `n` and read nucleotide `m` are in the order (A,C,T,G) , and `q` has the range from 0 to 40 (41 in total). For example, the first 41 entries are estimated transition probabilities for A->A when observed quality score q is in [0:40]; Then the 42nd - 82nd entries are estimated transition probabilities for A->C; the 165th - 205th entries are estimated transition probabilities for C->A.... In our recommended workflow the error profile should only be used with the dataset from which it was estimated. If you apply AmpliCI estimates to other datasets or use other estimates with AmpliCI, you should consider the following:

- AmpliCI encodes nucleotides in the order of (A, C, T, G), which is different from the commonly used alphabetic order (A, C, G, T).

- Not all quality scores will be observed, especially for quality scores < 3. To avoid extrapolation, error rates for quality scores outside the range of observed quality scores will not be estimated by LOESS regression. Instead, we just assume the error rates dictated by [Phred quality scores](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf) with equal probability of each possible nucleotide substitution when there is an error.

## Estimating haplotypes.

When run to estimate haplotypes and their abundances with argument `--outfile <output_base_filename>` or `--outfile <fasta_output_file> <information_output_file>`, there will be two output files:

***1.`output_base_filename.fa` or `fasta_output_file`***

FASTA-formatted file (will be used in the downstream analysis) containing denoised sequences (or haplotypes).  For each sequence, we also provide `size` (scaled true abundance), `DiagP` (diagnostic probability), `ee` (mean expected number of errors in reads), useful for chimera detection and *post hoc* filtering. For example for the first haplotype, the FASTA header might look like:

```
>H0;size=516.000;DiagP=0.00e+00;ee=0.405;
```

- `size`: scaled true abundance (expected number of error free reads) estimated for each selected haplotype, required for the subsequent chimera detection with UCHIME3.

- `DiagP`: diagnostic probability, which could be used as a criterion to check false positives. We suggest post hoc removal of haplotypes with `DiagP` > 1e-40 when applying AmpliCI on real datasets with more than 1 million reads to reduce false positives. The diagnostic probability may contain an allowance for contaminating sequences (see option `--contaminants`). For further information of the diagnostic probability and contamination screening, please see [our paper](https://www.biorxiv.org/content/10.1101/2020.02.23.961227v1).

- `ee`: mean expected number of errors per read. Edgar and Flyvbjerg ([Edgar and Flyvbjerg, 2015](https://academic.oup.com/bioinformatics/article/31/21/3476/194979)) suggested a strategy to filter reads according to their expected number of errors.  For example, you could remove haplotypes with `ee` > 1. Though this strategy works for some mock datasets, we have observed `ee` > 1 for several true haplotypes with very low abundance when analyzing a specific mock dataset (stag1, see [our paper](https://www.biorxiv.org/content/10.1101/2020.02.23.961227v1)). You can read [more about `ee`](https://www.drive5.com/usearch/manual/exp_errs.html).


***2.`output_base_filename.out` or `information_output_file`***

A text file with the following information provided as key: value pairs, one per line.  The keys are:

- `K`: Number of haplotypes selected by AmpliCI.

- `assignments`: AmpliCI-assigned haplotype by posterior probability for each read in FASTQ-determined input order.  Haplotypes are numbered 0, 1, ..., and match the sequences H0, H1, ... in the output FASTA file of haplotypes.  NA is output if the read's maximum conditional log likelihood (given the source haplotype) does not exceed a user-defined threshold (option `-ll`; default -100).  These assignments are not based on alignment of reads to the haplotypes, so some reads, particularly indel errors, may not be assigned (NA).  See option `--haplotypes` for more careful read assignment.

- `cluster sizes`: Number of reads assigned to each haplotype.

- `pi`: Estimated $\boldsymbol{\pi}$ from AmpliCI.  Each read is assigned to a haplotype by maximum transition probability $\Pr(r_i|h_k)$ (distinct from posterior probability used for assignments) and $\pi_k$ is the proportion of reads assigned to haplotype $k$.

- `reads ll`: For each read, the maximum conditional log likelihood (given the source haplotype), $\ln \pi_k + \ln \Pr(r_i|h_k)$.

- There is also a fasta listing of the haplotypes reported in this file.

- `ee`: For each read, the mean expected number of errors. See discussion on `ee` in ***`output_base_filename.fa`*** above.

- `uniq seq id`: The index of each selected haplotype in the unique sequence list, ordered from highest abundance to lowest.  If the haplotypes were selected in observed abundance order, then these will be increasing integers from 0.  If any unique sequence was discarded, some integers will be skipped.  For example, this line is `0   1   2   3   4   5   6   7  10  42  45` for test file `test/sim3.8.1.fastq`, indicating that the first 8 most observed sequences were selected as haplotypes, but the 9th and 10th most observed sequences were discarded, and so on.

- `scaled true abun`: The estimated scaled true abundances of each selected haplotype (expected number of error free reads).

- `obser abun`: The observed abundance of each selected haplotype.

- `Estimated common ancestor` (value on next two lines in FASTA format): The final, estimated common ancestor of all the haplotypes used in the BIC calculation.

- `Evolution_rate`: The estimated evolutionary time separating each haplotype from the ancestor.

- `log likelihood from JC69 model`: The log likelihood of the JC69 hierarchical model computed on the final, fitted model.

- `Diagnostic Probability threshold`: The threshold used to reject candidate haplotypes in the contamination test.  This is the value input through option `--diagnostic` divided by the number of possible candidate haplotypes.

- `aic`: The estimated [Akaike Information Criterion](https://en.wikipedia.org/wiki/Akaike_information_criterion) value from the final fitted model.

- `bic`: The estimated [Bayesian Information Criterion](https://en.wikipedia.org/wiki/Bayesian_information_criterion) value from the final fitted model.

When run with option `--haplotypes` to reassign reads to the user-provided **haplotype set** (a FASTA-formatted file), AmpliCI will output a read assignment file `<output_assignment_filename>` in text format. The keys are

- `assignments`: See the description above for outfile `output_base_filename.out`.  There should be fewer NA assignments because reads with low log likelihood are aligned to the haplotypes to detect indel sequencing errors.

- `cluster sizes`: See the description above for outfile `output_base_filename.out`.  The sizes should be higher if more reads are successfully assigned to haplotypes.

# Downstream Analysis <a name="downstream" />

The output FASTA file contains denoised raw haplotype sequences, which may include chimeric sequences.  The first step of any downstream analysis should be to remove chimeric sequences.

## **Chimera Detection** <a name="chimera" />

The AmpliCI-outputted FASTA file is in acceptable format to input into the [uchime3_denovo](https://www.drive5.com/usearch/manual/cmd_uchime3_denovo.html) method implemented in [usearch](https://drive5.com/usearch/).

Haplotype sorting by abundance
```sh
./usearch -sortbysize <input_fasta_file> -fastaout <output_sorted_fasta_file>
```

Chimera detection
```sh
./usearch -uchime3_denovo <input_sorted_fasta_file> -uchimeout <uchime_outfile> -chimeras <chimera_fasta_outfile> -nonchimeras <nonchimera_fasta_outfile>
```

You may also use other chimera detection algorithms to remove chimeras.

## **Generate Amplicon Sequence Variant (ASV or sOTU) Table** <a name = "otu" />

We have provided an [R script](https://github.com/DormanLab/AmpliCI/tree/master/script/Make_ASV_Tables.R) to help to generate the ASV (sOTU) table, where scaled true abundances (see `size`) per sample per ASVs/sOTUs are reported.

## **Taxa Assignment** <a name = "taxa" />

You may want to identify the non-chimeric haplotypes detected in your sample.
There are multiple methods.

1. [DECIPHER](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html) contains [IDTAXA](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0521-5), a novel approach for taxonomic classification.

2. [RDP classifier](http://rdp.cme.msu.edu), a Naive Bayesian Classifier.

## **Futher Analysis**	<a name = "further" />

[mothur](https://www.mothur.org), [qiime2](https://qiime2.org), [LefSE](https://huttenhower.sph.harvard.edu/galaxy/), [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html), ....


# Troubleshooting <a name = "troubleshooting" />

The algorithm may stop if your:

- quality scores are not in the typical range for Illumina datasets [33,73]

- reads contain ambiguous nucleotides

- reads are not in the right FASTQ input format, for example reads and quality scores cannot contain newline characters

- there are too few reads or reads are so noisy that there are no sequences observed more than the lower bound number of times (option `--abundance`, default 2.0)

# Detailed options <a name = "options" />

Main options:

- `--fastq` The fastq input file.  [REQUIRED]

- `--outfile` Output file(s) for haplotype discovery, estimated error profile (when used with --error), or cluster assignments (when used with --haplotypes).  [REQUIRED]

- `--profile` The input error profile. If none, convert quality score to Phred error probability.  [DEFAULT: none]

- `--error` Estimate the error profile. [Used in error estimation only]

- `--haplotypes` FASTA file with haplotypes. [Used in reads assignment only]

Options for sensitivity:

- `--abundance` Lower bound for scaled true abundance during haplotype reconstruction (should be >= 2.0).  [DEFAULT: 2.0]

- `--contaminants` Baseline count abundance of contaminating or noise sequences.  [DEFAULT: 1]

- `--indel` Indel sequencing error rate.  Cannot also use options --insertion or --deletion.  [DEFAULT: 0.00006]

- `--diagnostic`  Threshold of diagnostic probability in the diagnostic/contamination test.  [DEFAULT: 0.001 / number_candidates]

Other important options:

- `--align`  Align all reads to haplotypes (slow).  [DEFAULT: none]

- `--log_likelihood`  Lower bound for screening reads during cluster assignment.  This is the minimum log assignment likelihood, $\ln \pi_k + \ln \Pr(r_i|h_k)$. [DEFAULT: -100.000000]


# C library <a name = "library" />

AmpliCI provides a static C library for users to call function ```amplici_wfile()``` to cluster amplicon sequences from another program. The library `libamplici.a` will appear in the ```src``` directory when you compile AmpliCI.

**Input**

- `fastq_file`: The fastq input file. [REQUIRED]

- `error_profile_name`: The input error profile. If `NULL`, convert quality score to Phred error probability.

- `low_bound`: Allowed lowest abundance. See the description of option `--abundance`. [REQUIRED]

**Output**

- `seeds`: Estimated haplotypes.

- `seeds_length`: Lengths of Estimated haplotypes.

- `cluster_id`: See the description of `assignments` above for outfile `output_base_filename.out`. Note ```amplici_wfile()``` does not filter reads with maximal conditional log likelihood under the given threshold. Instead, it assigns all reads to its closest haplotypes with the maximum likelihood.

- `cluster sizes`: Number of reads assigned to each haplotype.

- `K`: Number of estimated haplotypes.

- `sample_size`: Number of reads in the fastq input file

- `max_read_length`: Maximum read length l. The kth (in [0,1,2,...K-1]) haplotype starts at seeds[k*l].

- `abun`: See the description of `scaled true abun` above for outfile `output_base_filename.out`.

- `ll`: See the description of `reads ll` above for outfile `output_base_filename.out`.


An example to call function ```amplici_wfile()``` is provided in [example_wfile.c](https://github.com/DormanLab/AmpliCI/tree/master/example_wfile.c). You can compile the source file with the C library libamplici.a (in the ```src``` directory):

```
gcc -o myprog example_wfile.c -lamplici -lRmath -lm -I ./src/ -L ./src/
```

Use -I to provide path to header file of the library libamplici.h and -L to provide path to the library libamplici.a. You may need to add additional path to Rmath library and header files if needed. Note example_wfile.c needs two more header files in the ```src``` directory, which are not required by the library libamplici.a.


# Acknowledgments <a name = "acknowledgements" />

- AmpliCI contains LOESS regression for error estimation, the original file is available at [https://www.netlib.org/a/dloess](https://www.netlib.org/a/dloess).  However, we modified and used related code from R, which derives from the above.

- We used the hash implemented in [uthash.h](https://troydhanson.github.io/uthash/userguide.html).

# Citation <a name = "citation" />

- Peng, X. and Dorman, K. (2020) ‘AmpliCI: A High-resolution Model-Based Approach for Denoising Illumina Amplicon Data’, Bioinformatics. doi: [10.1093/bioinformatics/btaa648](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btaa648/5875058).

# Contact <a name = "contact" />

If you have any problems with AmpliCI, please contact:

Xiyu Peng (xiyupeng@iastate.edu)

Karin Dorman (kdorman@iastate.edu)
