DAUMI
=====

DAUMI, Amplicon Clustering Inference with UMI information, is part of [AmpliCI](https://github.com/DormanLab/AmpliCI) v2.0+.
It takes advantage of UMI information to denoise amplicon sequencing data.
The UMIs contribute to a more accurate error profile and better clustering results.
DAUMI greatly enhances the accuracy of detecting rare sequences and provides deduplicated abundance estimation, eliminating sequencing errors, many PCR errors and PCR bias.

# Table of Contents
1. [Installation](#installation)
1. [Preparing input](#input)
	1. [Paired-end reads](#paired)
1. [Usage and tutorial](#usage)
1. [Choosing rho](#parameter)
1. [Output files](#output)
1. [Command-line options](#options)
1. [Acknowledgments](#acknowledgments)
1. [Citation](#citation)
1. [Data used in DAUMI publication](#data)
1. [Contact](#contact)


# Installation <a name = "installation" />

DAUMI is part of the AmpliCI software. Please check the [AmpliCI](https://github.com/DormanLab/AmpliCI) instructions for installation.

In addition, you may find the following software useful to implement the DAUMI pipeline:

- Automatic selection of command-line option [rho](#parameter) requires the [Fastest Fourier Transform of the West](https://tlo.mit.edu/technologies/fftw-fastest-fourier-transform-west) library.
- Processing of the FASTQ files in the pipeline may be easiest with [seqkit](https://github.com/shenwei356/seqkit).

We use both software packages in the tutorial and demonstrations below.


# Preparing input <a name="input" />

The input to AmpliCI is a FASTQ file. 
Besides the common practice for [preprocessing](https://github.com/DormanLab/AmpliCI#input) raw FASTQ files, additional preprocessing steps are required for DAUMI.
For each sample, users need to prepare three FASTQ files:

- `FILENAME.umi.fq`
- `FILENAME.trim.fq`
- `FILENAME.fq`


Each FASTQ file contains subsequences of the original FASTQ file: UMIs (`FILENAME.umi.fq`), biological sequences (`FILENAME.trim.fq`) or UMIs + biological sequences (`FILENAME.fq`).
Subsequences from the same read share the same seq ID and are in the same order in all three files.
In `FILENAME.fq`, UMIs should be at the front of biological sequences.
It is easy to prepare the three FASTQ files using [seqkit](https://github.com/shenwei356/seqkit) with subcommands `subseq` and `concat`.

Below is an example demonstrating how to prepare the other two files from `FILENAME.fq`, when the total read length is 250nt with a UMI length of 9nt.

```
seqkit subseq -r 1:9 FILENAME.fq > FILENAME.umi.fq
seqkit subseq -r 10:250 FILENAME.fq > FILENAME.trim.fq
```

or generate `FILENAME.fq` from the other two files

```
seqkit concat FILENAME.umi.fq FILENAME.trim.fq > FILENAME.fq 
```

## Paired-end reads <a name = "paired" />

If you have paired-end reads, you may simply concatenate them and follow the regular [Usage pipeline](#usage).
Here we justify this advice.

If your read pairs overlap, then we would not recommend merging overlapping reads prior to analysis, as most read merging tools do not properly update the quality scores, which we rely on to detect sequencing errors.
A better and easier solution with UMI-tagged paired-end reads is to simply to concatenate the reads, even if they overlap, and treat the concatenated reads as a single sampled sequence.
It is true that some of the base calls are reads of the same true nucleotide, and we are not using the full information available in these replicate reads to estimate the original molecule, but there is little harm done in ignoring the information.
If there is a PCR error in the overlap region, then two sites will register that change, which could lead to more false positives than a post-merge solution.
However, our method (and other methods) only accidentally handle PCR errors, so it seems foolhardy to recommend a step (merging reads) that is likely to disrupt the signal we do model (sequencing errors) to partially overcome a signal we do not model (PCR errors).

If your paired-end reads do not overlapped, then we would recommend the same concatenation strategy without the caveats.
The error model does not utilize the read position to predict errors or assume any dependence in the errors between sites, so there is no problem with concatenating reads, where you will lose information about the sequencing cycle of each nucleotide.



# Usage <a name="usage" />

DAUMI consists of three major steps.
1. Obtain candidate UMI sequences. We use `AmpliCI cluster` for this purpose. You may use other software to obtain a FASTA file with candidate UMIs.
1. Obtain candidate haplotype sequences (sample sequences without the UMI) in two steps using AmpliCI. You may use other software to obtain a FASTA file with candidate haplotypes.
	1. Estimate an error profile using `AmpliCI error`.
	1. Obtain candidate haplotype sequences using `AmpliCI cluster`.
1. Estimate the true haplotypes and their abundance using `AmpliCI daumi`.

We provide more details and demonstrate each of these steps below.

1. Cluster UMI sequences (the executable is called `run_AmpliCI`):

	```sh
	./run_AmpliCI cluster --umi --fastq FILENAME.umi.fq --outfile FILENAME.umi
	```

	The input of this command is the FASTQ file with UMIs only from [preparing data](#input).
	The output of this command are the found UMIs in `FILENAME.umi.fa` and the clusterings, among other information, in `FILENAME.umi.out`.
	You will need both files moving forward.
	This command can fail if there is no replication in the UMIs.
	If your data falls in this scenario then UMI-based deduplication will fail (or you have misidentified the UMIs).

	An example (from the ```src``` directory):

	```sh
	./run_AmpliCI cluster --umi --fastq ../test/sim2.umi.fq --outfile ../test/sim2.umi
	```

1. Obtain candidate haplotype sequences.
	1. Estimate error profile based on partition induced by umi clusters:

		```sh
		./run_AmpliCI error --fastq FILENAME.trim.fq --outfile FILENAME.trim.err --partition FILENAME.umi.out --exclude
		```

		The input of this command is the FASTQ file without UMIs from [preparing data](#input).
		The output of this command is the error profile stored in `FILENAME.trim.err` or filename of your choice.
		This command will fail if there are no replicated reads in the input file.
		AmpliCI relies on data replication to work.
		It assumes PCR amplification and low error rates guarantee that every true UMI and molecule is observed at least twice in the input file.

		An example:

		```sh
		./run_AmpliCI error --fastq ../test/sim2.trim.fq  --outfile ../test/sim2.trim.err --partition ../test/sim2.umi.out --exclude
		```

	1.  Cluster umi-tagged sequences to get initial haplotype set ([seqkit](https://github.com/shenwei356/seqkit) required for truncation and deduplication):

		```sh
		./run_AmpliCI cluster --fastq FILENAME.fq --outfile FILENAME --profile FILENAME.trim.err --trim UMI_LENGTH --nJC69
		cat FILENAME.fa | seqkit subseq -r START_IDX:END_IDX | seqkit rmdup -s | seqkit seq -w 0 > FILENAME.trim.fa
		```

		The input of this command are the FASTQ file `FILENAME.fq` with UMI and sampled read from [preparing data](#input) and the error profile `FILENAME.trim.err` from step 1.
		The output of the AmpliCI command are the estimated candidate haplotypes with UMIs in file `FILENAME.fa` and additional information in file `FILENAME.out`.
		The second step shown above removes the UMI tags and the resulting duplicated haplotypes from the candidate haplotypes in `FILENAME.fa`.
		The software [seqkit](https://github.com/shenwei356/seqkit) is useful for this purpose.
		`START_IDX` and `END_IDX` are start and end position of the sampled biological sequences, 1-based and inclusive.
		In other wrods, the UMI length is `START_IDX - 1` and the biological sequence (candidate haplotype) length is `END_IDX - START_IDX + 1`.
		The final output of this step is the candidate haplotype file `FILENAME.trim.fa`. All other intermediate files from this step can be discarded.
		Note: Compared to the published paper, we now recommend to add additional option --nJC69, which slightly increases the precision of the result.

		An example (total read length 250nt with UMI length 9nt):

		```sh
		./run_AmpliCI cluster --fastq ../test/sim2.fq --outfile ../test/sim2 --profile ../test/sim2.trim.err --trim 9 --nJC69
		cat ../test/sim2.fa | seqkit subseq -r 10:250 | seqkit rmdup -s | seqkit seq -w 0 > ../test/sim2.trim.fa
		```

1. Estimate deduplicated abundance of each haplotype. We describe how to [select parameter rho](#parameter) in the following section.

	```sh
	./run_AmpliCI daumi --fastq FILENAME.fq --umifile FILENAME.umi.fa --haplotype FILENAME.trim.fa --profile FILENAME.trim.err --rho RHO -umilen UMI_LENGTH --outfile FILENAME
	```

	The input of this command is the FASTQ file `FILENAME.fq` with UMI and sampled read from [preparing data](#input), the candidate UMI file `FILENAME.umi.fa` from Step 1, the error profile `FILENAME.trim.err` and candidate haplotype file `FILENAME.trim.fa` from Step 2, and the choice of rho (see [Choosing Rho](#parameter)).
	The output of this command is described in [Output Files](#output).

	An example:

	```sh
	./run_AmpliCI daumi --fastq ../test/sim2.fq --umifile ../test/sim2.umi.fa -haplotype ../test/sim2.trim.fa -umilen 9 --outfile ../test/sim2 --profile ../test/sim2.trim.err -rho 46
	```

You can run the whole pipeline on this example (from the ```src``` directory):

```sh
bash ../script/ana.bash
```

There are many possible modifications of this pipeline.
One possible alternative to the second step for initializing the haplotype set is to directly run AmpliCI on non-UMI-tagged FASTQ files.

```sh
./run_AmpliCI --fastq FILENAME.trim.fq --outfile FILENAME.trim --profile FILENAME.err
```
The alternative step costs less time but risks missing some very similar haplotypes.
We recommend running the alternative step on massive datasets with sequences of moderate similarity, e.g., 16S rRNA gene sequences.

You can always get more help with:

```sh
./run_AmpliCI --help
```

# Choosing Rho <a name = "parameter" />

Rho should be chosen as the expected number of copies of a molecule in the dataset that best separates error molecules from true molecules.
One can examine the histogram of observed UMI abundances to choose a good value of Rho.
In the cleanest situation, there will be a peak near 1 consisting of the errored UMIs and another well-separated peak (or multiple peaks) representing the true UMIs that have been amplified and sequenced multiple times.
Unfortunately, several forces collaborate to smear and overlap the two distributions of errored and true UMIs.
In particular, PCR amplification is not 100% perfect, so not every molecule is copied in every PCR cycle.
Furthermore, only a subset of amplified molecules will be sampled for sequencing and some of those will be sequenced with error.
These forces reduce the number of observed copies of true molecules.
Furthermore, errors during PCR amplification can amplify errored molecules to high levels, making some errored molecules look like true molecules.
Together, these forces cause the two distributions to overlap making it difficult to choose a threshold for some datasets.

We have written a separate program to help you make a good choice for parameter rho given the UMI abundance histogram.
Be warned that this program cannot do magic, and if the errored and true abundance distributions are highly overlapped, just like you, it will have a hard time finding a good threshold.
You can find the program in the ```script``` directory.
The following demonstration shows how to compile and run the script.

1. **Code Preparation**

First, you need to install the [Fastest Fourier Transform of the West](https://tlo.mit.edu/technologies/fftw-fastest-fourier-transform-west).
Then, the necessary code can be compiled with
```
cd ../script
gcc -Wall -pedantic -Wextra -g -o bp_pmf_mix bp_pmf_mixture.c fft.c -lfftw3 -lm
```

2. **Data Preparation**

Run the histogram command available in AmpliCI to obtain the UMI abundance table.
```
./run_AmpliCI histogram --fastq FILENAME.umi.fq --outfile FILENAME.umi.hist
```

3. **Model Fitting**

The input file is the UMI raw abundance distribution in the sample, that is the number in the ith row is the number of unique UMIs with abundance i.
You can find examples of this input file with `*.txt` extension under the ```script``` directory.
The first step is to fit our proposed model (details in the paper) to data (the executable is called `bp_pmf_mix`).
```
./bp_pmf_mix -f FILENAME.umi.hist -t TRUNCATION_POSITION > OUPUT_CSV_FILE
```
You can use option `-t` to select a truncation position to truncate the long right tail of the distribution, which may be contributed by unmodeled UMI collision.

4. **Refit Model**

You can also use the estimated values of parameters as the input to the program, in order to obtain the desired distribution without reestimating parameters.
Below we use the input file from a HIV dataset as an example. 
```
./bp_pmf_mix -f hiv1_raw_abun.txt -t 100 --efficiency 0.6 --ncycles 11 --epsilon 0.004 --delta 0.01 > hiv1.csv
```
Options `--efficiency`, `--ncycles`, `--epsilon` and `--delta` are four parameters in the model, that represent PCR efficiency, number of PCR cycles, PCR error rate and sequence error rate. You can find more details about the model in the paper. 

5. **Select Rho**

We recommend choosing rho as the argmax {Pr(X <= rho | Z = 1) < 0.05}, which means at most five percent of true variants will have observed abundance below or equal to rho.
You can find this probability in the 8th column of the output csv file.
Examine the lines in the output file until you find the last line with 8th column less than 0.05.
Then, the corresponding rho is found in the first column.

You can find more information about selecting rho in the [publication](#citation).

# Output Files <a name = "output" />

The final step of the pipeline with argument `--outfile FILENAME` generates two output files:

***1. `FILENAME.fa`***

FASTA-formatted file containing denoised sequences (or haplotypes) as well as their deduplicated abundances. 
For example, for the first haplotype, the FASTA header might look like:

```
>H0;Deduplicated Abundance=18.000;
```

For each haplotype, deduplicated abundance is the estimated number of unamplified sample sequences.


***2.  `FILENAME.out`***

A text file with the following information provided as key: value pairs, one per line.  The keys are:

- `log likelihood`: Penalized log likelihood of the model. See the paper for more details.

- `K`: Number of haplotypes (total number of sequences in the input haplotype set <input_hap_fasta_file> )

- `assignments`: model-assigned haplotype for each read in FASTQ-determined input order.  Haplotypes are numbered 0, 1, ..., and match the sequences H0, H1, ... in the output FASTA file of haplotypes.

- `cluster sizes`: Number of reads assigned to each haplotype.

- `UMI K`: Number of unique UMIs (total number of sequences in the input UMI set <input_umi_fasta_file> )

- `UMI assignment`: model-assigned UMI for each read in FASTQ-determined input order.

- `UMI cluster sizes`: Number of reads assigned to each UMI.

- `reads ll`: The posterior log likelihood for each read.

- `Eta`: Relative abundance of each unique UMI

- `Gamma`: A UMI K X K matrix used to indicate dependence between UMI and haplotypes. The deduplicated abundance of the kth haplotype is the kth column sum of the Gamma.

- There is also a list of haplotypes with related UMIs reported in this file.

The example below showed one haplotype with related 3 UMIs. Thus the deduplicated abundance of the haplotype is 3.  

```
AGGATTGATTAAATATTATTGTCCTATTGAAGTGTTCTCTCAATTTTTCACTTACTCTTTGTAAAGTTTTCTCCCATTTAGTTTTATTAATGTTACAATGTGCTTGTCTTATATCTCCTATTATGTCTCCTGTTGCATAGAATGTTTGTCCTGGTCCTATTCTTACACTTTTTCTTGTATTATTGTTGGGTCTTATACAATTAATCTCTACAGATTCATTGAGATGTACTATTATTGTTTT
TTTTAAAAC GTAAATAGT CGCTAATGA
```

# Options <a name = "options" />

Options of AmpliCI can be found in [here](https://github.com/DormanLab/AmpliCI#options). Below we list options specific for DAUMI.

- `--partition`: Use a partition file when estimating errors. The partition file contains cluster assignment of each read. Reads assigned with the same number are in the same group. We can use either true partition or UMI-induced partition file to generate a better error profile. 

- `--exclude`: Exclude small clusters during error estimation (set threshold with option --abundance).

- `--abundance` under `error` estimation mode: Lower bound on observed abundance for inclusion of seeded cluster during error estimation.

- `--umi`: Used for clustering UMIs. Compared to the default, we set the gap score -20, and band width 2, used in Needleman Welch alignment. We also disable JC69 model since UMIs are random sequences.

- `--trim`:Ignore first # nucleotides in JC69 model. Since UMIs are random sequences, they should be ignored when fit the JC69 model.

- `--ncollision`: Assume NO UMI collision, that same UMI CANNOT be attached to two different original haplotypes. This option can also be used when estimate errors based on UMI-induced partition file when there is no UMI collision.

- `--umifile`: FASTA file with UMIs.

- `--rho`:  Tuning parameter that control the sparsity of the transition matrix `Gamma`. We have described how to select `rho` above.

- `--umilen`: UMI Length.
 

# Acknowledgments <a name = "acknowledgments" />

- See acknowledgments of [AmpliCI](https://github.com/DormanLab/AmpliCI#acknowledgements)


# Citation <a name = "citation" />

- Peng, X. and Dorman, K. S. (2020) 'AmpliCI: A High-resolution Model-Based Approach for Denoising Illumina Amplicon Data', Bioinformatics. doi: [10.1093/bioinformatics/btaa648](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btaa648/5875058).

- Peng, X. and Dorman, K. S. (2022) 'Accurate estimation of molecular counts from amplicon sequence data with unique molecular identifiers', *BioRxiv*. doi: [10.1101/2022.06.12.495839v1](https://www.biorxiv.org/content/10.1101/2022.06.12.495839v1). 

- Peng, X. and Dorman, K. S. (2023) 'Accurate estimation of molecular counts from amplicon sequence data with unique molecular identifiers', *Bioinformatics*. [Advanced access](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btad002/6971842).

# Data used in DAUMI publication <a name = "data" />

- [HIV-1 V1V3 Data, SRR2241783](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2241783) 
- [HIV-1 V3 Data, SRR5105420](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5105420)
- [scRNA-seq Data, E-MTAB-10372](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10372/sdrf)

# Contact <a name = "contact" />

If you have any problems, please contact:

Xiyu Peng (pansypeng124@gmail.com)

Karin Dorman (kdorman@iastate.edu)
