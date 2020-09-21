/**
 * @file options.c
 * @author Karin S. Dorman
 *
 * Command-line options handling.
 *
 * Note about formatting.  Line widths are at 80 characters, not because we live
 * in the 60's but to help force good coding and to reduce complexity.  Function
 * predeclarations may break this rule so that the entire prototype can be
 * found with a simple grep on the source code.
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>



#include "options.h"
#include "cmdline.h"
#include "amplici.h"

/**
 * Setup options object.
 */
int make_options(options **opt) {
	options *op;
	*opt = malloc(sizeof **opt);

	if (*opt == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION, "options object");

	op = *opt;

	op->info = SILENT; //DEBUG_III;//DEBUG_I;  // control DEBUG message 
	op->use_curses = 0;
	op->wp = NULL;
	op->active_fp = NULL;


	/* ampliCI */
	op->fastq_file = NULL;
	op->offset_file = NULL;
	op->outfile_base = NULL;
	op->outfile_fasta = NULL;
	op->outfile_info = NULL;
	op->initialization_file = NULL;
	op->trans_matrix = NULL;
	op->partition_file = NULL;

	op->run_amplici =ALGORITHM_AMPLICI; 
	op->low_bound = 2.0;
	op->contamination_threshold = 1;
	op->associate_zc = 0; /* contamination_threshold = low_bound -1 */
	op->deletion_error = 1;

	/* model */
	op->JC69_model = 1;
	op->convergence_amplici = 1;
	op->check_false_positive = 1;
	op->use_aic = 0;
	op->nw_align = ALIGNMENT_HAPLOTYPES ; //ALIGNMENT_UNIQ_SEQ; //ALIGNMENT_HAPLOTYPES; // NO_ALIGNMENT
	op->indel_model = INDEL_PER_READ;  // consider a indel model
	
	
	/* K  number of clusters */
	op->K_max = 10000;
	op->K = 0;	/* invalid value */
	op->estimate_K =1;
	op->K_space=100;
	op->K_fix_err = 0;   
	op->filter_reads = 0; 


	/* error profile estimation */
	op->error_estimation = 0;  // estimate error profile in current run, default is 0
	op->min_cosdist = log(0.999);
	op->use_error_profile = 0;	
	op->err_encoding = XY_ENCODING;   // Indicate that how error encoded in error profiles STD_ENCODING or XY_ENCODING
	op->error_profile_name = NULL; 
	
	
	/* error */ 
	op->insertion_error = 0.00004;   
	op->deletion_error = 0.00002;	
	op->indel_error = op->insertion_error + op->deletion_error;   // per site per read 
	op->indel_error_set = 0;
	
	
	/* running criterion */
	op->n_iter_amplici = 1000;
	op->ll_cutoff = -100.0;
	op->epsilon_aln = -100; // Currently just as the epsilon for the alignment result.
	op->epsilon = 1e-6;
	op->alpha = 0.001;
	op->p_threshold = 1e-40;   /* will be overwritten by alpha / M or alpha */
	op->per_candidate = 1;
	op->most_abundant = 0;

	/* for sequence alignment */
	op->score[0][0] = 2; op->score[0][1] = -3; op->score[0][2] = -3;op->score[0][3] = -2;
	op->score[1][0] = -3; op->score[1][1] = 2; op->score[1][2] = -2;op->score[1][3] = -3;
	op->score[2][0] = -3; op->score[2][1] = -2; op->score[2][2] = 2;op->score[2][3] = -3;
	op->score[3][0] = -2; op->score[3][1] = -3; op->score[3][2] = -3;op->score[3][3] = 2;
	op->gap_p = -5;
	op->band =  20;  // maybe need to be changed ? change from 10 to 20
	op->ends_free = 0;  //default 0
	op->max_offset = 2;


	/* UMI information */
	op->UMI_length = 0;
	op->initialization_UMI = NULL;
	op->K_UMI = 0;
	op->topN = 10;
	op->trans_penalty = MPLE;
	op->rho = 1.01;
	op->omega = 1e-20;
	op->threshold_UMI = 1;    // allowed minimal UMI abundance 
	op->threshold_hap = 0;    // allowed minimal deduplicated abundance of haplotypes

	return NO_ERROR;
} /* make_options */

/**
 * Free options object.
 */
void free_options(options *opt)
{
	if (opt)
		free(opt);
} /* free_options */

/**
 * Parse command-line.
 */
int parse_options(options *opt, int argc, const char **argv)
{
	int i, j;
	size_t n;
	int err = NO_ERROR;
	char a;

	for (i = 1; i < argc; i++) {
		n = strlen(argv[i]);
		if (n < 2)
			usage_error(argv, i, (void *)opt);
		j = 0;
		while ((a = argv[i][++j]) == '-' && j < (int) n);
		switch(a) {
		case 'a':
			if (!strncmp(&argv[i][j], "ali", 3)) {
				opt->nw_align = ALIGNMENT_UNIQ_SEQ;
				break;
			}
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strncmp(&argv[i][j], "abu", 3)) {
				opt->low_bound = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} else {
				opt->alpha = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			}
			break;
		case 'c':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else {
				if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57)
					opt->contamination_threshold
						= read_uint(argc, argv, ++i,
								(void *)opt);
				opt->associate_zc = 0;
			}
			break;
		case 'z':
			opt->nw_align = ALIGNMENT_UNIQ_SEQ;
			break;
		case 'n':
			if (!strcmp(&argv[i][j], "nJC69")) {
				opt->JC69_model = 0;
			}else{
				opt->nw_align = NO_ALIGNMENT;
			}
			break;
		case 'k':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strcmp(&argv[i][j], "kmax")) {
				if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57)
					opt->K_max = read_uint(argc,
						argv, ++i, (void *)opt);
			} else if (argv[i+1][0] >= 48
					&& argv[i+1][0] <= 57) {
					opt->K = read_uint(argc, argv, ++i,
						(void *)opt);
					opt->estimate_K = 0;  // if K is provided, not estimate K anymore, 
					opt->K_fix_err = 1;   // fix the maximum number of clusters, K to estimated error profile.
				//	opt->run_amplici = 0;  // not run_amplici to select K
			} 
			if (errno)
				goto CMDLINE_ERROR;
			break;
		case 'r':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}else{
				opt->rho = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			}
			break;
		case 'l':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strcmp(&argv[i][j], "lb")) {
				opt->low_bound = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} else if (!strcmp(&argv[i][j], "ll")
					|| !strncmp(&argv[i][j], "log", 3)) {
				opt->ll_cutoff = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} 
			break;
		case 'e':
		/* KSD: Overriding XP, but making back compatible. */
			if (i + 1 < argc && !strcmp(argv[i + 1], "error")) {
				i++;
				opt->error_estimation = 1;
			} else {
				opt->error_estimation = 1;
			}
			break;
		case 'd':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strncmp(&argv[i][j], "dia", 3)) {
				opt->alpha = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} else if (!strncmp(&argv[i][j], "del", 3)) {
				if (opt->indel_error_set) {
					err = mmessage(INFO_MSG,
						INVALID_USER_INPUT, "Cannot "
						"set --deletion and --indel\n");
					goto CMDLINE_ERROR;
				}
				opt->deletion_error = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			}
			break;
		case 'i':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strcmp(argv[i+1], "amplici")) {
				opt->run_amplici = ALGORITHM_AMPLICI;
				++i;
			} else if (!strncmp(&argv[i][j], "ind", 3)) {
				opt->indel_error = read_cmdline_double(argc,
							argv, ++i, (void *)opt);
				opt->indel_error_set = 1;
				opt->insertion_error = opt->deletion_error
							= opt->indel_error / 2;
			} else if (!strncmp(&argv[i][j], "ins", 3)) {
				if (opt->indel_error_set) {
					err = mmessage(INFO_MSG,
						INVALID_USER_INPUT, "Cannot "
						"set --insertion and --indel\n");
					goto CMDLINE_ERROR;
				}
				opt->insertion_error = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} else {	/* deprecated: use --haplotype */
				opt->initialization_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Haplotype set: "
					"%s\n", opt->initialization_file);
				opt->run_amplici = 0; 
			}
			if (errno)
				goto CMDLINE_ERROR;
			break;
		case 'u':
			if (!strcmp(&argv[i][j], "umi")) {  /* Parameter set sepcific for clustering UMIs */
				opt->gap_p = -20;
				opt->band = opt->max_offset;   // need further investigation. currently set as 0.
				//opt->ends_free = 0;   // not counting the offset the begining 
				// opt->nw_align = NO_ALIGNMENT;
				opt->JC69_model = 0;
				mmessage(INFO_MSG, NO_ERROR, "Cluster UMIs .... \n");
			}else if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}else{
				opt->initialization_UMI = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "UMI set: "
					"%s\n", opt->initialization_UMI);
			}
			break;
		case 'x':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (argv[i+1][0] >= 48
					&& argv[i+1][0] <= 57) {
				opt->UMI_length =  read_uint(argc, argv, ++i,
						(void *)opt);
			}
			if (errno)
				goto CMDLINE_ERROR;
			break;
		case 'f':
			if (!strcmp(&argv[i][j], "filter")) {
				opt->filter_reads = 1;
				i++;
			}else if (i == argc - 1){
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}else{
				opt->fastq_file = argv[++i];
			}
			break;
		case 'm':	/* hidden option: --most */
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}
			opt->most_abundant = atoi(argv[++i]);
			break;
		case 'o':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}

			/* single option argument */
			if (i + 1 == argc - 1 || argv[i + 2][0] == '-') {
				opt->outfile_base = argv[++i];

				/* check for file extension */
				j = 0;
				n = strlen(argv[i]);
				while ((a = argv[i][j++]) != '.' && j < n);

				/* guess the user has named a fasta file */
				if (j < n && (!strcmp(&argv[i][j], "fsa")
					|| !strcmp(&argv[i][j], "fasta")
					|| !strcmp(&argv[i][j], "ffn")
					|| !strcmp(&argv[i][j], "frn")
					|| !strcmp(&argv[i][j], "faa")
					|| !strcmp(&argv[i][j], "fna")
					|| !strcmp(&argv[i][j], "fa")
					|| !strcmp(&argv[i][j], "fas"))) {
					opt->outfile_fasta = argv[i];
					opt->outfile_base = NULL;
				}
			/* >1 option arguments */
			} else {
				opt->outfile_fasta = argv[++i];
				opt->outfile_info = argv[++i];
			}
			break;
		case 'p':
			if (!strncmp(&argv[i][j], "per", 3)) {
				opt->per_candidate = 1;
				break;
			}
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}
			if (!strcmp(&argv[i][j], "partition")){
				opt->partition_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Partition file: "
					"%s\n", opt->partition_file);
			}else if (!strcmp(&argv[i][j], "pdiag")) {
				opt->alpha = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
				opt->per_candidate = 0;
			} else {
				opt->error_profile_name = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Error profile: "
					"%s\n", opt->error_profile_name);
			}
			break;
		case 'w':
			opt->use_curses = 1;
			break;
		case 's':
			if (i + 3 >= argc) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}
			int as = read_int(argc, argv, ++i, (void *)opt);
			opt->score[0][0] = opt->score[1][1] = opt->score[2][2]
							= opt->score[3][3] = as;
			as = read_int(argc, argv, ++i, (void *)opt);
			opt->score[0][3] = opt->score[3][0] = opt->score[1][2]
							= opt->score[2][1] = as;
			if (i + 2 < argc && strlen(argv[i+2]) > 1
				&& argv[i + 2][1] >= 49 && argv[i + 2][1] <= 57)
				as = read_int(argc, argv, ++i, (void *)opt);
			opt->score[0][1] = opt->score[1][0]
				= opt->score[0][2] = opt->score[2][0]
				= opt->score[1][3] = opt->score[3][1]
				= opt->score[2][3] = opt->score[3][2] = as;
			as = read_int(argc, argv, ++i, (void *)opt);
			opt->gap_p = as;
			mmessage(INFO_MSG, NO_ERROR, "Alignment parameters: "
				"match = %d, transition = %d, transversion = %d"
				", gap = %d\n", opt->score[0][0],
				opt->score[0][3], opt->score[0][1], opt->gap_p);
			break;
		case 't':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else {
				opt->trans_matrix = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Output transition matrix: "
					"%s\n", opt->trans_matrix);
			}
			break;
		case 'h':
			if (!strncmp(&argv[i][j], "hap", 3)) {
				if (i == argc - 1) {
					err = INVALID_CMD_OPTION;
					goto CMDLINE_ERROR;
				}
				opt->initialization_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Haplotype set: "
					"%s\n", opt->initialization_file);
				opt->run_amplici = 0; 
			} else {
				fprint_usage(stderr, argv[0], opt);
				free_options(opt);
				exit(EXIT_SUCCESS);
			}
			break;
		default:
			err = INVALID_CMD_OPTION;
			goto CMDLINE_ERROR;
		}
	}


	/* check compatibility of options */
	if (!opt->fastq_file)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"No input fastq file! (see help -h) \n");

	if (!opt->outfile_base && !opt->most_abundant && !opt->outfile_fasta)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"No output file specified! (see help -h) \n");

	if (opt->error_profile_name && opt->error_estimation)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Cannot combine options --error and --profile\n");

	if (opt->error_estimation && !opt->outfile_base)
		opt->outfile_base = opt->outfile_fasta;

	if (opt->initialization_file && !opt->outfile_base)
		opt->outfile_base = opt->outfile_fasta;

	//if(opt->initialization_file && opt->estimate_K)
	//	err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
	//		"Please provide number of haplotypes (K) in your haplotype set (-k)\n");

	/* If the user want to find a true K */
	/* [KSD] AmpliCI always estimates true K? Nope [XY] if K is provided, opt->estimate_K = 0 */
	if (opt->estimate_K)
		opt->K = opt->K_space;
	else
		opt->K_max = opt->K;
	

	if (opt->error_profile_name)
		opt->use_error_profile = 1;

	/* update indel error rates */
	if (!opt->indel_error_set)
		opt->indel_error = opt->insertion_error + opt->deletion_error;


	/* check contamination threshold */
	if (opt->associate_zc) /* no input for contamination threshold. Use the default */
		opt->contamination_threshold = (unsigned int) opt->low_bound - 1;

	if(opt->contamination_threshold > opt->low_bound){
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Contamination threshold should be set below low_bound \n");
	}


	return err;

CMDLINE_ERROR:
	if (err == NO_ERROR) {
		err = INVALID_CMD_ARGUMENT;
		i--;
	}
	usage_error(argv, i, (void *)opt);
	return err;
} /* parse_options */

/**
 * Usage statement.
 */
void fprint_usage(FILE *fp, const char *cmdname, void *obj)
{
	options *opt = (options *) obj;
	size_t start = strlen(cmdname) - 1;

	while (cmdname[start] != '/' && start) start--;
	if (cmdname[start] == '/') start++;

	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	//fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "(v1.0.1)\n");
	fprintf(fp, "\nNAME\n\t%s - Amplicon Clustering Inference\n",
		&cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s [--error | --profile <pstr>]"
		" [--abundance <adbl> --log_likelihood <lldbl> --diagnostic <ddbl> --per_candidate --align --scores <smint> <smmint> <sgint>]"
		" [--haplotypes <hstr>]"
		" --fastq <fstr> --outfile <ostr>\n",
		&cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\tProgram %s clusters amplicon sequences presented "
		"as reads in fastq file <fstr> in two steps.\n"
		"\n\tStep 1:\n\t\tEstimate error profile and write to file <ostr> given input fastq file <fstr> with option -e error."
		"\n\t\tExample: \n\t\t\trun_ampliCI --error --fastq <fstr> --outfile <ostr>"
		"\n\tStep 2:\n\t\tCluster amplicon sequences in fastq file <fstr> with estimated error profile from file <pstr>, using lower bound 2, and output results in file <ostr>."
		"\n\t\tExample: \n\t\t\trun_ampliCI --profile <pstr> --abundance <adbl> --fastqs <fstr> --outfile <ostr>"

		"\n\n\tBy default, %s will automatically select K true haplotypes with estimated scaled "
		"true abundance greater than or equal to <adbl> using an alignment-free strategy. "
		"A Diagnostic test is used to screen false positives with provided threshold.  The "
		"default diagnostic test threshold is quite liberal.  "
		"Reads with log likelihood higher than <lldbl> under the current model are assigned "
		"to clusters.\n", &cmdname[start], &cmdname[start]);
	fprintf(fp, "\n\n\t%s can also be used to assign reads to user-provided haplotypes by using the "
		"--haplotypes option.  This can be helpful if there are known haplotypes in the sample, "
		" but it is also useful for careful abundance estimation.", &cmdname[start]);
		
	/* KSD: Maybe we should make AmpliCI-cons the default. */
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t--abundance | -lb <adbl>\n\t\tLower bound for scaled true abundance during haplotype reconstruction.  [DEFAULT: %f]\n", opt->low_bound);
	fprintf(fp, "\t--align | -z \n\t\tAlign all reads to haplotypes (slow).  [DEFAULT: no]\n");	 /* KSD:  --align | -a */
	fprintf(fp, "\t--contaminants | -c <ctuint>\n\t\tBaseline count abundance of contaminating or noise sequences.  [DEFAULT: %i]\n", opt->contamination_threshold);
	fprintf(fp, "\t--deletion <deldbl>\n\t\tSequencing deletion error rate.  [DEFAULT: %f]\n", opt->deletion_error);
	fprintf(fp, "\t--diagnostic | -a <ddbl>\n\t\tThreshold for diagnostic probability in the diagnostic/contamination test.  [DEFAULT: %f].\n", opt->alpha);
	fprintf(fp, "\t--error | -e\n\t\tEstimate the error profile.\n");
	fprintf(fp, "\t--fastq | -f <fstr>\n\t\tThe fastq input file.  [REQUIRED]\n");
	fprintf(fp, "\t--haplotypes | -i <hstr>\n\t\tFASTA file with haplotypes.  [DEFAULT: none]\n");
	fprintf(fp, "\t--help | -h\n\t\tThis help.\n");
	fprintf(fp, "\t--indel <inddbl>\n\t\tSequencing indel error rate.  Cannot also use options --insertion or --deletion.  [DEFAULT: %f]\n", opt->indel_error);
	fprintf(fp, "\t--insertion <insdbl>\n\t\tSequencing insertion error rate.  [DEFAULT: %f]\n", opt->insertion_error);
	fprintf(fp, "\t--kmax <kuint>\n\t\tSet maximum number of clusters K.  [DEFAULT: %i]\n", opt->K_max);
	fprintf(fp, "\t--log_likelihood | -ll <lldbl>\n\t\tLower bound for screening reads during cluster assignment.  This is the minimum log assignment likelihood, ln pi_k + ln Pr(r_i|h_k). [DEFAULT: %f]\n", opt->ll_cutoff);
	fprintf(fp, "\t--outfile | -o <ostr>|<ostr1> <ostr2>\n\t\tOutput file(s) for haplotype discovery, estimated error profile (when used with --error), or cluster assignments (when used with --haplotypes).  [REQUIRED]\n");
	fprintf(fp, "\t\tBy default, provide base name of file to output haplotypes (extension .fa) and information (extension .out) or provide names for both files, FASTA first.\n");
	fprintf(fp, "\t\tWith --error, provide name of file to output error profile.\n");
	fprintf(fp, "\t\tWith --haplotypes, provide name of file to output read cluster assignments.\n");
	fprintf(fp, "\t--per_candidate | --pdiag <pdbl>\n\t\tAdjust diagnostic threshold (--diagnostic) to %f / number_candidates.  [DEFAULT: %s]\n", opt->alpha, opt->per_candidate ? "yes" : "no");
	fprintf(fp, "\t--profile | -p <estr>\n\t\tThe input error profile. If none, convert quality score to Phred error probability.  [DEFAULT: none]\n");
	fprintf(fp, "\t--scores <match> <mismatch> [<transversion_mismatch>] <gap>\n\t\tSet scores of the Needleman-Wunsch aligner.  [DEFAULT: %d %d %d %d]\n", opt->score[0][0], opt->score[0][3], opt->score[0][1], opt->gap_p);
//	fprintf(fp, "\t-k <kuint>\n\t\tNumber of haplotypes in the haplotype set (used with -i <hstr>).  [DEFAULT: %i]\n", opt->K);	/* KSD: get rid of this option */
	/* fprintf(fp, "\t--most <mint>\n\t\tReport top m-most abundant sequences and quit. [DEFAULT: %i]\n", opt->most_abundant); */
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */
