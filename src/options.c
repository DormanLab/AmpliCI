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

	op->info = SILENT;  // control DEBUG message 
	op->use_curses = 0;
	op->wp = NULL;
	op->active_fp = NULL;


	/* ampliCI */
	op->fastq_file = NULL;
	op->offset_file = NULL;
	op->outfile = NULL;
	op->initialization_file = NULL;
	op->run_amplici =ALGORITHM_AMPLICI; 
	op->low_bound = 2.0;
	op->contamination_threshold = 1;
	op->associate_zc = 1; /* contamination_threshold = low_bound -1 */
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

	
	/* error profile estimation */
	op->error_estimation = 0;  // estimate error profile in current run, default is 0
	op->min_cosdist = log(0.999);
	op->use_error_profile = 1;	/* [KSD] bug setup: initialized mismatched */
	op->err_encoding = XY_ENCODING;   // Indicate that how error encoded in error profiles STD_ENCODING or XY_ENCODING
	op->error_profile_name = NULL; 
	
	
	/* error */ 
	op->insertion_error = 0.00004;   
	op->deletion_error = 0.00002;	
	op->indel_error = op->insertion_error + op->deletion_error;   // per site per read 
	
	
	/* running criterion */
	op->n_iter_amplici = 1000;
	op->ll_cutoff = -100.0;
	op->epsilon_aln = -100; // Currently just as the epsilon for the alignment result.
	op->epsilon = 1e-6;
	op->alpha = 0.001;
	op->p_threshold = 1.;   /* default = 1, or 1e-40 (conserved); depends on alpha */
	op->most_abundant = 0;

	/* for sequence alignment */
	op->score[0][0] = 2; op->score[0][1] = -3; op->score[0][2] = -3;op->score[0][3] = -2;
	op->score[1][0] = -3; op->score[1][1] = 2; op->score[1][2] = -2;op->score[1][3] = -3;
	op->score[2][0] = -3; op->score[2][1] = -2; op->score[2][2] = 2;op->score[2][3] = -3;
	op->score[3][0] = -2; op->score[3][1] = -3; op->score[3][2] = -3;op->score[3][3] = 2;
	op->gap_p = -5;
	op->band =  20;  // maybe need to be changed ? change from 10 to 20

	return NO_ERROR;
} /* make_options */

/**
 * Free options object.
 */
void free_options(options *opt)
{
	if(opt) {
		free(opt);
	}
} /* free_options */

/**
 * Parse command-line.
 */
int parse_options(options *opt, int argc, const char **argv)
{
	int i, j;
	int err = NO_ERROR;
	char a;

	for (i = 1; i < argc; i++) {
		if (strlen(argv[i]) < 2)
			usage_error(argv, i, (void *)opt);
		j = 1;
		a = argv[i][j];
		while (a == '-' && ++j < (int) strlen(argv[i]))
			a = argv[i][j];
		switch(a) {
		case 'a':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else {
				opt->alpha = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			}
			break;
		case 'c':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strcmp(&argv[i][j], "cont")){
				if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57)
				opt->contamination_threshold = read_uint(argc,
						argv, ++i, (void *)opt);
				opt->associate_zc = 0;
			}
			break;
		case 'z':
			opt->nw_align = ALIGNMENT_UNIQ_SEQ;
			break;
	//	case 'n':  // not used anymore, moved
	//		opt->estimate_K = 0;  // default change to 1, use -n to turn off
	//		break;
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
				//	opt->run_amplici = 0;  // not run_amplici to select K
			} 
			if (errno)
				goto CMDLINE_ERROR;
			break;
		case 'l':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strcmp(&argv[i][j], "lb")) {
				opt->low_bound = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} else if (!strcmp(&argv[i][j], "ll")) {
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
			} else if (!strcmp(&argv[i][j], "deletion")) {
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
			} else if (!strncmp(&argv[i][j], "ins", 3)) {
				opt->insertion_error = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} else {
				opt->initialization_file = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Haplotype set: "
					"%s\n", opt->initialization_file);
					opt->run_amplici = 0; 
				}
			if (errno)
				goto CMDLINE_ERROR;
			break;
		case 'f':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}
			opt->fastq_file = argv[++i];
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
			opt->outfile = argv[++i];
			break;
		case 'p':
			if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			} else if (!strcmp(&argv[i][j], "pdiag")) {
				opt->p_threshold = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
			} else {
				opt->error_profile_name = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Error profile: "
					"%s\n", opt->error_profile_name);
			}
			break;
		case 'w':
			opt->use_curses = 1;
			break;
		case 'h':
			fprint_usage(stderr, argv[0], opt);
			free_options(opt);
			exit(EXIT_SUCCESS);
		default:
			err = INVALID_CMD_OPTION;
			goto CMDLINE_ERROR;
		}
	}


	/* check compatibility of options */
	if (!opt->fastq_file)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"No input fastq file! (see help -h) \n");

	if (!opt->outfile && !opt->most_abundant)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"No output file specified! (see help -h) \n");

	if (opt->error_profile_name && opt->error_estimation)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Cannot combine options --error and --profile\n");

	if(opt->initialization_file && opt->estimate_K)
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Please provide number of haplotypes (K) in your haplotype set (-k)\n");

	/* If the user want to find a true K */
	/* [KSD] AmpliCI always estimates true K? Nope */
	if (opt->estimate_K)
		opt->K = opt->K_space;
	else
		opt->K_max = opt->K;
	

	if (!opt->error_profile_name)
		opt->use_error_profile = 0;

	/* update indel error rates */
	opt->indel_error = opt->insertion_error + opt->deletion_error;


	/* check contamination threshold */
	if(opt->associate_zc) /* no input for contamination threshold. Use the default */
		opt->contamination_threshold = (unsigned int) opt->low_bound - 1;

	if (opt->contamination_threshold < 1){
		err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
			"Contamination threshold could not be set below 1\n");
	}else if(opt->contamination_threshold >= opt->low_bound){
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
	fprintf(fp, "(v1.0.0)\n");
	fprintf(fp, "\nNAME\n\t%s - Amplicon Clustering Inference\n",
		&cmdname[start]);
	fprintf(fp, "\nSYNOPSIS\n\t%s [--error --profile <estr> "
		"-lb <lbdbl> -ll <lldbl> -a <adbl> --pdiag <pdbl> -i <hstr> -k <kuint> -z ] "
		"-f <fstr> -o <ostr>\n",
		&cmdname[start]);
	fprintf(fp, "\nDESCRIPTION\n\tProgram %s clusters amplicon sequences presented "
		"as reads in fastq file <fstr> in two steps.\n"
		"\n\tStep 1:\n\t\tEstimate error profile and write to file <ostr> given input fastq file <fstr> with option -e error."
		"\n\t\tExample: \n\t\t\trun_ampliCI --error -f <fstr> -o <ostr>"
		"\n\tStep 2:\n\t\tCluster amplicon sequences in fastq file <fstr> with estimated error profile from file <estr>, using lower bound 2, and output results in file <ostr>."
		"\n\t\tExample: \n\t\t\trun_ampliCI --profile <estr> -lb 2 -f <fstr> -o <ostr>"

		"\n\n\tBy default, %s will automatically select K true haplotypes with estimated scaled "
		"true abundance greater than or equal to <lbdbl> using an alignment-free strategy. "
		"A Diagnostic test is used to screen false positives with provided threshold."
		"Reads with log likelihood higher than <lldbl> under the current model are assigned "
		"to clusters.\n", &cmdname[start], &cmdname[start]);
		
	fprintf(fp, "\nOPTIONS\n");
	fprintf(fp, "\t--error | -e\n\t\tEstimate the error profile.\n");
	fprintf(fp, "\t-f <fstr>\n\t\tThe fastq input file. [REQUIRED]\n");
	fprintf(fp, "\t--profile | -p <estr>\n\t\tThe input error profile. If none, treat quality score literally. [DEFAULT: none]\n");
	fprintf(fp, "\t-a <adbl>\n\t\tThreshold of diagnostic probability in the diagnostic test will be given as [DEFAULT: %f/number_candidates].\n", opt->alpha);
	fprintf(fp, "\t--pdiag <pdbl>\n\t\tOverride option -a and give a hard threshold of diagnostic probability.\n");
	fprintf(fp, "\t-i <hstr> \n\t\tThe fasta input haplotype set. [DEFAULT:none]\n");	/* KSD: why keep this option; there is no other possibility, right? XY: I will do reads assignment with given haplotype set in the future*/
	fprintf(fp, "\t-k <kuint>\n\t\tNumber of haplotypes in the haplotype set (used with -i <hstr>). [DEFAULT: %i]\n", opt->K);
	fprintf(fp, "\t--kmax <kuint>\n\t\tSet maximum number of clusters K. [DEFAULT: %i]\n", opt->K_max);
	fprintf(fp, "\t-lb <lbdbl>\n\t\tLower bound for scaled true abundance during haplotype reconstruction.  [DEFAULT: %f]\n", opt->low_bound);
	fprintf(fp, "\t-cont <ctuint>\n\t\t baseline count abundance of contaminating or noise sequences.  [DEFAULT: %i]\n",opt->contamination_threshold);
	fprintf(fp, "\t-ll <lldbl>\n\t\tLower bound for reads maximum posterior assignment probability screening during reads assignment. [DEFAULT: %f]\n", opt->ll_cutoff);
	fprintf(fp, "\t-insertion <indbl>\n\t\tInsertion error rate.  [DEFAULT: %f]\n", opt->insertion_error);
	fprintf(fp, "\t-deletion <dedbl>\n\t\tDeletion error rate.  [DEFAULT: %f]\n", opt->deletion_error);
	fprintf(fp, "\t-o <ostr>\n\t\tOutput file to record best clustering solution or the estimated error profile.  [REQUIRED]\n");
	/* fprintf(fp, "\t--most <mint>\n\t\tReport top m-most abundant sequences and quit. [DEFAULT: %i]\n", opt->most_abundant); */
	fprintf(fp, "\t-z \n\t\tTurn off the alignment-free.  [DEFAULT: none]\n");	 /* KSD: And turn on what? If not alignment-free, what do you get? XY: pairwise alignment of each reads and haplotypes*/
	fprintf(fp, "\t-h\n\t\tThis help.\n");
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(cmdname); ++i)
		fputc(toupper(cmdname[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */
