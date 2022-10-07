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
	op->ignor_nc = 0;
	
	
	/* K  number of clusters */
	op->K_max = 10000;
	op->K = 0;	/* invalid value */
	op->estimate_K =1;
	op->K_space=100;
	op->K_fix_err = 0;   
	op->filter_reads = 0; 


	/* error profile estimation */
	op->error_estimation = 0;  // estimate error profile in current run, default is 0
	op->seed_min_observed_abundance = 2;  // set lower for low deduplicated sample
	op->exclude_low_abundance_seeds = 0;
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
	op->topN = 10;      // should > than 1 and < K_UMI * K
	op->trans_penalty = MPLE;
	op->rho = 1.01;
	op->omega = 1e-20;
	op->threshold_UMI = 1;    // allowed minimal UMI abundance 
	op->threshold_hap = 0;    // allowed minimal deduplicated abundance of haplotypes
	op->umicollision = 1;   // consider UMI collision by default

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
	char const *cmd = "cluster";	/* default command */
	char const *user_cmd = NULL;
	char a;

	for (i = 1; i < argc; i++) {
		/* new command format */
		if (i == 1 && argv[i][0] != '-') {
			cmd = argv[i];
			if (!strcmp(cmd, "cluster")) {
				user_cmd = cmd;
				mmessage(INFO_MSG, NO_ERROR,
							"Command: cluster\n");
			} else if (!strcmp(cmd, "error")) {
				user_cmd = cmd;
				opt->error_estimation = 1;
				mmessage(INFO_MSG, NO_ERROR,
							"Command: error\n");
			} else if (!strcmp(cmd, "assignment")) {
				user_cmd = cmd;
				opt->run_amplici = 0; 
				mmessage(INFO_MSG, NO_ERROR,
						"Command: assignment\n");
			} else if (!strcmp(cmd, "cluster_wumi")){
				user_cmd = cmd;
				opt->run_amplici = 0; 
				mmessage(INFO_MSG, NO_ERROR,
						"Command: cluster_wumi\n");
			}else {
				mmessage(ERROR_MSG, INVALID_CMDLINE,
					"Command '%s' not recognized\n", cmd);
				goto CMDLINE_ERROR;
			}
			continue;
		}
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
				if (!strcmp(cmd, "cluster")) {
					opt->low_bound = read_cmdline_double(
						argc, argv, ++i, (void *)opt);
					mmessage(INFO_MSG, NO_ERROR,
						"Lower bound: %f.\n",
						opt->low_bound);
				} else if (!strcmp(cmd, "error")) {
					opt->seed_min_observed_abundance = 
						strtoul(argv[++i], NULL, 0);
					mmessage(INFO_MSG, NO_ERROR,
						"Minimum abundance: %u.\n",
						opt->seed_min_observed_abundance);
				} else {
					err = INVALID_CMD_OPTION;
					goto CMDLINE_ERROR;
				}
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
		case 'v':
			if (i + 1 < argc && argv[i+1][0] != '-') {
				opt->info = strtoul(argv[++i], NULL, 0);
			} else {
				++opt->info;
			}
			mmessage(INFO_MSG, NO_ERROR, "Verbosity set to %d.\n",
								opt->info);
			break;
		case 'n':
			if (!strcmp(&argv[i][j], "nJC69")) {
				opt->JC69_model = 0;
			}else if (!strncmp(&argv[i][j], "ncol", 4)){
				opt->umicollision = 0;
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
			} else if(!strncmp(&argv[i][j], "ex", 2)) {
				opt->exclude_low_abundance_seeds = 1;
				mmessage(INFO_MSG, NO_ERROR, "Excluding low "
					"abundance seeds (set --abundance).\n");
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
			if (!strcmp(&argv[i][j], "umi")) {  /* Parameter set sepcific for clustering UMIs --umi */
				opt->gap_p = -20;
				opt->band = opt->max_offset;   // need further investigation. 
				//opt->ends_free = 0;   // not counting the offset the begining 
				// opt->nw_align = NO_ALIGNMENT;
				opt->JC69_model = 0;
				opt->use_aic = 1;
				// opt->per_candidate = 0;
				mmessage(INFO_MSG, NO_ERROR, "Cluster UMIs .... \n");
			}else if (i == argc - 1) {
				err = INVALID_CMD_OPTION;
				goto CMDLINE_ERROR;
			}else if(!strncmp(&argv[i][j], "umilen",6)){
				if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57)
					opt->UMI_length =  read_uint(argc, argv, ++i,
						(void *)opt);
			}else{
				opt->initialization_UMI = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "UMI set: "
					"%s\n", opt->initialization_UMI);
			}
			if(errno)
				goto CMDLINE_ERROR;
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
			if (!strcmp(&argv[i][j], "omega")) {
				opt->omega = read_cmdline_double(argc,
					argv, ++i, (void *)opt);
				break; 
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
			}
			if (!strcmp(&argv[i][j], "trim")) {
				if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57)
					opt->ignor_nc = read_uint(argc,
						argv, ++i, (void *)opt);
				mmessage(INFO_MSG, NO_ERROR, "ignore first "
					"%i nucleotides in JC69 model\n", opt->ignor_nc);
			}else if (!strcmp(&argv[i][j], "topT")) {
				if (argv[i + 1][0] >= 48
						&& argv[i + 1][0] <= 57)
					opt->topN = read_uint(argc,
						argv, ++i, (void *)opt);
				if(opt->topN < 1){
					err = mmessage(ERROR_MSG, INVALID_USER_INPUT,
						"topN is set below 1 \n");
					goto CMDLINE_ERROR;
				}
			} else {
				opt->trans_matrix = argv[++i];
				mmessage(INFO_MSG, NO_ERROR, "Output "
						"transition matrix: %s\n",
							opt->trans_matrix);
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
				fprint_usage(stderr, argv[0], user_cmd, opt);
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
void fprint_usage(FILE *fp, const char *exe_name, const char *command, void *obj)
{
	options *opt = (options *) obj;
	size_t start = strlen(exe_name) - 1;


	while (exe_name[start] != '/' && start) start--;
	if (exe_name[start] == '/') start++;

	for (size_t i = start; i < strlen(exe_name); ++i)
		fputc(toupper(exe_name[i]), fp);
	//fprintf(fp, "(%d)\n", 1);
	fprintf(fp, "(v1.99)\n");
	/* default command is to cluster */
	if (!command) {
		command = "cluster";
		fprintf(fp, "\nNAME\n\t%s - Amplicon Clustering Inference\n",
			&exe_name[start]);
		fprintf(fp, "\nSYNOPSIS\n\t%s [--error | --profile <pstr>]"
			" [--abundance <adbl> --log_likelihood <lldbl> --diagnostic <ddbl> --per_candidate --align --scores <smint> <smmint> <sgint>]"
			" [--haplotypes <hstr>] [ --umifile <hstr> --umilen <uuint> --rho <rdbl>] [--umi --nJC69 --nNW]"
			" --fastq <fstr> --outfile <ostr>\n",
			&exe_name[start]);
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
			"to clusters.\n", &exe_name[start], &exe_name[start]);
		fprintf(fp, "\t%s can also be used to assign reads to user-provided haplotypes by using the "
			"--haplotypes option.  This can be helpful if there are known haplotypes in the sample, "
			" but it is also useful for careful abundance estimation.\n", &exe_name[start]);
		fprintf(fp, "\n\tPlease check %s subcommand cluster, error, assignment, "
				"cluster-wumi with option -h for more information.\n",&exe_name[start]);
	} else if (!strcmp(command, "cluster")) {
		fprintf(fp, "\nNAME\n\t%s-%s - Cluster Reads\n", &exe_name[start], command);
		fprintf(fp, "\nSYNOPSIS\n\t%s cluster [--profile FILE] [--abundance FLOAT] [--log_likelihood FLOAT] [--contaminants UINT --diagnostic FLOAT --per_candidate] [--align --scores INT INT INT] [--umi --nJC69 --nNW] --fastq FILE --outfile FILE\n", &exe_name[start]);
		fprintf(fp, "\nDESCRIPTION\n\tCluster reads, finding an estimated K true haplotypes (bounded above if --kmax set) with estimated scaled true abundance greater than a lower bound (--abundance), using, by default, an alignment-free strategy (unless --align). You can provide an error profile (--profile) or use the PHRED defaults, but the latter is likely to produce many false haplotypes. Candidate haplotypes can be screened for possible contamination (--contaminants, --diagnostic, --per_candidate).\n");
		fprintf(fp, "\n\tThough the method is alignment-free, alignments to existing haplotypes are verified when considering candidate haplotypes to rule out rare indel errors. You may set the alignment scoring system (--scores) and the indels probabilities (--insertion, --deletion, --indel).");
		fprintf(fp, "\n\tAfter haplotypes are defined, a final clustering of all reads is attempted. You may limit which reads get clustered in the final clustering by placing a lower bound on the log likelihood (--log_likelihood), for example to leave unclustered likely contaminants.\n");
		fprintf(fp, "\n\tWe recommend to use --umi option for clustering UMIs.\n");
	} else if (!strcmp(command, "error")) {
		fprintf(fp, "\nNAME\n\t%s-%s - Estimate Error Profile\n", &exe_name[start], command);
		fprintf(fp, "\nSYNOPSIS\n\t%s error [--partition FILE] [--abundance INT --exclude] [--ncollision] --fastq FILE --outfile FILE\n", &exe_name[start]);
		fprintf(fp, "\nDESCRIPTION\n\tEstimates the error profile from a FASTQ file, either by alternating AmpliCI clustering with estimation, or given an existing partition of the reads. Since AmpliCI clustering may be used to find the partition, all the options for the cluster command are also active. Type '%s cluster --help' for more information.\n", &exe_name[start]);
		fprintf(fp, "\n\tGiven a partition of the reads, either obtained by AmpliCI clustering or the partition, we fit an error model to the counts of errors across all clusters in the partition. For an externally provided partition, perhaps given by an independent clustering of UMIs, possible collisions within clusters may be detected (NOT detected with --ncollision) by checking for other high abundance members in a cluster and splitting by Hamming distance to the candidate members. A high abundance member reaches a minimum observed abundance and has at least (max. abundance + 1)/2 the abundance of the most abundant member. The user can set the minimum observed abundance threshold using the --abundance option. The user may also choose to exclude clusters that have no member meeting the minimum observed abundance threshold (--exclude). Without replication, it is hard to be sure the cluster represents a real variant or what is the true center of the cluster, so the counts of errors from such clusters may be less reliable. If there is high replication in the dataset, there may be plenty of information to estimate errors from only the clusters around high abundance centers.\n");
		fprintf(fp, "\n\tLike DADA2, we use LOESS regression to estimate the relationship between quality score and each of the error probabilities.\n");
	} else if (!strcmp(command, "assignment")) {
		fprintf(fp, "\nNAME\n\t\t%s-%s -- Assign Reads to Clusters\n", &exe_name[start], command);
		fprintf(fp, "\nSYNOPSIS\n\t%s assignment [--log_likelihood FLOAT] [--profile FILE --nNW] --fastq FILE --haplotype FILE --outfile FILE\n", &exe_name[start]);
		fprintf(fp, "\nDESCRIPTION\n\tAssigns reads to provided haplotypes (--haplotype) with possible lower bound on log likelihood to exclude poorly explained reads. You may use --nNW to avoid time-consuming sequence alignment.\n");
	} else if (!strcmp(command, "cluster_wumi")){
		fprintf(fp, "\nNAME\n\t\t%s-%s -- cluster reads tagged with UMIs\n", &exe_name[start], command);
		fprintf(fp, "\nSYNOPSIS\n\t%s cluster_wumi [--profile FILE] [--ncollision] --fastq FILE --haplotype FILE --umifile FILE --umilen UINT --outfile FILE --rho FLOAT  \n", &exe_name[start]);
		fprintf(fp, "\nDESCRIPTION\n\t Cluster reads tagged with UMIs. Assume UMIs with length (--umilen) are at the beginning of each read in the fastq file. You need to provide two fasta files for initialization, one for haplotypes (--haplotype), one for umis (--umilen). Parameter rho controls the sparsity of transition matrix gamma and needs to be preselected (see the github page for how to select rho). You can provide an error profile (--profile) or use the PHRED defaults.\n");
		fprintf(fp, "\n\tThe algorithm assumes UMI collision by default. The deduplicated abundances of haplotypes would be adjusted if you assume NO UMI collision (--ncollision).\n");
	}
		
	/* KSD: Maybe we should make AmpliCI-cons the default. */
	fprintf(fp, "\n\nOPTIONS\n");
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--abundance | -lb <adbl> \n\t\tLower bound for scaled true abundance during haplotype reconstruction.  [DEFAULT: %f]\n", opt->low_bound);
	else if (!strcmp(command, "error"))
		fprintf(fp, "\t--abundance <adbl> \n\t\tLower bound on observed abundance for inclusion of seeded cluster during error estimation.  [DEFAULT: %u]\n", opt->seed_min_observed_abundance);
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--align | -z \n\t\tAlign all reads to haplotypes (slow).  [DEFAULT: no]\n");	 /* KSD:  --align | -a */
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--contaminants | -c <ctuint>\n\t\tBaseline count abundance of contaminating or noise sequences.  [DEFAULT: %i]\n", opt->contamination_threshold);
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--deletion <deldbl>\n\t\tSequencing deletion error rate (see also --insertion or --indel).  [DEFAULT: %f]\n", opt->deletion_error);
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--diagnostic | -a <ddbl>\n\t\tThreshold for diagnostic probability in the diagnostic/contamination test.  [DEFAULT: %f].\n", opt->alpha);
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--error | -e\n\t\tEstimate the error profile.\n");
	if(!strcmp(command,"error"))
		fprintf(fp, "\t--partition <pstr> \n\t\tPartition file used for a better error profile. [DEFAULT: none] \n");
	if (!strcmp(command, "cluster")){
		//fprintf(fp, "\t--n  \n\t\t Disnable sequence alignment during clustering. Use it when there are no indel errors.  [DEFAULT: no]\n");
		fprintf(fp, "\t--nNW\n\t\tDo NOT use Needleman-Wunsch alignment to align candidate haplotypes to the haplotype set to detect indel errors  [DEFAULT: %s]\n", opt->nw_align  ? "use" : "don't use");
    }	else if(!strcmp(command, "assignment")){
		fprintf(fp, "\t--nNW  \n\t\t Do NOT use Needleman-Wunsch alignment to assign reads. [DEFAULT: %s]\n", opt->nw_align  ? "use" : "don't use");
	}
	if(!strcmp(command, "cluster"))
		fprintf(fp, "\t--umi \n\t\tCluster UMIs. [DEFAULT: no] \n");
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--trim <tuint>\n\t\tIgnore first <tuint> nucleotides in the JC69 model.  [DEFAULT: %i]\n", opt->ignor_nc);
	if (!strcmp(command, "error"))
		fprintf(fp, "\t--exclude\n\t\tExclude small clusters during error estimation (set threshold with option --abundance). [DEFAULT: %s]\n", opt->exclude_low_abundance_seeds ? "yes" : "no");
	fprintf(fp, "\t--fastq | -f <fstr>\n\t\tThe fastq input file.  [REQUIRED]\n");
	if (!strcmp(command, "assignment") || !strcmp(command, "cluster_wumi"))
		fprintf(fp, "\t--haplotypes | -i <hstr>\n\t\tFASTA file with haplotypes.  [REQUIRED]\n");
	if(!strcmp(command, "cluster_wumi"))
		fprintf(fp, "\t--umifile | -u <ustr>\n\t\tFASTA file with UMIs.  [REQUIRED]\n");
	if(!strcmp(command, "cluster_wumi"))
		fprintf(fp, "\t--rho <rdbl> \n\t\t Tunning parameter that control the sparsity of the transition matrix gamma.  [DEFAULT: %f]\n", opt->rho);
	if(!strcmp(command, "cluster_wumi"))
		fprintf(fp, "\t--umilen <uuint> \n\t\t Length of UMI [REQUIRED]\n");
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--indel <inddbl>\n\t\tSequencing indel error rate.  Cannot also use options --insertion or --deletion.  [DEFAULT: %f]\n", opt->indel_error);
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--insertion <insdbl>\n\t\tSequencing insertion error rate (see also --deletion or --indel).  [DEFAULT: %f]\n", opt->insertion_error);
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--kmax <kuint>\n\t\tSet maximum number of clusters K.  [DEFAULT: %i]\n", opt->K_max);
	if (!strcmp(command, "cluster"))
	//fprintf(fp, "\t--nJC69 \n\t\tDisable JC69 model. It should be used when biological sequences are unrelated. [DEFAULT: no]\n");
		fprintf(fp, "\t--nJC69\n\t\tDo NOT use JC69 model for haplotypes. The JC69 model reduces the number of model parameters and increases sensitivity.  [DEFAULT: %s]\n", opt->JC69_model ? "use" : "don't use");
	if (!strcmp(command, "cluster") || !strcmp(command, "assignment"))
		fprintf(fp, "\t--log_likelihood | -ll <lldbl>\n\t\tLower bound for screening reads during cluster assignment.  This is the minimum log assignment likelihood, ln pi_k + ln Pr(r_i|h_k). [DEFAULT: %f]\n", opt->ll_cutoff);
	if (!strcmp(command, "cluster")) {
		fprintf(fp, "\t--outfile, -o FILE | FILE1 FILE2\n\t\tOutput file(s) for haplotype discovery. [REQUIRED]\n");
		fprintf(fp, "\t\tBy default, provide base name of file to output haplotypes (extension .fa) and information (extension .out) or provide names for both files, FASTA first.\n");
	} else if (!strcmp(command, "error")) {
		fprintf(fp, "\t--outfile, -o FILE\n\t\tOutput file for estimated error profile.  [REQUIRED]\n");
	} else if (!strcmp(command, "assignment")) {
		fprintf(fp, "\t--outfile, -o FILE\n\t\tOutput file for cluster assignments.  [REQUIRED]\n");
	}else if (!strcmp(command, "cluster_wumi")){
		fprintf(fp, "\t--outfile, -o FILE \n\t\tOutput file(s) for haplotypes with estimated deduplicated abundance. [REQUIRED]\n");
	}
	if (!strcmp(command, "error"))
		fprintf(fp, "\t--ncollision \n\t\t Assume NO UMI collision when estimating errors based on UMI-induced partition file. [DEFAULT: %s]\n", opt->umicollision ? "collision" : "no collision");
	if (!strcmp(command, "cluster_wumi"))
		fprintf(fp, "\t--ncollision \n\t\t Assume NO UMI collision, that same UMI CANNOT be attached to two different original haplotypes. [DEFAULT: %s]\n", opt->umicollision ? "collision" : "no collision");
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--per_candidate | --pdiag <pdbl>\n\t\tAdjust diagnostic threshold (--diagnostic) to %f / number_candidates.  [DEFAULT: %s]\n", opt->alpha, opt->per_candidate ? "yes" : "no");
	if (!strcmp(command, "cluster") || !strcmp(command, "cluster_wumi") || !strcmp(command, "assignment"))
		fprintf(fp, "\t--profile | -p <estr>\n\t\tThe input error profile. If none, convert quality score to Phred error probability.  [DEFAULT: none]\n");
	if (!strcmp(command, "cluster"))
		fprintf(fp, "\t--scores <match> <mismatch> [<transversion_mismatch>] <gap>\n\t\tSet scores of the Needleman-Wunsch aligner.  [DEFAULT: %d %d %d %d]\n", opt->score[0][0], opt->score[0][3], opt->score[0][1], opt->gap_p);
	fprintf(fp, "\t--verbose INT\n\t\tVerbosity level; set to 8+ for debugging.  [DEFAULT: %d]\n", opt->info);
	fprintf(fp, "\t--help | -h\n\t\tThis help.\n");
//	fprintf(fp, "\t-k <kuint>\n\t\tNumber of haplotypes in the haplotype set (used with -i <hstr>).  [DEFAULT: %i]\n", opt->K);	/* KSD: get rid of this option */
	/* fprintf(fp, "\t--most <mint>\n\t\tReport top m-most abundant sequences and quit. [DEFAULT: %i]\n", opt->most_abundant); */
	fprintf(fp, "\n");
	for (size_t i = start; i < strlen(exe_name); ++i)
		fputc(toupper(exe_name[i]), fp);
	fprintf(fp, "(%d)\n", 1);
} /* fprint_usage */
