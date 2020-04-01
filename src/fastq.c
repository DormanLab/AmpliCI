#include <stdlib.h>

#include "constants.h"
#include "fastq.h"
#include "align.h"
#include "lmath.h"
#include "error.h"

/**
 * XY used penultimate 2 bits to encode 4 nucleotides.
 * IUPAC re-encodes using 4 bits to represent 17 characters, ignoring U as not
 *   distinct from T.
 *
 * char ASCII   binary xy_t iupac_t    binary
 *    A    65  1000001    0       1  00000001
 *    C    67  1000011    1       2  00000010
 *    G    71  1000111    3       4  00000100
 *    T    84  1010100    2       8  00001000
 *    U    85  1010101            8  00001000
 *    R    82  1010010            5  00000101
 *    Y    89  1011001           10  00001010
 *    S    83  1010011            6  00000110
 *    W    87  1010111            9  00001001
 *    K    75  1001011           12  00001100
 *    M    77  1001101            3  00000011
 *    B    66  1000010           14  00001110
 *    D    68  1000100           13  00001101
 *    H    72  1001000           11  00001011
 *    V    86  1010110            7  00000111
 *    N    78  1001110    -      15  00001111
 *    X    88  1011000            0  00000000
 */

/**
 * Number of nucleotides consistent with each IUPAC symbol. [NOT USED]
 */
const unsigned char popcnt[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};

/**
 * Convert xy_t to char for human consumption.
 */
unsigned char const xy_to_char[NUM_NUCLEOTIDES] = {'A', 'C', 'T', 'G'};

/**
 * Convert iupac_t to char for human consumption.
 */
unsigned char const iupac_to_char[NUM_IUPAC_SYMBOLS] = {
	'-', 'A', 'C', 'M', 'G', 'R', 'S',
	'V', 'T', 'W', 'Y', 'H', 'K',
	'D', 'B', 'N'
};

/**
 * Convert xy_t to iupac_t.
 */
iupac_t const xy_to_iupac[NUM_NUCLEOTIDES] = {
	IUPAC_A,
	IUPAC_C,
	IUPAC_T,
	IUPAC_G
};

/**
 * Convert xy_t to reverse complement xy_t.
 */
xy_t const xy_to_rc[NUM_NUCLEOTIDES] = {
	XY_T,
	XY_G,
	XY_A,
	XY_C
};

/**
 * Convert iupac_t to reverse complement iupac_t.
 */
iupac_t const iupac_to_rc[NUM_IUPAC_SYMBOLS] = {
	15, 8, 4, 2, 1, 10, 5, 9, 6, 3, 12, 1, 2, 4, 8, 0
};

/**
 * Convert iupac_t to xy_t: only meaningful for A, C, G, T.
 */
xy_t const iupac_to_xy[NUM_IUPAC_SYMBOLS] = {
	0, XY_A, XY_C, 0, XY_G, 0, 0, 0, XY_T, 0, 0, 0, 0, 0, 0, 0
};

/**
 * Convert iupac_t to standard nucleotide order A=0, C=1, G=2, T=3, which is
 * NOT xy_t.
 */
unsigned char const iupac_to_std[NUM_IUPAC_SYMBOLS] = {
	0, STD_A, STD_C, 0, STD_G, 0, 0, 0, STD_T, 0, 0, 0, 0, 0, 0, 0
};

/**
 * Standard order: A, C, G, T
 * XY order: A, C, T, G
 */
unsigned char const xy_to_std[NUM_NUCLEOTIDES] = {0, 1, 3, 2};
xy_t const std_to_xy[NUM_NUCLEOTIDES] = {0, 1, 3, 2};

/**
 * Convert char-encoded nucleotide - MIN_NUCLEOTIDE_ASCII to iupac_t
 * type.  Handles 'A' thru 'Y' and converts unknown to 0 aka 'X'.
 */
iupac_t const nuc_to_iupac[NUCLEOTIDE_ALPHABET_SIZE] = {
	1,14,2,13,0,0,4,11,0,0,12,0,3,15,0,0,0,5,6,8,8,7,9,0,10
	/*                      1                    2
        0  1 2  3 4 5 6  7 8 9  0 1 2  3 4 5 6 7 8 9 0 1 2 3  4 */
};

/**
 * Convert char-encoded nucleotide - MIN_NUCLEOTIDE_ASCII to xy_t
 * type.  Handles 'A' thr 'Y' but converts unknown to 0 aka 'A'.
 */
xy_t const nuc_to_xy[NUCLEOTIDE_ALPHABET_SIZE] = {
	0,0,1,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0
	/*                  1                   2
        0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 */
};

/**
 * Macro forward to the next newline or EOF.
 */
#define fforward(f, c, t) do {                                                 \
	(c) = fgetc((f));                                                      \
} while ((c) != EOF && (c) != (t));

/**
 * Macro forward to the next newline or EOF, and count characters consumed.
 */
#define fforward_cnt(f, c, t, a)                                               \
while (((c) = fgetc((f))) != EOF && (c) != (t)) {                              \
	(a)++;                                                                 \
}

/**
 * Validate human-readable nucleotide characters as IUPAC symbols or one
 * of standard nucleotides: A, C, G, T.
 */
extern int valid_iupac(unsigned char c);
extern int valid_nucleotide(unsigned char c);

/**
 * Convert quality score to probability.
 */
extern double error_prob(fastq_data *fqd, unsigned char c);

/**
 * Print all observed quality scores as probabilities.
 */
extern void fprint_error_probs(FILE *fp, fastq_data *fqd);

/**
 * Write one read of given length to file in R table format (space-separated
 * integers).
 */
extern void write_read_in_table(FILE *fp, unsigned char *read, unsigned int len);

/**
 * Return length of requested read.
 */
extern unsigned int read_length(fastq_data *fqd, unsigned int i);

int read_read(FILE *fp, fastq_data *fqd, unsigned int *len, unsigned char *nptr, unsigned char *qptr);
extern unsigned int number_nucleotide(iupac_t c);
extern const iupac_t *nucleotide_list(iupac_t c);

/**
 * Open fastq file and count number of reads.
 *
 * @param filename	fastq filename
 * @return		number of reads (0 for failure)
 */
unsigned int cnt_reads(char const * const filename)
{
	int fxn_debug = ABSOLUTE_SILENCE;	//DEBUG_I;	//
	FILE *fp = fopen(filename, "r");

	debug_msg(DEBUG_I, fxn_debug, "opening fastq file '%s'\n", filename);

	if (fp == NULL) {
		mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);
		return 0;
	}

	unsigned int cnt = fcnt_reads(fp);

	fclose(fp);

	return cnt;
} /* cnt_reads */

unsigned int fcnt_reads(FILE *fp)
{
	int fxn_debug = ABSOLUTE_SILENCE;	//DEBUG_I;	//
	int file_type = FASTQ_FILE;
	unsigned int nread = 0;
	char c;

	debug_msg(DEBUG_I, fxn_debug, "entering\n");

	c = fgetc(fp);
	if (c != '@' && c != '>') {
		mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "not fast[aq] file\n");
		return 0;
	}

	if (c == '>')
		file_type = FASTA_FILE;

	do {

		if (nread) {
			c = fgetc(fp);
			if (c == EOF)
				return nread;
			if (c != '@' && c != '>') {
				mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					"not fast[aq] file\n");
				return 0;
			}
		}

		/* fast forward through name */
		fforward(fp, c, '\n');

		if (c == EOF) {
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"premature end-of-file\n");
			return 0;
		}

		/* fast forward through read */
		fforward(fp, c, '\n');

		if (c == EOF) {
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"premature end-of-file\n");
			return 0;
		}

		if (file_type == FASTA_FILE)
			continue;

		/* fast forward through divider */
		fforward(fp, c, '\n');

		if (c == EOF) {
			mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"premature end-of-file\n");
			return 0;
		}

		/* fast forward through quality score */
		fforward(fp, c, '\n');

		nread++;
	} while (c != EOF);

	return nread;

} /* fcnt_reads */

/**
 * Index fast[aq] file.  Record byte index of the start of each read.
 *
 * @param fp	fast[aq] open file handle
 * @param fqd	fastq_data object pointer, with fastq_data.n_reads set
 * @return	error status
 */
int findex_reads(FILE *fp, fastq_data *fqd)
{
	int fxn_debug = ABSOLUTE_SILENCE;	//DEBUG_II;	//
	unsigned int i = 0;
	size_t nbytes = 0;
	char c;

	debug_msg(DEBUG_I, fxn_debug, "entering\n");

	c = fgetc(fp);
	if (c != '@' && c != '>')
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
			"not fast[aq] file\n");

	fqd->file_type = FASTQ_FILE;
	if (c == '>')
		fqd->file_type = FASTA_FILE;

	fqd->index = malloc(fqd->n_reads * sizeof *fqd->index);

	if (!fqd->index)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data.index");

	do {

		if (i) {
			c = fgetc(fp);
			if (c == EOF)
				return NO_ERROR;
			if (c != '@' && c != '>')
				return mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
					"not fast[aq] file\n");
		}

		fqd->index[i++] = nbytes;

		debug_msg(DEBUG_II, fxn_debug, "Read %u at byte %u.\n", i,
			nbytes);

		++nbytes;

		/* fast forward through name */
		fforward_cnt(fp, c, '\n', nbytes);

		if (c == EOF)
			return FASTQ_EOF;
		++nbytes;	/* for newline */

		/* fast forward through read */
		fforward_cnt(fp, c, '\n', nbytes);

		if (c == EOF)
			return FASTQ_EOF;
		++nbytes;

		if (fqd->file_type == FASTA_FILE)
			continue;

		/* fast forward through divider */
		fforward_cnt(fp, c, '\n', nbytes);

		if (c == EOF)
			return FASTQ_EOF;
		++nbytes;

		/* fast forward through quality score */
		fforward_cnt(fp, c, '\n', nbytes);
		++nbytes;

	} while (c != EOF);

	return NO_ERROR;

} /* findex_reads */

/**
 * Open, read and store fastq file in newly allocated memory.
 *
 * @param filename	fastq filename
 * @param in_fqd	unallocated fastq object
 * @param fqo		fastq options pointer
 * @return		error code
 */
int read_fastq(const char *filename, fastq_data **in_fqd, fastq_options *fqo)
{
	int fxn_debug = ABSOLUTE_SILENCE;	//SILENT;	//
	int err = NO_ERROR;
	FILE *fp = fopen(filename, "r");

	if (fxn_debug)
		fprintf(stderr, "%s:%d: opening fastq file '%s'\n", __func__,
			__LINE__, filename);

	if (fp == NULL)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename);

	err = fread_fastq(fp, in_fqd, fqo);

	fclose(fp);

	return err;
} /* read_fastq */

int fread_fastq(FILE *fp, fastq_data **in_fqd, fastq_options *fqo)
{
	int fxn_debug = ABSOLUTE_SILENCE;	//SILENT;	//DEBUG_III;	//
	unsigned char *rptr, *qptr;
	unsigned int n_reads, n_bytes, *uiptr;
	unsigned char c, elen = 1;
	fastq_data *fqd = *in_fqd;
	int err = NO_ERROR;

	debug_msg(DEBUG_I, fxn_debug, "entering\n");

	c = fgetc(fp);
	if (c != '@' && c != '>')
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "not fast[aq] file\n");

	rewind(fp);

	if (fqd == NULL) {
		*in_fqd = malloc(sizeof **in_fqd);
		if (*in_fqd == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"fastq_data reads");
		fqd = *in_fqd;
		fqd->file_type = FASTQ_FILE;
		fqd->n_lengths = NULL;
		fqd->reads = NULL;
		fqd->quals = NULL;
		fqd->reference_seq = NULL;
		fqd->index = NULL;
	}

	if (c == '>')
		fqd->file_type = FASTA_FILE;

	debug_msg(DEBUG_I, fxn_debug, "file type = %d\n", fqd->file_type);

	fqd->read_encoding = fqo->read_encoding == DEFAULT_ENCODING
		? XY_ENCODING : fqo->read_encoding;
	fqd->min_quality = MAX_ASCII_QUALITY_SCORE;
	fqd->max_quality = MIN_ASCII_QUALITY_SCORE;
	fqd->n_reads = 0;
	fqd->n_min_length = (unsigned int) -1;	/* maximum unsigned int */
	fqd->empty = 1;

	n_bytes = 0;
	n_reads = 0;
	do {
		err = read_read(fp, fqd, &fqd->n_max_length, NULL, NULL);

		debug_msg(DEBUG_III, fxn_debug, "err=%d (%s), length=%u\n",
			err, fastq_error_message(err), fqd->n_max_length);

		/* valid read, increment read count and cumulative length */
		if (err == NO_ERROR || (err == FASTQ_EOF && fqd->n_max_length > 0)) {

			n_bytes += fqd->n_max_length;
			if (fqd->n_max_length < fqd->n_min_length)
				fqd->n_min_length = fqd->n_max_length;
			fqd->n_reads++;
			++n_reads;
			debug_msg(DEBUG_III, fxn_debug, "read sequence %d of "
				"length %d\n", fqd->n_reads, fqd->n_max_length);

		/* tolerate ambiguous characters: drop these reads */
		} else if (err == FASTQ_AMBIGUOUS_READ_CHAR &&
			fqo->drop_invalid_reads) {

			mmessage(WARNING_MSG, err, "%s in read %u: dropping "
				"read\n", fastq_error_message(err), n_reads);
			++n_reads;

		/* tolerate ambiguous characters: change to iupac encoding */
		} else if (err == FASTQ_AMBIGUOUS_READ_CHAR
			&& fqo->read_encoding != XY_ENCODING) {

			if (fqd->read_encoding != IUPAC_ENCODING)
				mmessage(WARNING_MSG, err, "%s in read %u: "
					"using iupac encoding\n",
					fastq_error_message(err),
					fqd->n_reads + 1);
			fqd->read_encoding = IUPAC_ENCODING;

		/* cannot tolerate ambiguous characters: abort */
		} else if (err == FASTQ_AMBIGUOUS_READ_CHAR) {
			mmessage(ERROR_MSG, err, "%s in read %u: cannot encode "
				"using XY encoding\n", fastq_error_message(err),
				fqd->n_reads + 1);
			break;

		/* tolerate invalid characters: drop read */
		} else if (err == FASTQ_INVALID_READ_CHAR || err == FASTQ_INVALID_QUALITY_CHAR)

			mmessage(WARNING_MSG, err, "%s: read %u will be "
				"discarded\n", fastq_error_message(err),
				fqd->n_reads + 1);

		/* other return codes indicate EOF or irrecoverable error */
		else
			break;
	} while (err != FASTQ_EOF);

	if (err && err != FASTQ_EOF)
		return err;

	debug_msg(DEBUG_I, fxn_debug, "Found %u sequences.\n", n_bytes);

	fqd->reads = malloc(n_bytes * sizeof *fqd->reads);
	if (fqd->reads == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data reads");

	fqd->quals = malloc(n_bytes * sizeof *fqd->quals);
	if (fqd->quals == NULL) {
		free(fqd->reads);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data quals");
	}

	fqd->n_lengths = malloc(fqd->n_reads * sizeof *fqd->n_lengths);
	if (fqd->n_lengths == NULL) {
		free(fqd->reads);
		free(fqd->quals);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data.n_lengths");
	}

	rewind(fp);

	rptr = fqd->reads;
	qptr = fqd->quals;
	uiptr = fqd->n_lengths;
	n_reads = 0;
	do {
		err = read_read(fp, fqd, &n_bytes, rptr, qptr);

		/* increment pointers for next read */
		if (err == NO_ERROR || (err == FASTQ_EOF && n_bytes > 0)) {
			if (!n_reads) fqd->n_max_length = n_bytes;
			else if (n_bytes != fqd->n_min_length) elen = 0;
			if (n_bytes > fqd->n_max_length)
				fqd->n_max_length = n_bytes;
			n_reads++;
			if (fxn_debug)
				fprintf(stderr, "%s:%d: reread sequence %d of length %d: %.*s (%d)\n", __func__, __LINE__, n_reads, n_bytes, n_bytes, display_sequence(rptr, n_bytes, fqd->read_encoding), elen);
			*uiptr = n_bytes;
			rptr += n_bytes;
			qptr += n_bytes;
			uiptr++;

		/* tolerable errors */
		} else if (err != FASTQ_INVALID_READ_CHAR
			&& err != FASTQ_AMBIGUOUS_READ_CHAR
			&& err != FASTQ_INVALID_QUALITY_CHAR)
			break;
	} while (err != FASTQ_EOF);

	if (err && err != FASTQ_EOF)
		return err;

	if (n_reads != fqd->n_reads) {
		free(fqd->reads);
		free(fqd->quals);
		free(fqd->n_lengths);
		fprintf(stderr, "n_reads = %u; fastq::n_reads = %u\n", n_reads,
			fqd->n_reads);
		return mmessage(ERROR_MSG, FILE_FORMAT_ERROR, "fastq file");
	}

	if (fxn_debug >= DEBUG_III) {
		debug_msg(DEBUG_III, fxn_debug, "First 10 reads:\n");
		rptr = fqd->reads;
		for (unsigned int i = 0; i < 10; ++i) {
			debug_msg(DEBUG_III, fxn_debug, "Read %2u: %.*s\n", i + 1,
				fqd->n_lengths[i], display_sequence(rptr,
				fqd->n_lengths[i], fqd->read_encoding));
			rptr += fqd->n_lengths[i];
		}
	}

	if (fqd->file_type == FASTQ_FILE) {
		debug_msg(TALKATIVE, fxn_debug, "Minimum quality score: %c (%d)\n",
			fqd->min_quality, (int) fqd->min_quality);
		debug_msg(TALKATIVE, fxn_debug, "Maximum quality score: %c (%d)\n",
			fqd->max_quality, (int) fqd->max_quality);
	}
	debug_msg(TALKATIVE, fxn_debug, "Minimum read length: %u\n",
		fqd->n_min_length);
	debug_msg(TALKATIVE, fxn_debug, "Maximum read length: %u\n",
		fqd->n_max_length);

	qptr = fqd->quals;
	for (unsigned int i = 0; i < fqd->n_reads; ++i)
		for (unsigned int j = 0; j < fqd->n_lengths[i]; ++j) {
			*qptr = *qptr - fqd->min_quality;
			qptr++;
		}

	if (elen) {
		free(fqd->n_lengths);
		fqd->n_lengths = NULL;
	}
	fqd->empty = 0;

	return NO_ERROR;
} /* fread_fastq */

/**
 * Process a single read from a fastq file.
 *
 * @param fp	open fastq file handle
 * @param fqd	allocated fastq object
 * @param len	pointer to memory to store length of current read
 * @param nptr	pointer to memory to store read base characters
 * @param qptr	pointer to memory to store quality characters
 *
 * @return	error code
 */
int read_read(FILE *fp, fastq_data *fqd, unsigned int *len, unsigned char *nptr,
							 unsigned char *qptr)
{
	int fxn_debug = ABSOLUTE_SILENCE;	//SILENT;	//DEBUG_II;	//
	int err = NO_ERROR;
	unsigned int qlen = 0;
	int v_iupac, v_nuc;

	(*len) = 0;	/* signal incomplete read */

	/* check for @ at start of read record */
	char c = fgetc(fp);

	debug_msg(DEBUG_I, fxn_debug, "First character: '%c'\n", c);

	if (c != '@' && c != EOF && fqd->file_type == FASTQ_FILE)
		return FASTQ_FILE_FORMAT_ERROR;
	else if (c != '>' && c != EOF && fqd->file_type == FASTA_FILE)
		return FASTQ_FILE_FORMAT_ERROR;
	else if (c == EOF)
		return FASTQ_EOF;

	/* fast forward through name */
	fforward(fp, c, '\n');

	if (c == EOF)
		return FASTQ_EOF;

	/* read or skip the nucleotide sequence */
	debug_msg(DEBUG_II, fxn_debug, "");
	while ((c = fgetc(fp)) != '\n' && c != EOF) {

		debug_msg_cont(DEBUG_II, fxn_debug, "%c", c);

		v_iupac = err == FASTQ_INVALID_READ_CHAR ? 0 : valid_iupac(c);
		v_nuc = v_iupac
			? err == FASTQ_AMBIGUOUS_READ_CHAR ? 0 : valid_nucleotide(c)
			: 0;

		/* record valid read */
		if (nptr && v_iupac && fqd->read_encoding == IUPAC_ENCODING) {

			*nptr = nuc_to_iupac[c - 'A'];
			nptr++;

		} else if (nptr && v_nuc) {

			*nptr = (c >> 1) & 3;	/* unique encoding of A, C, G, T, but no others! */
			nptr++;

		/* invalid read: invalid character */
		} else if (!v_iupac) {

			err = FASTQ_INVALID_READ_CHAR;
			break;

		/* invalid read: cannot encode ambiguous nucleotide */
		} else if (!v_nuc && fqd->read_encoding == XY_ENCODING) {

			err = FASTQ_AMBIGUOUS_READ_CHAR;
			break;

		/* ambiguous nucleotide: keep parsing */
		} else if (!v_nuc)

			err = FASTQ_AMBIGUOUS_READ_CHAR;

		(*len)++;
	}

	debug_msg_cont(DEBUG_II, fxn_debug, " (length = %u)\n", *len);

	/* invalid read: skip rest */
	if (err == FASTQ_INVALID_READ_CHAR) {

		mmessage(WARNING_MSG, FASTQ_INVALID_READ_CHAR, "%s '%c'\n",
			fastq_error_message(FASTQ_INVALID_READ_CHAR), c);
		fforward(fp, c, '\n');

	/* ambiguous nucleotide: issue warning */
	} else if (err == FASTQ_AMBIGUOUS_READ_CHAR &&
		fqd->read_encoding == XY_ENCODING) {

		mmessage(WARNING_MSG, FASTQ_AMBIGUOUS_READ_CHAR, "%s",
			fastq_error_message(FASTQ_AMBIGUOUS_READ_CHAR));
		fprintf(stderr, " '%c'", c);
		fforward(fp, c, '\n');
		fprintf(stderr, "\n");
	}

	if (c == EOF) {
		if (fqd->file_type == FASTQ_FILE)
			*len = 0;	/* error: quality scores missing */
		return FASTQ_EOF;
	}

	if (fqd->file_type == FASTA_FILE)
		return NO_ERROR;

	/* fast forward through divider */
	fforward(fp, c, '\n');

	if (c == EOF) {
		*len = 0;
		return FASTQ_EOF;
	}

	/* read or skip the quality score sequence */
	if (qptr)
		while ((c = fgetc(fp)) != '\n' && c != EOF) {
			*qptr = (unsigned char) c;
			qptr++;
			qlen++;
			if (c < fqd->min_quality) fqd->min_quality = (unsigned char) c;
			if (c > fqd->max_quality) fqd->max_quality = (unsigned char) c;
		}
	else
		fforward_cnt(fp, c, '\n', qlen);

	if (c == EOF) {
		/* incomplete quality score string */
		if (qlen < *len)
			*len = 0;
		return FASTQ_EOF;
	}

	return err;
} /* read_read */

int allocate_empty_fastq(fastq_data **in_fqd, fastq_options *fqo,
			unsigned int nreads, unsigned int read_length)
{
	int fxn_debug = ABSOLUTE_SILENCE;	//SILENT;	//DEBUG_III;	//
	fastq_data *fqd = *in_fqd;

	debug_msg(DEBUG_I, fxn_debug, "entering\n");

	if (fqd == NULL) {
		*in_fqd = malloc(sizeof **in_fqd);
		if (*in_fqd == NULL)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
				"fastq_data reads");
		fqd = *in_fqd;
		fqd->file_type = FASTQ_FILE;
		fqd->n_lengths = NULL;
		fqd->reads = NULL;
		fqd->quals = NULL;
		fqd->reference_seq = NULL;
		fqd->index = NULL;
	}

	fqd->empty = 1;
	fqd->read_encoding = fqo->read_encoding == DEFAULT_ENCODING
		? XY_ENCODING : fqo->read_encoding;
	fqd->min_quality = MIN_ASCII_QUALITY_SCORE;
	fqd->max_quality = MAX_ASCII_QUALITY_SCORE;
	fqd->n_reads = nreads;
	fqd->n_min_length = fqd->n_max_length = read_length;

	fqd->reads = malloc(nreads * read_length * sizeof *fqd->reads);
	if (fqd->reads == NULL)
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data reads");

	fqd->quals = malloc(nreads * read_length * sizeof *fqd->quals);
	if (fqd->quals == NULL) {
		free(fqd->reads);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data quals");
	}

	fqd->n_lengths = NULL;

	return NO_ERROR;
}/* allocate_empty_fastq */

unsigned char const * display_sequence(unsigned char const * const in_str, unsigned int len, int encoding) {
	unsigned int j;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION, "string to store read");
		return NULL;
	}

	for (j = 0; j < len; ++j)
		str[j] = encoding == IUPAC_ENCODING
			? iupac_to_char[(int) in_str[j]]
			: xy_to_char[(int) in_str[j]];
	str[j] = '\0';

	return str;
} /* display_sequence */

unsigned char const * display_quals(unsigned char const * const in_str, unsigned int len, unsigned char min) {
	unsigned int j;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"string to store qualities");
		return NULL;
	}

	for (j = 0; j < len; ++j)
		str[j] = in_str[j] + min;
	str[j] = '\0';

	return str;
} /* display_quals */

unsigned char const * display_reverse_complement(unsigned char const * const in_str, unsigned int len, int encoding) {
	unsigned int j, l;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"string to store qualities");
		return NULL;
	}

	for (j = len - 1, l = 0; l < len; --j, ++l)
		str[l] = encoding == IUPAC_ENCODING
			? iupac_to_char[(int) iupac_to_rc[(int) in_str[j]]]
			: xy_to_char[(int) xy_to_rc[(int) in_str[j]]];
	str[l] = '\0';

	return str;
} /* display_reverse_complement */

unsigned char const * display_reverse_quals(unsigned char const * const in_str, unsigned int len, unsigned char min) {
	unsigned int j, l;
	unsigned char *str = malloc((len + 1) * sizeof *str);

	if (str == NULL) {
		mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"string to store qualities");
		return NULL;
	}

	for (j = len - 1, l = 0; l < len; --j, ++l)
		str[l] = in_str[j] + min;
	str[l] = '\0';

	return str;
} /* display_reverse_quals */

int write_fastq(fastq_data *fqd, fastq_options *fqo)
{
	unsigned int i;
	unsigned char *reads = fqd->reads;
	unsigned char *quals = fqd->quals;
	FILE *fp = fopen(fqo->outfile, fqo->append ? "a" : "w");

	if (!fp)
		return(mmessage(ERROR_MSG, FILE_OPEN_ERROR, fqo->outfile));

	for (i = 0; i < fqd->n_reads; ++i) {
		fprintf(fp, "%c%u\n", fqo->fasta ? '>' : '@', i);
		fprintf(fp, "%s\n", !fqo->reverse_complement
			? display_sequence(reads, read_length(fqd, i),
				fqd->read_encoding)
			: display_reverse_complement(reads, read_length(fqd, i),
				fqd->read_encoding)
			);
		if (!fqo->fasta) {
			fprintf(fp, "+\n");
			fprintf(fp, "%s\n", !fqo->reverse_complement
				? display_quals(quals, read_length(fqd, i),
					fqd->min_quality)
				: display_reverse_quals(quals,
					read_length(fqd, i), fqd->min_quality)
				);
		}
		reads += read_length(fqd, i);
		if (!fqo->fasta)
			quals += read_length(fqd, i);
	}
	fclose(fp);

	return NO_ERROR;
} /* write_fastq */

int write_fastq_marked(fastq_data *fqd, fastq_options *fqo, unsigned int *id,
	unsigned int selected_id)
{
	unsigned int i;
	unsigned char *reads = fqd->reads;
	unsigned char *quals = fqd->quals;
	FILE *fp = fopen(fqo->outfile, fqo->append ? "a" : "w");

	if (!fp)
		return(mmessage(ERROR_MSG, FILE_OPEN_ERROR, fqo->outfile));

	for (i = 0; i < fqd->n_reads; ++i) {
		if (id[i] == selected_id) {
			fprintf(fp, "%c%u\n", fqo->fasta ? '>' : '@', i);
			fprintf(fp, "%s\n", !fqo->reverse_complement
				? display_sequence(reads, read_length(fqd, i),
					fqd->read_encoding)
				: display_reverse_complement(reads,
					read_length(fqd, i), fqd->read_encoding)
				);
			if (!fqo->fasta) {
				fprintf(fp, "+\n");
				fprintf(fp, "%s\n", !fqo->reverse_complement
					? display_quals(quals,
						read_length(fqd, i),
						fqd->min_quality)
					: display_reverse_quals(quals,
						read_length(fqd, i),
						fqd->min_quality)
					);
			}
		}
		reads += read_length(fqd, i);
		if (!fqo->fasta)
			quals += read_length(fqd, i);
	}

	return NO_ERROR;
} /* write_fastq_marked */

/**
 * Write fastq data as R-style data table.
 *
 * @param fqd		fastq_data object pointer
 * @param filename	name of output file
 * @return		error status
 */
int write_table(fastq_data *fqd, char const *filename)
{
	unsigned int i, len;
	unsigned char *reads = fqd->reads;
	FILE *fp = fopen(filename, "w");

	if (!fp)
		return(mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename));

	for (i = 0; i < fqd->n_reads; ++i) {
		len = read_length(fqd, i);
		write_read_in_table(fp, reads, len);
		reads += len;
	}
	fclose(fp);

	return NO_ERROR;
} /* write_table */

/**
 * Free fastq object.
 *
 * @param fqd	fastq object pointer
 */
void free_fastq(fastq_data *fqd) {
	if (fqd) {
		if (fqd->reads) free(fqd->reads);
		if (fqd->quals) free(fqd->quals);
		if (fqd->n_lengths) free(fqd->n_lengths);
		free(fqd);
	}
} /* free_fastq */

/**
 * Return human-friendly error message for fastq error code.
 *
 * @param err_no	error number
 * @return		error message string
 */
const char *fastq_error_message(int err_no) {
	if (err_no == FASTQ_INVALID_READ_CHAR)
		return "illegal nucleotide";
	else if (err_no == FASTQ_AMBIGUOUS_READ_CHAR)
		return "ambiguous nucleotide";
	else if (err_no == FASTQ_INVALID_QUALITY_CHAR)
		return "illegal quality score";
	else if (err_no == FASTQ_INCOMPLETE_READ)
		return "incomplete read";
	else if (err_no == FASTQ_FILE_FORMAT_ERROR)
		return "invalid file format";
	else if (err_no == FASTQ_EOF)
		return "end-of-file";
	else
		return "No error";
} /* fastq_error_message */

/**
 * Check if two reads are equal.  This function assumes there is biological
 * homology at read positions.
 *
 * This may be most useful outside fastq_data where the homology
 * relationship among reads is known.
 *
 * @param fqd	fastq object pointer
 * @param i	index of first read
 * @param j	index of second read
 * @return	<0, 0, >0
 */
int read_compare(fastq_data *fqd, unsigned int i, unsigned int j) {
	/* same index */
	if (i == j)
		return 0;

	unsigned int len = fqd->n_max_length;
	size_t i_idx = len * i, j_idx = len * j;

	if (fqd->n_lengths) {

		/* assume reads of different lengths not equal */
		if (fqd->n_lengths[i] < fqd->n_lengths[j])
			return -1;
		else if (fqd->n_lengths[i] > fqd->n_lengths[j])
			return 1;

		/* they share a common length */
		len = fqd->n_lengths[i];

		/* find index of ith and jth read */
		i_idx = j_idx = 0;
		for (unsigned int l = 0; l < MIN(i, j); ++l) {
			i_idx += fqd->n_lengths[l];
			j_idx += fqd->n_lengths[l];
		}
		if (i < j)
			for (unsigned int l = i; l < j; ++l)
				j_idx += fqd->n_lengths[l];
		else
			for (unsigned int l = j; l < i; ++l)
				i_idx += fqd->n_lengths[l];
	}

	for (unsigned int l = 0; l < len; ++l) {
		if (fqd->reads[i_idx + l] < fqd->reads[j_idx + l])
			return -1;
		else if (fqd->reads[i_idx + l] > fqd->reads[j_idx + l])
			return 1;
	}
	return 0;
} /* read_compare */

int pw_align_reads(fastq_data *fqd, char const * const rfile) {
	int err = NO_ERROR;
	FILE *fp = fopen(rfile, "r");	/* fasta format */

	if (!fp)
		return mmessage(ERROR_MSG, FILE_OPEN_ERROR, rfile);

	char c = fgetc(fp);
	if (c != '>') {
		fclose(fp);
		return FILE_FORMAT_ERROR;
	}

	fforward(fp, c, '\n');

	unsigned int len = 0;
	while ((c = fgetc(fp)) != EOF) {
		if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
			len++;
		else if (c != '\n' && c != ' ' && c != '\t') {
			fclose(fp);
			return mmessage(ERROR_MSG, INVALID_USER_INPUT,
				"Reference sequence contains ambiguous "
				"character: %c\n", c);
		}
	}

	fqd->reference_seq = malloc(len * sizeof *fqd->reference_seq);

	if (!fqd->reference_seq) {
		fclose(fp);
		return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
			"fastq_data.reference_seq");
	}

	rewind(fp);
	fforward(fp, c, '\n');

	for (unsigned int i = 0; i < len; ++i) {
		c = fgetc(fp);
		if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
			fqd->reference_seq[i] = (c >> 1) & 3;
		else
			--i;
	}

	fclose(fp);

	mmessage(INFO_MSG, NO_ERROR, "Reference sequence: ");
	for (unsigned int i = 0; i < len; ++i)
		fprintf(stderr, "%c", xy_to_char[(int) fqd->reference_seq[i]]);
	fprintf(stderr, "\n");

	unsigned char *rptr = fqd->reads;
	int score[NUM_NUCLEOTIDES][NUM_NUCLEOTIDES] = {{2, -3, -3, -2},
		{-3, 2, -2, -3}, {-3, -2, 2, -3}, {-2, -3, -3, 2}};
	double const perr[] = {0.999401, 0.992814, 0.993413, 0.997006, 0.996407, 0.994910, 0.994311, 0.742627, 0.993713, 0.995808, 0.992216, 0.997305, 0.997006, 0.997904, 0.968593, 0.993114, 0.994611, 0.995808, 0.997305, 0.998204, 0.975150, 0.998503, 0.998204, 0.994611, 0.984731, 0.967365, 0.996707, 0.997904, 0.915868, 0.984431, 0.987126, 0.997605, 0.979042, 0.993114, 0.994311, 0.989820, 0.985629, 0.993114, 0.985329, 0.977246, 0.995808, 0.997305, 0.986527, 0.996108, 0.997006, 0.988623, 0.989532, 0.970659, 0.944346, 0.998204, 0.989820, 0.996707, 0.996707, 0.960778, 0.982934, 0.986527, 0.998204, 0.986527, 0.998503, 0.981737, 0.991916, 0.991018, 0.995210, 0.985329, 0.991916, 0.978789, 0.972156, 0.970060, 0.994963, 0.990719, 0.992814, 0.994611, 0.991916, 0.985917, 0.976700, 0.990419, 0.996707, 0.991018, 0.982635, 0.985329, 0.997305, 0.986228, 0.978144, 0.997006, 0.994311, 0.994012, 0.996108, 0.985928, 0.971856, 0.962874, 0.980838, 0.986228, 0.988323, 0.994311, 0.915194, 0.971512, 0.994311, 0.968862, 0.977545, 0.981437, 0.985329, 0.997006, 0.995210, 0.994012, 0.993450, 0.987183, 0.993413, 0.991916, 0.996379, 0.998503, 0.985629, 0.993513, 0.998503, 0.987725, 0.975449, 0.981138, 0.979641, 0.961770, 0.996108, 0.994910, 0.981437, 0.981437, 0.985329, 0.994311, 0.965015, 0.990075, 0.953293, 0.979341, 0.990120, 0.991617, 0.993114, 0.995210, 0.985928, 0.996707, 0.997904, 0.985329, 0.998204, 0.979341, 0.991916, 0.979641, 0.981437, 0.984132, 0.991617, 0.997006, 0.959880, 0.991916, 0.996707, 0.989521, 0.997006, 0.961386, 0.993413, 0.972156, 0.995509, 0.973653, 0.990120, 0.996707, 0.976647, 0.988323, 0.997605, 0.988922, 0.985629, 0.965269, 0.998503, 0.977545, 0.971557, 0.973952, 0.981737, 0.992814, 0.986527, 0.981737, 0.995509, 0.994311, 0.981737, 0.997006, 0.973054, 0.965569, 0.994311, 0.978144, 0.972455, 0.994311, 0.990719, 0.988623, 0.997305, 0.964072, 0.988623, 0.986527, 0.991018, 0.995210, 0.996407, 0.981437, 0.971788, 0.959581, 0.982335, 0.992515, 0.993713, 0.991617, 0.993114, 0.988623, 0.986826, 0.994611, 0.987126, 0.967504, 0.997305, 0.936483, 0.994910, 0.996707, 0.982335, 0.996407, 0.997305, 0.947006, 0.985940, 0.994910, 0.963473, 0.950299, 0.953593, 0.994311, 0.972156, 0.995210, 0.989222, 0.920475, 0.997305, 0.941018, 0.988024, 0.971557, 0.960479, 0.989222, 0.994910, 0.988323, 0.977246, 0.996707, 0.967365, 0.912167, 0.948165, 0.995509, 0.979940, 0.985050, 0.944311, 0.973353, 0.947305, 0.990719, 0.987126, 0.970958, 0.975150, 0.997006, 0.992814, 0.932335, 0.948802, 0.933832, 0.937504, 0.991317, 0.982335, 0.991617, 0.967365, 0.956287, 0.996108, 0.960180, 0.985329, 0.994311, 0.971257, 0.994611, 0.994311, 0.937126, 0.994311, 0.995210, 0.997904, 0.997305, 0.904790, 0.923653, 0.926048, 0.901982, 0.998503, 0.986826, 0.964970, 0.997605};

	for (unsigned int i = 0; i < fqd->n_reads; ++i) {
		size_t alen;
		unsigned char **aln = nwalign(fqd->reference_seq, rptr, len,
			read_length(fqd, i), score, -1, -1, 1, perr, &err,
			&alen);
		fprintf(stderr, "Read %u alignment length %lu\n", i, alen);
		size_t ngap1 = 0, ngap2 = 0;
		for (size_t j = 0; j < alen; ++j) {
			if (aln[0][j] == '-') ngap1++;
			fprintf(stderr, "%c", aln[0][j] == '-'
				? '-' : xy_to_char[(int) aln[0][j]]);
		}
		fprintf(stderr, "\n");
		for (size_t j = 0; j < alen; ++j) {
			if (aln[1][j] == '-') ngap2++;
			fprintf(stderr, "%c", aln[1][j] == '-'
				? '-' : xy_to_char[(int) aln[1][j]]);
		}
		fprintf(stderr, "\nGaps: %lu %lu\n", ngap1, ngap2);
		if (ngap1 != ngap2) exit(0);
		rptr += read_length(fqd, i);
	}

	return NO_ERROR;
} /* pw_align_reads */

double read_distance(fastq_data *fqd, unsigned int i, unsigned int j)
{
	if (i == j)
		return 0;

	unsigned char *qptr1 = NULL, *qptr2 = NULL;
	unsigned char *rptr1 = NULL, *rptr2 = NULL;
	unsigned char *rptr = fqd->reads;
	unsigned char *qptr = fqd->quals;
	unsigned int len;

	for (unsigned int n = 0; n < fqd->n_reads; ++n) {
		if (i == n) {
			rptr1 = rptr;
			qptr1 = qptr;
		}
		if (j == n) {
			rptr2 = rptr;
			qptr2 = qptr;
		}
		if ((i < j && j == n) || (j < i && i == n))
			break;

		len = read_length(fqd, n);
		rptr += len;
		qptr += len;
	}

	len = fqd->n_lengths ? MIN(fqd->n_lengths[i], fqd->n_lengths[j])
		: fqd->n_max_length;
	return read_distance_ptr(fqd, len, rptr1, rptr2, qptr1, qptr2);
} /* read_distance */

double read_distance_ptr(fastq_data *fqd, unsigned int len, unsigned char *rptr1,
	unsigned char *rptr2, unsigned char *qptr1, unsigned char *qptr2)
{
	if (qptr1 == qptr2) {
fprintf(stderr, "here!\n");
		return 0;
	}
	double dis = 0;
	for (unsigned int n = 0; n < len; ++n) {
/*
		fprintf(stderr, "%c %c %f %f", fqd->read_encoding == IUPAC_ENCODING
                        ? iupac_to_char[(int) *rptr1]
                        : xy_to_char[(int) *rptr1],
			fqd->read_encoding == IUPAC_ENCODING
                        ? iupac_to_char[(int) *rptr2]
                        : xy_to_char[(int) *rptr2],
			error_prob(fqd, *qptr1), error_prob(fqd, *qptr2));
*/
		double ldis = 0;
		if (*rptr1 == *rptr2) {
			ldis += error_prob(fqd, *qptr1) * (1 - error_prob(fqd, *qptr2));
			ldis += error_prob(fqd, *qptr2) * (1 - error_prob(fqd, *qptr1));
			ldis += 2*error_prob(fqd, *qptr2) * error_prob(fqd, *qptr1)/3;
		} else {
			ldis += error_prob(fqd, *qptr1)/3 * (1 - error_prob(fqd, *qptr2));
			ldis += error_prob(fqd, *qptr2)/3 * (1 - error_prob(fqd, *qptr1));
			ldis += 2*error_prob(fqd, *qptr2) * error_prob(fqd, *qptr1) / 9;
		}
		dis += ldis;
//fprintf(stderr, ": %f (%f)\n", ldis, dis);
		qptr1++;
		qptr2++;
		rptr1++;
		rptr2++;
	}
	return dis;
} /* read_distance_ptr */

int make_fastq_options(fastq_options **opt)
{
	fastq_options *op;
	*opt = malloc(sizeof **opt);
	if (!*opt)
		return(mmessage(ERROR_MSG, MEMORY_ALLOCATION, "fastq_data"));
	op = *opt;

	op->outfile = NULL;
	op->read_encoding = DEFAULT_ENCODING;
	op->reverse_complement = 0;
	op->fasta = 0;
	op->drop_invalid_reads = 0;

	return NO_ERROR;
} /* make_fastq_options */

void free_fastq_options(fastq_options *opt) {
	free(opt);
	UNUSED(opt);
} /* free_fastq_options */
