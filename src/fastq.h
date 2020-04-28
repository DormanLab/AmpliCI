/**
 * @file fastq.c
 * @author Karin S. Dorman
 *
 * FASTA/FASTQ library.
 *
 * DATA ENCODING:
 * Because NGS data are big data, it can be important to consider how they
 * are internally represented.  The four nucleotides, A, C, G, and T, can be
 * stored in 2 bits, the IUPAC nucleotide codes in 4 (if we treat U and T as
 * equivalent).  FASTA/FASTQ files represent the nucleotides (and quality
 * scores) as 8 bit characters even though all 8 bits are not necessary, so
 * some of the characters do not correspond to valid codes.  Because the
 * printable ASCII characters are in the mid-range of the 8-bit integers, it is
 * also possible to subtract 'A' (the integer interpretation of character 'A')
 * from the FASTA/FASTQ characters to encode the nucleotides as integers in the
 * range 'A' - 'A' (0) to 'Y' - 'A' (24).
 *
 * When reading FASTA/FASTQ files, we convert to the 4-bit encoding (iupac_t),
 * or the 2-bit encoding (xy_t).  Only the former can encode ambiguous
 * nucleotides such as R, Y, or N.
 *
 * Beware: There must be internal consistency in the orderings encoding.  For
 * example, which bit of the 4-bit encoding represents A?  I use the order (high
 * bit) TGCA (low bit).
 * The xy_t encoding orders the nucleotides A=0, C=1, T=2, G=3, because it makes
 * for convenient conversion between the FASTA/FASTQ encoding and the 2-bit
 * encoding, but many people assume the alphabetic ordering A, C, G, T, or the
 * biochemical ordering G, A, C, T we used to use on sequencing gels (Why did we
 * do that exactly?).  The encodings are shown in gory detail in fastq.c.
 *
 * MORE EFFICIENT ENCODING FOR FUTURE:
 * Despite all this fancy stuff about encoding, the current library STILL STORES
 * the nucleotides and quality scores in 8 bits, thus wasting at least 4 bits of
 * unused memory.  The reason is that it is substantially more cumbersome to
 * work with 2 or 4-bit units rather than 8 bit units, and I decided not to
 * undertake that challenge yet.  IDEA: store nucleotides (2 bits) in the low
 * bits, quality scores (6 bits) in the high bits.
 *
 * VARIABLE READ LENGTHS:
 * This implementation handles variable length reads (fastq_data:n_lengths), but
 * does not provide pointers to the start of each read.  There is a slight
 * saving of memory because read lengths are assumed to be rather short
 * (unsigned int), while pointers are large.  Unless you store pointers to the
 * start of each read in fastq_data:reads or qualities in fastq:quals, you will
 * have to use pointer arithmetic to find individual read data.  If reads are
 * all the same length, then fastq_data:n_lengths is NULL and
 * fastq_data:n_max_length is the shared length.  This inconsistent interface
 * to the read length is a nuisance, so see read_length() for access.
 *
 * DATA TYPES AND LIMITS:
 * This implementation assumes the number and lengths of reads can be stored in
 * unsigned int, and there is no verification of this truth.  Be cautious.
 */

#ifndef __H_FASTQ_DATA__
#define __H_FASTQ_DATA__

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "error.h"

/**
 * Types of nucleotide encodings.
 */
enum {
	DEFAULT_ENCODING,	/*!< default: xy_t unless need iupac_t */
	IUPAC_ENCODING,		/*!< iupac_t */
	XY_ENCODING,		/*!< xy_t */
	NUC_ENCODING,		/*!< nuc_t: should not use */
	STD_ENCODING        /*!< std_t */
};

/**
 * IUPAC symbols encoded in 4 lower bits.
 */
typedef unsigned char iupac_t;

/**
 * A, C, G, T as encoded by Xin Yin in 2 lower bits.
 */
typedef unsigned char xy_t;

/**
 * IUPAC symbols and 'X' encoded (with some gaps) as chars 0 to 24.
 */
typedef unsigned char nuc_t;

/**
 * Types of files.
 */
enum {
	FASTQ_FILE,
	FASTA_FILE,
	NUM_FILE_TYPES
};

/**
 * Generated errors.
 */
enum {
	FASTQ_INVALID_READ_CHAR		/*!< invalid character in read */
		= NUM_ERRORS + 1,
	FASTQ_AMBIGUOUS_READ_CHAR,	/*!< ambiguous character in read */
	FASTQ_INVALID_QUALITY_CHAR,	/*!< invalid quality score (NOT USED) */
	FASTQ_INCOMPLETE_READ,		/*!< file ends before end of read */
	FASTQ_FILE_FORMAT_ERROR,	/*!< other file format error */
	FASTQ_EOF,			/*!< end of fast[qa] file */
	FASTQ_NUM_ERRORS		/*!< number of possible errors */
};

/**
 * Number of nucleotides: A, C, G, T.
 */
#define NUM_NUCLEOTIDES	4
#define NUM_IUPAC_SYMBOLS	16

/**
 * A, C, G, T as xy_t, iupac_t, and nuc_t.
 */
enum {
	XY_A = 0,
	XY_C = 1,
	XY_G = 3,
	XY_T = 2
};
enum {
	IUPAC_A = 1,
	IUPAC_C = 2,
	IUPAC_G = 4,
	IUPAC_T = 8
};
enum {
	STD_A,	/* 0 */
	STD_C,	/* 1 */
	STD_G,	/* 2 */
	STD_T	/* 3 */
};

/**
 * Minimum allowed nucleotide sequence character in ASCII is char 'A'.
 */
#define MIN_NUCLEOTIDE_ASCII	'A'

/**
 * Maximum number of letters in ASCII nucleotide alphabet: 'A' to 'Y'
 */
#define NUCLEOTIDE_ALPHABET_SIZE	('Z' - 'A')

/**
 * Convert xy_t, iupac_t to display char.
 */
extern unsigned char const xy_to_char[NUM_NUCLEOTIDES];
extern unsigned char const iupac_to_char[NUM_IUPAC_SYMBOLS];

/**
 * Convert among iupac_t, xy_t, std, and nuc_t.
 */
extern unsigned char const iupac_to_xy[NUM_IUPAC_SYMBOLS];
extern unsigned char const iupac_to_std[NUM_IUPAC_SYMBOLS];

extern unsigned char const xy_to_std[NUM_NUCLEOTIDES];
extern unsigned char const xy_to_iupac[NUM_NUCLEOTIDES];
extern unsigned char const std_to_xy[NUM_NUCLEOTIDES];
extern iupac_t const nuc_to_iupac[NUCLEOTIDE_ALPHABET_SIZE];
extern xy_t const nuc_to_xt[NUCLEOTIDE_ALPHABET_SIZE];

/**
 * Reverse complement the codes.
 */
extern xy_t const xy_to_rc[NUM_NUCLEOTIDES];
extern iupac_t const iupac_to_rc[NUM_IUPAC_SYMBOLS];

#define NUM_IUPAC_SYMBOLS	16
extern unsigned char iupac_symbols[NUM_IUPAC_SYMBOLS];
extern const unsigned char popcnt[NUM_IUPAC_SYMBOLS];

/**
 * Minimum and maximum quality score.
 * [TODO] we assume a particular type of Illumina data
 */
#define MAX_QUALITY_SCORE       40      /*<! 73 in ASCII (modern: 93, 126 ASCII) */
#define MIN_QUALITY_SCORE       0       /*<! 33 in ASCII */
#define MAX_ASCII_QUALITY_SCORE	73	/*<! ~ */
#define MIN_ASCII_QUALITY_SCORE	33	/*<! ! */

typedef struct _fastq_data fastq_data;
typedef struct _fasta_data fasta_data;
typedef struct _fastq_options fastq_options;

/**
 * FASTQ data structure.
 */
struct _fastq_data {
	int empty;			/*!< whether there is data loaded */
	int file_type;			/*!< type of file */
	int read_encoding;		/*!< read encoding used */
	unsigned int n_reads;		/*!< number of reads */
	unsigned int n_max_length;	/*!< length of reads if identical */
	unsigned int n_min_length;	/*!< length of shortest read */
	unsigned char max_quality;	/*!< maximum observed quality score */
	unsigned char min_quality;	/*!< minimum observed quality score */
	size_t *index;			/*!< byte pointer of reads in file */
	unsigned int *n_lengths;	/*!< length of reads if not identical */
	unsigned char *reads;		/*!< sequence reads */
	unsigned char *quals;		/*!< quality score strings */
	unsigned char *reference_seq;	/*!< optional reference sequence */
};

/**
 * FASTQ options.
 */
struct _fastq_options {
	char drop_invalid_reads;/*!< drop invalid reads */
	char const *outfile;	/*!< outfile */
	int append;		/*!< append to outfile */
	int read_encoding;	/*!< which read encoding to request */
	int reverse_complement;	/*!< reverse complement the data */
	int fasta;		/*!< convert to fasta format */
};

/* read fasta/fastq files */
int allocate_empty_fastq(fastq_data **in_fqd, fastq_options *fqo, unsigned int nreads, unsigned int read_length);
int fread_fastq(FILE *fp, fastq_data **fqd, fastq_options *fqo);
int read_fastq(const char *filename, fastq_data **fqd, fastq_options *fqo);
int read_read(FILE *fp, fastq_data *fqd, unsigned int *len, unsigned char *nptr, unsigned char *qptr);
unsigned int cnt_reads(char const * const filename);
unsigned int fcnt_reads(FILE *fp);
int findex_reads(FILE *fp, fastq_data *fqd);
int make_fastq_options(fastq_options **opt);

/* do stuff with reads */
int read_compare(fastq_data *fqd, unsigned int i, unsigned int j);
int pw_align_reads(fastq_data *fqd, char const * const rfile);
double read_distance(fastq_data *fqd, unsigned int i, unsigned int j);
double read_distance_ptr(fastq_data *fqd, unsigned int len, unsigned char *rptr1, unsigned char *rptr2, unsigned char *qptr1, unsigned char *qptr2);

/* output */
unsigned char const * display_sequence(unsigned char const * const in_str, unsigned int len, int encoding);
unsigned char const * display_reverse_complement(unsigned char const * const in_str, unsigned int len, int encoding);
unsigned char const * display_quals(unsigned char const * const in_str, unsigned int len, unsigned char min);
unsigned char const * display_reverse_quals(unsigned char const * const in_str, unsigned int len, unsigned char min);
int write_fastq(fastq_data *fqd, fastq_options *fqo);
int write_fastq_marked(fastq_data *fqd, fastq_options *fqo, unsigned int *id, unsigned int selected_id);
int write_table(fastq_data *fqd, char const *filename);

const char *fastq_error_message(int err_no);

void free_fastq(fastq_data *fqd);
void free_fastq_options(fastq_options *opt);


/**
 * Check whether character is valid IUPAC symbol.  Incoming character is
 * human-readable letter and output is 0/1 to indicate if it is one of the
 * allowed IUPAC symbols.
 *
 * @param c	human-readable character
 * @return	0|1
 */
inline int valid_iupac(unsigned char c) {
	for (size_t i = 0; i < NUM_IUPAC_SYMBOLS; ++i)
		if (iupac_to_char[i] == c)
			return 1;
	if (c == 'X')
		return 1;
	return 0;
} /* valid_iupac */

/**
 * Check whether character is valid nucleotide.  Incoming character is
 * human-readable letter and output is 0/1 to indicate if it is one of
 * A, C, G, or T.
 *
 * @param c	human-readable character
 * @return	1|0
 */
inline int valid_nucleotide(unsigned char c) {
	return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U';
} /* valid_nucleotide */

/**
 * Convert quality code to error probability.
 *
 * @param q	quality code encoded as ASCII
 * @param fqd	fastq_data object pointer
 * @return	probability
 */
inline double error_prob(fastq_data *fqd, char q) {
	return exp(- (q + fqd->min_quality - 33) / 10. * log(10.));
} /* error_prob */

/**
 * Print out true error probabilities.
 *
 * @param fp	FILE pointer
 * @param fqd	fastq_data object pointer
 */
inline void fprint_error_probs(FILE *fp, fastq_data *fqd)
{
	unsigned int n_quality = fqd->max_quality - fqd->min_quality + 1;
	for (unsigned int q = 0; q < n_quality; ++q)
		fprintf(fp, " %f", error_prob(fqd, q));
} /* fprint_error_probs */

/**
 * Return length of given read.
 *
 * @param fqd	fastq_data object pointer
 * @param i	index of read
 * @return	length of read
 */
inline unsigned int read_length(fastq_data *fqd, unsigned int i)
{
	return fqd->n_lengths ? fqd->n_lengths[i] : fqd->n_max_length;
} /* read_length */

/**
 * Write one read of a given length as a row in R table format.
 *
 * @param fp	file handle pointer
 * @param read	read to write
 * @param len	length of read
 */
inline void write_read_in_table(FILE *fp, unsigned char *read, unsigned int len)
{
	for (unsigned int j = 0; j < len; ++j) {
		if (j)
			fprintf(fp, " ");
		fprintf(fp, "%u", (unsigned int) read[j]);
	}
	fprintf(fp, "\n");
} /* write_read_in_table */

#endif
