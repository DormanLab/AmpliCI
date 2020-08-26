/**
 * @file hash.h
 * @ Author Xiyu Peng
 * implement of hash table for DNA sequence
 */
#ifndef __H_HASH__
#define __H_HASH__

#include <string.h>
#include <stddef.h>
#include <stdlib.h>

#include "uthash.h"

typedef struct _hash_dna hash;

struct _hash_dna {
    unsigned char *sequence;	/*<! sequence */
    unsigned int count;		/*<! sequence abundance */
    size_t idx;			/*<! index of first sequence in data struct */
    size_t *idx_array;		/*<! indices of all reads */
    UT_hash_handle hh;
};

int add_sequence(hash **seq_count, unsigned char *seq, unsigned int length,size_t idx,int* err);
int add_seq_idx(hash *seq_count, unsigned char *seq, unsigned int length, size_t idx);
void delete_all(hash **seq_count);
unsigned int count_sequences(hash *seq_count, unsigned int k);
int count_sort(hash *a, hash *b);
void sort_by_count(hash **seq_count);
int store_index(hash *seq_count, unsigned int* uniq_id, size_t *index, unsigned int length, size_t sample_size);
int store_count(hash *seq_count, unsigned int *count, unsigned int length);
int find_index(hash *seq_count, unsigned char *seq, unsigned int length, size_t** idx_array);

#endif

