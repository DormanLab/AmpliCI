/**
 * @file hash.c
 * @author Xiyu Peng
 * implemention of hash table for DNA sequences
 */

#include <string.h>
#include <stddef.h>
#include <stdlib.h>

#include "uthash.h"
#include "hash.h"
#include "error.h"

/**
 * Build hash table for sequences.
 *
 * @param seq_count	hash table
 * @param seq		sequence to add
 * @param length	length of sequence
 * @param idx		index of the sequence in data struct
 * 
 * @return		The first observation time
 */
int add_sequence(hash **seq_count, unsigned char *seq, unsigned int length,
								size_t idx, int *err)
{
	hash *new;
	int first = 0;

	HASH_FIND( hh, *seq_count, seq, length * sizeof *seq, new );

	if (!new) {
		new = (hash *) malloc(sizeof *new);
		if (!new) {
			*err = MEMORY_ALLOCATION;
			mmessage(ERROR_MSG, *err, "hash");
			return 0;
		}

		new->sequence = seq;
		new->count = 1;
		new->idx = idx;	
		new->idx_array = NULL;
		HASH_ADD_KEYPTR(hh, *seq_count, new->sequence,
						length * sizeof *seq, new);
		new->seeds = 0;
		first = 1;
	} else {
		new->count++;
	}

	return first;
}/* add_sequence */

/**
 * construct the index array in the hash table
 *
 * @param seq_count	hash table
 * @param seq		sequence to add
 * @param length	length of sequence
 * @param idx		index of the sequence in data struct
 * 
 * @return		error status
 */
int add_seq_idx(hash *seq_count, unsigned char *seq, unsigned int length,
								size_t idx)
{
	hash *unit;

	HASH_FIND(hh, seq_count, seq, length * sizeof *seq, unit);

	if (!unit)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Sequence not found "
					"in hash table ! Please check !");

	if (!unit->idx_array) {
		unit->idx_array = malloc(unit->count *sizeof *unit->idx_array);

		if (!unit->idx_array)
			return mmessage(ERROR_MSG, MEMORY_ALLOCATION,
							"hash::idx_array");
		
		unit->idx_array[0] = idx;		
		unit->count = 1; 
	} else {
		unit->idx_array[unit->count] = idx;
		unit->count++;
	}

	return NO_ERROR;
}/* add_seq_idx */

/**
 * Delete all entries in hash table and free memory at the same time.
 *
 * @param seq_count	hash table
 */
void delete_all(hash **seq_count){
	hash *tmp, *current,*hash_count;
	hash_count = *seq_count;

	HASH_ITER( hh, hash_count, current, tmp ) {
		if(current->idx_array) free(current->idx_array);
		HASH_DEL( hash_count, current );
		free(current);
	}
	*seq_count = NULL;
}/* delete_all */

/**
 * Count number of uniq sequences whose frequency are larger or equal than k.
 *
 * @param seq_count	hash table
 * @param k		count
 * @return		count of hash entries with count greater than k
 */
unsigned int count_sequences(hash *seq_count, unsigned int k) {
	unsigned int count = 0;
	hash *s;

	for (s = seq_count; s != NULL; s = s->hh.next)
		if (s->count >= k)
			count++;
	return count;
}/* count_sequences */

/**
 * compare the counts from two hash structs for sorting
 * 
 * @param a hash table
 * @param b hash table
 * 
 * @ return an int as the result of comparison
 * */
int count_sort(hash *a, hash *b)
{
	return (b->count - a->count);  // change for a order of decreasing
}/* count_sort */

/**
 * sort the hash tables based on the count
 * 
 * @param seq_count pointer to hash table
 * 
 * */
void sort_by_count(hash **seq_count)
{
	HASH_SORT(*seq_count, count_sort);
}/* sort_by_count */


/**
 * Fill an array with the index of first observed instance of each unique
 * sequence in hash table.
 * And also fill an array with index of each reads in unique sequence table.
 * 
 * 
 * @param seq_count hash table
 * @param index     pointer to uniq_seq_idx array
 * @param uniq_id   pointer to reads_uniq_id array
 * @param length    length of uniq_seq_isx array
 * @param sample_size length of reads_uniq_id array
 * 
 * @return		error status
 * 
 * */
int store_index(hash *seq_count, unsigned int* uniq_id, size_t *index, unsigned int length, size_t sample_size){
	hash *s;
	unsigned int i = 0;
	unsigned int idx;

	for (s = seq_count; s != NULL; s = s->hh.next) {
		if (i == length)
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
						"exceed the length");
		index[i] = s->idx;   // store idx of uniq seq in dmat

		/* store idx of read in uniq seq array */
		for (unsigned j = 0; j < s->count; j++){
			idx = s->idx_array[j];
			if(idx >= sample_size)
				return mmessage(ERROR_MSG,INTERNAL_ERROR,"exceed the sample_size");
			uniq_id[idx] = i;
		}

		i++;
	}
	return NO_ERROR;
}/* store_index */

/**
 * Fill an array with count of each unique sequence in hash table.
 * 
 * @param seq_count	hash table
 * @param count		pointer to array
 * @param length	length of array
 * 
 * @return		error status
 * */
int store_count(hash *seq_count, unsigned int *count, unsigned int length)
{

	hash *s;
	unsigned int i = 0;

	for (s = seq_count; s != NULL; s = s->hh.next) {
		if (i == length)
			return mmessage(ERROR_MSG, INTERNAL_ERROR,
							"exceed the length");
		count[i] = s->count;
		i++;
	}
	return NO_ERROR;
}/* store_count */

/**
 * Find the index array of the given sequence in hash table.
 * 
 * @param seq_count	hash table
 * @param seq		sequence
 * @param length	length of the sequence
 * @idx_array       address storing the pointer of the array
 * 
 * @return		error status
 * */
int find_index(hash *seq_count, unsigned char *seq, unsigned int length, size_t** idx_array){
	
	hash *unit = NULL;

	HASH_FIND(hh, seq_count, seq, length * sizeof *seq, unit);
	if (!unit)
		return mmessage(ERROR_MSG, INTERNAL_ERROR, "Do not find the sequence in hash table");
	
	*idx_array = unit->idx_array;

	return NO_ERROR;
}/* find_index */

