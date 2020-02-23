#include "util.h"
#include "error.h"

int fread_size_ts(FILE *fp, size_t *idx, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		if (fscanf(fp, "%lu", &idx[i]) != 1)
			return(mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"expecting %u elements in index file\n", n));

	return NO_ERROR;
} /* fread_size_ts */

int fread_uints(FILE *fp, unsigned int *id, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		if (fscanf(fp, "%u", &id[i]) != 1)
			return(mmessage(ERROR_MSG, FILE_FORMAT_ERROR,
				"expecting %u elements in index file\n", n));

	return NO_ERROR;
} /* fread_uints */

int read_uints(char const * const filename, unsigned int *id, size_t nreads)
{
	FILE *fp = fopen(filename, "r");

	if (!fp)
		return(mmessage(ERROR_MSG, FILE_OPEN_ERROR, filename));

	int err = fread_uints(fp, id, nreads);

	fclose(fp);
	return err;
} /* read_uints */
