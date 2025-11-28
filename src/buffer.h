#ifndef __BUFFER_H__
#define __BUFFER_H__

//#include <stddef.h>
#include <stdlib.h>

#define allocate_buffer(buf, buf_len, err) do {                                \
	if (!(buf))                                                            \
		(buf) = malloc((buf_len) * sizeof(*(buf)));                    \
	if (!(buf))                                                            \
		(err) = MEMORY_ALLOCATION;                                     \
} while (0);

#define callocate_buffer(buf, buf_len, err) do {                               \
	if (!(buf))                                                            \
		(buf) = calloc((buf_len), sizeof(*(buf)));                     \
	if (!(buf))                                                            \
		(err) = MEMORY_ALLOCATION;                                     \
} while (0);

/**
 * (Alternative) Macro buffer expansion.
 */
#define reallocate_buffer(buf_ptr, buf_len, buf, err) do {                     \
	if ((buf_ptr) == (buf) + buf_len) {                                    \
		(buf_ptr) = realloc((buf), 2 * (buf_len) * sizeof(*buf));      \
		if ((buf_ptr)) {                                               \
			(buf) = (buf_ptr);                                     \
			(buf_ptr) += buf_len;                                  \
			(buf_len) *= 2;                                        \
		} else {                                                       \
			(err) = MEMORY_ALLOCATION;                             \
		}                                                              \
	}                                                                      \
} while (0);

#endif
