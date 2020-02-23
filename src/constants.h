/**
 * @file constants.h
 * @author Karin S. Dorman
 * @date Sun Jun 25 07:43:21 CDT 2017
 *
 * !!!!WARNING!!!!
 * You must go to the end of this file and set the typedef data_t to match
 * the type of data you would like to cluster.
 * !!!!WARNING!!!!
 */

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <inttypes.h>

#define UNUSED(x) (void)(x)

/* suss out C standard enforced */
#if defined(__STDC__)
#	define C89
#	if defined(__STDC_VERSION__)
#		define C90
#		if (__STDC_VERSION__ >= 199409L)
#			define C94
#		endif
#		if (__STDC_VERSION__ >= 199901L)
#			define C99
#		endif
#		if (__STDC_VERSION__ >= 201112L)
#			define C11
#		endif
#	endif
#endif

#ifndef C99
#define inline
#endif

/*
#define SIZE_T int
#define SIZE_T_FMT "%d"
*/
#include <stddef.h>
#define SIZE_T size_t
#define SIZE_T_FMT "%zu"

/* !!!!WARNING!!! Set this to match the data you would like to cluster. */
#ifndef RUN_KMODES
typedef uint8_t data_t;
#define AMPLICLUST_PRIu_data_t PRIu8
#define AMPLICLUST_SCNu_data_t SCNu8
#else
typedef uint32_t data_t;
#define AMPLICLUST_PRIu_data_t PRIu32
#define AMPLICLUST_SCNu_data_t SCNu32
#endif

#endif
