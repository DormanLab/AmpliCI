#ifndef _ERROR_H_
#define _ERROR_H_

#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#ifdef USE_CURSES 
#include <curses.h>
#endif
#include "error.h"

/** Types of messages.
 */
enum {	NO_MSG,		/*!< no message */
	INFO_MSG,	/*!< informative message */
	DEBUG_MSG,	/*!< debugging message */
	WARNING_MSG,	/*!< warning message */
	ERROR_MSG	/*!< error message */
};

/** Types of errors.
 */
enum {	NO_ERROR,		/*!< no error */
	CUSTOM_ERROR,		/*!< customer error */
	NO_DATA,		/*!< necessary data not provided */
	MEMORY_ALLOCATION,	/*!< memory allocation error */
	FILE_NOT_FOUND,		/*!< file not found */
	FILE_OPEN_ERROR,	/*!< file open error */
	END_OF_FILE,		/*!< premature end of file error */
	FILE_FORMAT_ERROR,	/*!< invalid file format error */
	INVALID_CMDLINE,	/*!< invalid command line, specific below */
	INVALID_CMD_OPTION,	/*!< invalid command-line option */
	INVALID_CMD_ARGUMENT,	/*!< invalid argument to command-line option */
	INVALID_USER_INPUT,	/*!< invalid user setup */
	INTERNAL_MISMATCH,	/*!< data does not match in some way it must */
	INTERNAL_ERROR,		/*!< internal error */
	CLUSTER_SIZE_OVERFLOW,	/*!< overflow in no. of clusters? */
	STATE_SPACE_OVERFLOW,	/*!< memory overflow due to state space size */
	OUT_OF_TIME,		/*!< ran out of time */
	MEMORY_USAGE_LIMIT,	/*!< ran out of memory */
	MEMCPY_ERROR,		/*!< error in call to memcpy() */
	EXCEED_ITERATIONS,	/*!< exceed some iteration limit */
	NUM_ERRORS		/*!< total number of errors */
};

/** Level of verbosity.
 */
enum {	ABSOLUTE_SILENCE,	/*!< only output through files */
	SILENT,			/*!< no verbosity; final output only */
	QUIET,			/*!< try to be quiet */
	MINIMAL,		/*!< minimal verbosity */
	RESTRAINED,		/*!< some output */
	TALKATIVE,		/*!< talkative output */
	VERBOSE,		/*!< verbose output */
	DEBUG_I,		/*!< debugging output */
	DEBUG_II,		/*!< debugging output */
	DEBUG_III,		/*!< debugging output */
	DEBUG_OVERRIDE		/*!< ignores global level */
};

extern int global_debug_level;

#define CHECK_TIME(start_time, time_lim, file_name, fxn_name) {                \
	errno = NO_ERROR;                                                      \
	if ((time_lim) && difftime(time(NULL), (start_time)) > (time_lim)) {   \
		message(stderr, (file_name), (fxn_name), __LINE__, ERROR_MSG,  \
			OUT_OF_TIME, "", (time_lim));                          \
		errno = OUT_OF_TIME;                                           \
	}                                                                      \
}

#define CHECK_MEMORY(s_info, mem_lim, mem_new, file_name, fxn_name) {         \
	errno = NO_ERROR;                                                     \
	if ((mem_lim) && !sysinfo(&(s_info))) {                               \
		unsigned long mem_b = ((s_info).totalram - (s_info).freeram)  \
			* (s_info).mem_unit + (mem_new);                      \
		if (mem_b > (mem_lim)) {                                      \
			errno = MEMORY_USAGE_LIMIT;                           \
			message(stderr, (file_name), (fxn_name), __LINE__,    \
				ERROR_MSG, errno, "request for %.1"           \
				DOUBLE_F_NFMT "Gb exceeds %.1" DOUBLE_F_NFMT  \
				"Gb limit\n", (DOUBLE) mem_b / 1e9,           \
				(DOUBLE) (mem_lim) / 1e9);                    \
		}                                                             \
	}                                                                     \
}


/**
 * Print a formatted message to stderr.
 */
#define mmessage(type, err, ...) message(stderr, __FILE__, __func__,  __LINE__, (type), (err),  __VA_ARGS__)
/**
 * Print a formatted message to curses window.
 */
#ifdef USE_CURSES 
#define mwmessage(wp, type, err, ...) wmessage((wp), __FILE__, __func__,  __LINE__, (type), (err),  __VA_ARGS__)
#endif
/**
 * Print a formatted message to either stderr or curses window.
 */
#ifdef USE_CURSES
#define AMPICLI_MESSAGE(wp, type, err, ...) (                                  \
	(wp) ? wmessage((wp), __FILE__, __func__, __LINE__, (type), (err),     \
		__VA_ARGS__) :                                                 \
	message(stderr, __FILE__, __func__, __LINE__, (type), (err),           \
		__VA_ARGS__))
#else
#define AMPLICI_MESSAGE(wp, type, err, ...) (                                  \
	message(stderr, __FILE__, __func__, __LINE__, (type), (err),           \
		__VA_ARGS__))
#endif

#ifdef USE_CURSES
#define AMPLICI_MESSAGE_CONT(wp, ...) (                                        \
	(wp) ? wprintw((wp), __VA_ARGS__) :                                    \
	fprintf(stderr, __VA_ARGS__))
#else
#define AMPLICI_MESSAGE_CONT(wp, ...) (                                        \
	fprintf(stderr, __VA_ARGS__))
#endif

/**
 * Conditionally print a formatted message to stderr.
 */
#define debug_msg(level, fxn_debug_level, ...) do {                                  \
		if ((level) <= (fxn_debug_level) || (level) <= global_debug_level)              \
		message(stderr, __FILE__, __func__, __LINE__, level >= DEBUG_I \
		? DEBUG_MSG : INFO_MSG, NO_ERROR, __VA_ARGS__);                \
} while (0)

/**
 * Conditionally print a continuing message to stderr.
 */
#define debug_msg_cont(level, fxn_debug_level, ...) do {                             \
	if ((level) <= (fxn_debug_level) || (level) <= global_debug_level)        \
		fprintf(stderr, __VA_ARGS__);                                  \
} while(0)

/**
 * Conditionally print a formatted message to stderr or curses window.
 */
#ifdef USE_CURSES
#define cc_msg(wp, condition, lvl, msg, ...) do {                              \
	if ((condition) || ((lvl) && (lvl) <= global_debug_level)) (           \
		(wp)                                                           \
		? wmessage((wp), __FILE__, __func__, __LINE__, (lvl) >= DEBUG_I\
			? DEBUG_MSG : INFO_MSG, NO_ERROR, msg, __VA_ARGS__) :  \
		message(stderr, __FILE__, __func__, __LINE__, (lvl)>=DEBUG_I   \
			? DEBUG_MSG : INFO_MSG, NO_ERROR, msg, __VA_ARGS__));  \
} while (0)
#else
#define cc_msg(wp, condition, lvl, msg, ...) do {                              \
	if ((condition) || ((lvl) && (lvl) <= global_debug_level)) (           \
		message(stderr, __FILE__, __func__, __LINE__, (lvl)>=DEBUG_I   \
			? DEBUG_MSG : INFO_MSG, NO_ERROR, msg, __VA_ARGS__));  \
} while (0)
#endif

/**
 * Conditionally print a continuing message to stderr or curses window.
 */
#ifdef USE_CURSES
#define cc_msg_cont(wp, condition, level, ...) do {                            \
	if ((condition) || ((level) && (level) <= global_debug_level)) (       \
		(wp)                                                           \
		? wprintw((wp), __VA_ARGS__) :                                 \
		fprintf(stderr, __VA_ARGS__));                                 \
} while(0)
#else
#define cc_msg_cont(wp, condition, level, ...) do {                            \
	if ((condition) || ((level) && (level) <= global_debug_level)) (       \
		fprintf(stderr, __VA_ARGS__));                                 \
} while(0)
#endif

int message(FILE *, const char *, const char *, int, int, int, const char *, ...);
#ifdef USE_CURSES 
int wmessage(WINDOW *, const char *, const char *, int, int, int, const char *, ...);
#endif 

/**
 * Conditionally call a function.
 */
#define debug_call(condition, level, fxn_call) do {                            \
	if ((condition) || ((level) && (level) <= global_debug_level))         \
		(fxn_call);                                                    \
} while (0)

#endif
