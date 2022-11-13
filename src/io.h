#ifndef __IO_H__
#define __IO_H__

#include <stdio.h>	/* FILE, fprintf() */
#ifdef USE_CURSES
#include <curses.h>	/* WINDOW, wprintw() */
#endif

/**
 * Macro forward to next instance of selected char or EOF.
 */
#define fforward(f, c, t) do {                                                 \
	(c) = fgetc((f));                                                      \
} while ((c) != EOF && (c) != (t));

/**
 * Macro forward to next instance of selected char or EOF, counting characters
 * consumed.
 */
#define fforward_cnt(f, c, t, a)                                               \
while (((c) = fgetc((f))) != EOF && (c) != (t)) {                              \
	(a)++;                                                                 \
}

/**
 * Macro read next max n chars up to selected char.
 */
#define next_word(f, s, n, t) {                                                \
	unsigned int i = 0;                                                    \
	while (i < (n) && ((s)[i] = fgetc((f))) != EOF && (s)[i++] != t);      \
        s[i] = 0;                                                              \
}


/**
 * Seamlessly handle presence/absence of ncurses window.
 */
#ifdef USE_CURSES
#define PRINT(wp, ...) if ((wp)) wprintw((wp), __VA_ARGS__); else fprintf(stderr, __VA_ARGS__);
#else
#define PRINT(wp, ...) fprintf(stderr, __VA_ARGS__)
#endif

/* print vectors */
void fprint_doubles(FILE *fp, double *v, size_t len, int precision, int newline);
void fprint_uints(FILE *fp, unsigned int *v, size_t len, int width, int newline);
void fprint_size_ts(FILE *fp, size_t *v, size_t len, int width, int newline);

#ifdef USE_CURSES
void wprint_doubles(WINDOW *wp, double *v, size_t len, int precision, int newline);
void wprint_uints(WINDOW *fp, unsigned int *v, size_t len, int width, int newline);
void wprint_size_ts(WINDOW *fp, size_t *v, size_t len, int width, int newline);
#endif

/* read vectors */

/* print matrices */
void fprint_vectorized_sq_matrix(FILE *fp, double *mat, size_t n, int row);
void fprint_vectorized_matrix(FILE *fp, double *mat, size_t n, size_t l, int row);
void fprint_vectorized_uintmatrix(FILE *fp, unsigned int *mat, unsigned int n, unsigned int l,int row);

#ifdef USE_CURSES
void wprint_vectorized_sq_matrix(WINDOW *fp, double *mat, size_t n, int row);
void wprint_vectorized_matrix(WINDOW *wp, double *mat, size_t n, size_t l, int row);
#endif

#endif
