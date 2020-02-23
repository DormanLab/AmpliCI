/**
 * @file cmdline.h
 * @author Karin Dorman, kdorman@iastate.edu
 * @date Mon Dec 10 16:20:59 CST 2012
 *
 * Header file for cmdline.c.
 */

#ifndef __H_CMDLINE__
#define __H_CMDLINE__

#include <errno.h>
#include <stdio.h>
#include <string.h>

int usage_error(const char **argv, int i, void *obj);

/* check single entries */
int is_numeric(const char *argv);

/* read single entries */
int read_int(int argc, const char **argv, int i, void *obj);
unsigned int read_uint(int argc, const char **argv, int i, void *obj);
long read_long(int argc, const char **argv, int i, void *obj);
unsigned long read_ulong(int argc, const char **argv, int i, void *obj);
double read_cmdline_double(int argc, const char **argv, int i, void *obj);
char read_char(int argc, const char **argv, int i, void *obj);

/* read multiple enries */
unsigned int read_cmdline_doubles(unsigned int argc, const char **argv, unsigned int i, double **dret, void *obj);
unsigned int read_ulongs(unsigned int argc, const char **argv, unsigned int i, unsigned long **sret, void *obj);
unsigned int read_cmdline_uints(unsigned int argc, const char **argv, unsigned int i, unsigned int **uiret, void *obj);
unsigned int read_cmdline_strings(unsigned int argc, const char ** argv, unsigned int i, char const ***sret, void *obj);

void print_usage(const char *);
void fprint_usage(FILE *fp, const char *program_name, void *obj);

#endif
