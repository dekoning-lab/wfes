#include "dkl.h"

#ifndef DKL_ARGS_H
#define DKL_ARGS_H

void __args_set_field(int argc, char **argv, int64_t n_keywords, char **keywords, bool required, void *field_ptr, char type);

bool dkl_args_parse_flag(int argc, char **argv, bool required, ...);

double dkl_args_parse_double(int argc, char **argv, bool required, ...);

int64_t dkl_args_parse_int(int argc, char **argv, bool required, ...);

char *dkl_args_parse_string(int argc, char **argv, bool required, ...);

#endif /* DKL_ARGS_H */
