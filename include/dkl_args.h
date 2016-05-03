#include "dkl.h"

#ifndef DKL_ARGS_H
#define DKL_ARGS_H

void __args_set_field(size_t argc, char **argv, size_t n_keywords,
                      char **keywords, bool required, void *field_ptr,
                      char type);

bool dkl_args_parse_flag(size_t argc, char **argv, bool required, ...);

double dkl_args_parse_double(size_t argc, char **argv, bool required, ...);

int64_t dkl_args_parse_int(size_t argc, char **argv, bool required, ...);

char *dkl_args_parse_string(size_t argc, char **argv, bool required, ...);

#endif /* DKL_ARGS_H */
