#ifndef DKL_ERROR_H
#define DKL_ERROR_H

int64_t dkl_errno = 0; 

#define DKL_RUNTIME_ERROR 1
#define DKL_ARG_PARSE_ERROR 2
#define DKL_PARAM_ERROR 3
#define DKL_HELP_EXIT 4
#define DKL_OPTION_NOT_FOUND 5

#define dkl_clear_errno() dkl_errno = 0

#endif /* DKL_ERROR_H */
