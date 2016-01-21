#include "dkl.h"

#ifndef DKL_MEMORY_H
#define DKL_MEMORY_H

#define dkl_alloc(SIZE, TYPE) malloc(sizeof(TYPE) * (SIZE))
#define dkl_new(TYPE) malloc(sizeof(TYPE))
#define dkl_size_alloc(SIZE, TYPESIZE) malloc((TYPESIZE) * (SIZE))
#define dkl_dealloc(PTR) free((PTR))
#define dkl_del(PTR) free((PTR))
#define dkl_realloc(PTR, SIZE, TYPE) realloc((PTR), (SIZE) * sizeof(TYPE))

#endif /* DKL_MEMORY_H */
