#include "dkl.h"

#ifndef DKL_MEMORY_H
#define DKL_MEMORY_H

void *__dkl_alloc(int64_t n_elements, size_t type_size);

#define dkl_alloc(NMEMB, TYPE) __dkl_alloc((NMEMB), sizeof(TYPE))
#define dkl_new(TYPE) __dkl_alloc(1, sizeof(TYPE))
#define dkl_size_alloc(NMEMB, TYPESIZE) __dkl_alloc((NMEMB), (TYPESIZE))

#define dkl_dealloc(PTR) free((PTR))
#define dkl_del(PTR) free((PTR))
#define dkl_realloc(PTR, SIZE, TYPE) realloc((PTR), (SIZE) * sizeof(TYPE))

#define check_mem(A)                                                           \
  if (!(A)) {                                                                  \
    error_print("Memory not allocated");                                       \
  }

#endif /* DKL_MEMORY_H */
