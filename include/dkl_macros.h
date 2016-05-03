#include "dkl.h"

#ifndef DKL_MACROS_H
#define DKL_MACROS_H

#define streq(A,B) strcmp((A), (B)) == 0
#define strneq(A,B) strcmp((A), (B)) != 0
#define println(A, ...) printf(A "\n", ##__VA_ARGS__)
#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define error_print(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)



// Macro definitions
#ifdef DEBUG
#define debug_assert(EXP)                       \
  do {                                          \
    assert((EXP));                              \
  } while (0)
#else
#define debug_assert(EXP)                       \
  do {                                          \
  } while (0)
#endif

#endif /* DKL_MACROS_H */
