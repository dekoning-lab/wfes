#include "dkl_memory.h"

void *__dkl_alloc(int64_t n_elements, size_t type_size) {
  void *obj = calloc(n_elements, type_size);
  if (!obj) {
    error_print("Failed to allocate %.3f MB of memory",
                (double)(n_elements * type_size) / 1024.0);
  }
  return obj;
}
