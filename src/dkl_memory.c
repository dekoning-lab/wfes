#include "dkl_memory.h"

void *__dkl_alloc(int64_t n_elements, size_t type_size) {
  void *obj = calloc(n_elements, type_size);
  if (!obj) {
    error_print("Failed to allocate %.3f MB of memory",
                (double)(n_elements * type_size) / 1024.0);
  }
  return obj;
}


void __dkl_free_v(double *d) { free(d); }


void print_double_matrix(double *A, DKL_INT size) {
  printf("[");
  for(DKL_INT i = 0; i < size; i ++) {
    for(DKL_INT j = 0; j < size - 1; j ++) {
      DKL_INT idx = (i * size) + j;
      printf("%e, ", A[idx]);
    }
    printf("%e\n", A[(i * size) + (size - 1)]);
  }
  printf("\n]");
}

void print_double_vector(double *A, DKL_INT size) {
  printf("[");
  for(DKL_INT j = 0; j < size - 1; j ++) {
    printf("%e, ", A[j]);
  }
  printf("%e]\n", A[size - 1]);
}
