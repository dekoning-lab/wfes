#include "dkl_sparse_matrix.h"

bool csr_sparse_is_correct(csr_sparse_matrix *A) {
  DKL_INT i, j;
  for (i = 0; i < A->nrows - 1; i++) {
    if (A->row_index[i] >= A->row_index[i + 1]) {
      printf("%" PRId64 " %" PRId64 "\n", A->row_index[i], A->row_index[i + 1]);
      printf("The row indexes are not sorted - failed at row %" PRId64 "\n", i);
      return false;
    }
    for (j = A->row_index[i]; j < A->row_index[i + 1] - 1; j++) {
      if (A->cols[j] >= A->cols[j + 1]) {
        printf("%" PRId64 " %" PRId64 " %g %g %g %g\n", A->cols[j],
               A->cols[j + 1], A->data[j], A->data[j + 1], A->data[j + 2],
               A->data[j + 3]);
        printf("The column indexes are not sorted - failed at col %" PRId64
               ", row %" PRId64 ", entry %" PRId64 "\n",
               A->cols[j], i, j);
        return false;
      }
    }
  }
  return true;
}
