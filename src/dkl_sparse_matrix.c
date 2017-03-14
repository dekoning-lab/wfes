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

double *csr_to_dense(csr_sparse_matrix *A) {
  // Copy the buffer in the CSR representation into Q
  // https://software.intel.com/en-us/node/520848

  double *buffer = dkl_alloc(A->nrows * A->ncols, double);
  DKL_INT info;
  DKL_INT *job = dkl_alloc(6, DKL_INT);
  job[0] = 1; // the rectangular matrix A is restored from the CSR format
  job[1] = 0; // zero-based indexing for the rectangular matrix A is used
  job[2] = 0; // zero-based indexing for the matrix in CSR format is used
  job[3] = 2; // `buffer` is a whole matrix A
  job[4] = A->nnz; // maximum number of the non-zero elements allowed if job[0]=0

  mkl_ddnscsr(job, &A->nrows, &A->ncols, buffer, &A->ncols, A->data, A->cols, A->row_index, &info);
  assert(info == 0);
  dkl_dealloc(job);

  return buffer;
}

// Copy xth column of A into y
void csr_matrix_column(csr_sparse_matrix *A, DKL_INT x, double* y) {
  for (DKL_INT i = 0; i < A->nrows; i++) {
    DKL_INT idx = A->row_index[i];          // Start of current row
    DKL_INT n = A->row_index[i + 1] - idx;  // Number of entries in current row
    bool found = false;                     // Is this element stored?
    for(DKL_INT j = 0; j < n; j++) {
      DKL_INT col = A->cols[idx + j];
      if (col == x) {
        found = true;
        y[i] = A->data[idx + j];
        break;
      }
    }
    if (!found) {
      y[i] = 0;
    }
  }
}
