#include "dkl.h"

#ifndef DKL_SPARSE_MATRIX_H
#define DKL_SPARSE_MATRIX_H

/**
 * csr_sparse_matrix - compressed sparse row matrix representation (3-array intel variant)
 * Implemented according to: https://software.intel.com/en-us/node/599882
 */
typedef struct csr_sparse_matrix {
  DKL_INT nrows;      // Number of rows
  DKL_INT ncols;      // Number of columns
  DKL_INT nnz;        // Number of non-zero elements
  double *data;       // CSR-stored values
  DKL_INT *cols;      // Column of each cell in data
  DKL_INT *row_index; // Index of the first element in each row
} csr_sparse_matrix;

/**
 * csr_sparse_check: assert that the sparse matrix is in the correct CSR format
 *representation
 *
 * @param[in] A sparse matrix structure pointer
 *
 * @return true if the matrix is correct, false otherwise
 *
 */
bool csr_sparse_is_correct(csr_sparse_matrix *A);

#endif /* DKL_SPARSE_MATRIX_H */
