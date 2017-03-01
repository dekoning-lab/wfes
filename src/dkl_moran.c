#include "dkl_moran.h"

void moran_row(double *row, DKL_INT i, double Ne, double u, double v, double w_AA, double w_Aa, double w_aa) {

  double p, q, p2, q2, pq, pq2, iu, iv, w_bar;
  double P_b_AA, P_b_Aa, P_b_aa, P_d_AA, P_d_Aa, P_d_aa, P_m2, P_m1, P_0, P_p1, P_p2;

  p = (double)(i / Ne);
  q = 1 - p;
  p2 = p * p;
  q2 = q * q;
  pq = p * q;
  pq2 = 2 * pq;
  iu = 1 - u;
  iv = 1 - v;
  w_bar = (w_AA * p2) + (w_Aa * pq2) + (w_aa * q2);

  P_b_AA = (pq2 * u) + (p2 * iv);
  P_b_Aa = (p2 * v) + (q2 * u) + (pq * iv) + (pq * iu);
  P_b_aa = (pq2 * v) + (q2 * iu);

  P_d_AA = w_AA * p2 / w_bar;
  P_d_Aa = w_Aa * pq2 / w_bar;
  P_d_aa = w_aa * q2 / w_bar;

  P_m2 = P_b_aa * P_d_AA;
  P_m1 = (P_b_Aa * P_d_AA) + (P_b_aa * P_d_Aa);
  P_p1 = (P_b_AA * P_d_Aa) + (P_b_Aa * P_d_aa);
  P_p2 = P_b_AA * P_d_aa;
  P_0 = 1 - P_m2 - P_m1 - P_p1 - P_p2;

  row[0] = P_m2;
  row[1] = P_m1;
  row[2] = P_0;
  row[3] = P_p1;
  row[4] = P_p2;
}

csr_sparse_matrix *moran_matrix_csr(wf_parameters *wf) {

  double Ne = (double)2 * wf->population_size;
  DKL_INT size = (2 * wf->population_size) - 1;

  csr_sparse_matrix *Q = dkl_alloc(1, csr_sparse_matrix);

  Q->nrows = size;
  Q->ncols = size;
  DKL_INT nnz = ((Ne - 3) * 5) + 10;
  Q->nnz = nnz;

  Q->data = dkl_alloc(nnz, double);
  Q->cols = dkl_alloc(nnz, DKL_INT);
  Q->row_index = dkl_alloc(Q->nrows + 1, DKL_INT);

  double *row = dkl_alloc(5, double);

  double u = wf->backward_mutation_rate;
  double v = wf->forward_mutation_rate;
  double w_AA = 1 + wf->selection;
  double w_Aa = 1 + (wf->selection * wf->dominance_coefficient);
  double w_aa = 1;

  // Head
  // [* . . . . .]
  // [* * * * . .]
  Q->cols[0] = 0;
  Q->cols[1] = 0;
  Q->cols[2] = 1;
  Q->cols[3] = 2;
  Q->cols[4] = 3;

  Q->row_index[0] = 0;
  Q->row_index[1] = 1;

  Q->data[0] = 1;
  moran_row(row, 1, Ne, u, v, w_AA, w_Aa, w_aa);
  row[2] = 1 - row[1] - row[3] - row[4];
  Q->data[1] = row[1];
  Q->data[2] = row[2];
  Q->data[3] = row[3];
  Q->data[4] = row[4];

  // Tail
  // [. . * * * *]
  // [. . . . . *]
  Q->cols[nnz - 1] = size - 1;
  Q->cols[nnz - 2] = size - 1;
  Q->cols[nnz - 3] = size - 2;
  Q->cols[nnz - 4] = size - 3;
  Q->cols[nnz - 5] = size - 4;

  Q->row_index[Q->nrows] = Q->nnz + 1; // Special CSR value
  Q->row_index[Q->nrows - 1] = Q->nnz - 1;
  Q->row_index[Q->nrows - 2] = Q->nnz - 5;

  Q->data[nnz - 1] = 1;
  moran_row(row, 1, Ne, u, v, w_AA, w_Aa, w_aa);
  row[2] = 1 - row[0] - row[1] - row[3];
  Q->data[nnz - 5] = row[0];
  Q->data[nnz - 4] = row[1];
  Q->data[nnz - 3] = row[2];
  Q->data[nnz - 2] = row[3];

  for(DKL_INT i = 2; i < size; i++) {
    DKL_INT l = (i - 1) * 5;
    Q->row_index[i] = l - 1;

    Q->cols[l] = i - 2;
    Q->cols[l + 1] = i - 1;
    Q->cols[l + 2] = i;
    Q->cols[l + 3] = i + 1;
    Q->cols[l + 4] = i + 2;

    moran_row(row, i, Ne, u, v, w_AA, w_Aa, w_aa);

    Q->data[l    ] = row[0];
    Q->data[l + 1] = row[1];
    Q->data[l + 2] = row[2];
    Q->data[l + 3] = row[3];
    Q->data[l + 4] = row[4];

  }
  return Q;
}
