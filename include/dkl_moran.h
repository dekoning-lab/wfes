#include "dkl.h"

#ifndef DKL_MORAN_H
#define DKL_MORAN_H

void moran_row(double *row, DKL_INT i, double Ne, double u, double v, double w_AA, double w_Aa, double w_aa);
csr_sparse_matrix *moran_matrix_csr(wf_parameters *wf);

#endif /* DKL_MORAN_H */
