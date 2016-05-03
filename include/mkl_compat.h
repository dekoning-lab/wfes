#include "dkl.h"
#ifndef MKL_COMPAT_H
#define MKL_COMPAT_H

#define MKL_INTERFACE_ILP64 1
#define MKL_THREADING_INTEL 0
#define MKL_PEAK_MEM_ENABLE 1
#define MKL_PEAK_MEM 2

void vdExp(const DKL_INT, const double *, double *);
void cblas_dscal(const DKL_INT, const double a, double *x, const DKL_INT incx);
void mkl_ddnscsr(const DKL_INT *job, const DKL_INT *m, const DKL_INT *n,
                 double *adns, const DKL_INT *lda, double *acsr, DKL_INT *ja,
                 DKL_INT *ia, DKL_INT *info);
void pardiso_64(void *pt, DKL_INT *maxfct, DKL_INT *mnum, DKL_INT *mtype,
                DKL_INT *phase, DKL_INT *n, void *a, DKL_INT *ia, DKL_INT *ja,
                DKL_INT *perm, DKL_INT *nrhs, DKL_INT *iparm, DKL_INT *msglvl,
                void *b, void *x, DKL_INT *error);
DKL_INT MKL_Peak_Mem_Usage(int mode);

#endif /* MKL_COMPAT_H */
