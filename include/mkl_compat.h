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

// PARDISO-specific parameters - refer to
// https://software.intel.com/en-us/node/521690

// Matrix types - mtype
#define MKL_PARDISO_MATRIX_TYPE_REAL_STRUCT_SYMMETRIC 1
#define MKL_PARDISO_MATRIX_TYPE_REAL_SYMMETRIC_POS_DEF 2
#define MKL_PARDISO_MATRIX_TYPE_REAL_SYMMETRIC_INDEF -2
#define MKL_PARDISO_MATRIX_TYPE_COMPLEX_STRUCT_SYMMETRIC 3
#define MKL_PARDISO_MATRIX_TYPE_COMPLEX_HERMITIAN_POS_DEF 4
#define MKL_PARDISO_MATRIX_TYPE_COMPLEX_HERMITIAN_INDEX -4
#define MKL_PARDISO_MATRIX_TYPE_COMPLEX_SYMMETRIC 6
#define MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC 11
#define MKL_PARDISO_MATRIX_TYPE_COMPLEX_UNSYMMETRIC 13

// Solver phases - phase
#define MKL_PARDISO_SOLVER_PHASE_ANALYSIS 11
#define MKL_PARDISO_SOLVER_PHASE_ANALYSIS_NUMERICAL_FACTORIZATION 12
#define MKL_PARDISO_SOLVER_PHASE_ANALYSIS_NUMERICAL_FACTORIZATION_SOLVE 13
#define MKL_PARDISO_SOLVER_PHASE_NUMERICAL_FACTORIZATION 22
#define MKL_PARDISO_SOLVER_PHASE_NUMERICAL_FACTORIZATION_SOLVE 23
#define MKL_PARDISO_SOLVER_PHASE_SOLVE_ITERATIVE_REFINEMENT 33
#define MKL_PARDISO_SOLVER_PHASE_SOLVE_FORWARD_SUBSTITUTION 331
#define MKL_PARDISO_SOLVER_PHASE_SOLVE_DIAGONAL_SUBSTITUTION 332
#define MKL_PARDISO_SOLVER_PHASE_SOLVE_BACKWARD_SUBSTITUTION 333
#define MKL_PARDISO_SOLVER_PHASE_RELEASE_MEMORY_MNUM 0
#define MKL_PARDISO_SOLVER_PHASE_RELEASE_MEMORY_ALL -1

// PARDISO IPARM - refer to https://software.intel.com/en-us/node/521691
#define MKL_PARDISO_NA 100
#define MKL_PARDISO_DEFAULT 0
#define MKL_PARDISO_FALSE 1
#define MKL_PARDISO_REPORT_ENABLE -1
#define MKL_PARDISO_REPORT_DISABLE 0

#define MKL_PARDISO_DEFAULT_SETTINGS 0

#define MKL_PARDISO_FILL_IN_REDUCING_ORDERING_OPTION 1
#define MKL_PARDISO_FILL_IN_REDUCING_ORDERING_MIN_DEGREE 0
#define MKL_PARDISO_FILL_IN_REDUCING_ORDERING_NESTED_DISSECTION_METIS 2
#define MKL_PARDISO_FILL_IN_REDUCING_ORDERING_NESTED_DISSECTION_OMP 3

#define MKL_PARDISO_PRECONDITIONED_CGS_CG 3
#define MKL_PARDISO_PRECONDITIONED_DEFAULT 0

#define MKL_PARDISO_PERMUTATION_OPTION 4
#define MKL_PARDISO_PERMUTATION_USE_PERM 1
#define MKL_PARDISO_PERMUTATION_RETURN 2

#define MKL_PARDISO_RETURN_OPTION 5
#define MKL_PARDISO_RETURN_OVERRIDE 1

#define MKL_PARDISO_ITERATIVE_STEPS_USED 6

#define MKL_PARDISO_ITERATIVE_REFINEMENT_MAX 7

#define MKL_PARDISO_PIVOTING_PERTURBATION 9

#define MKL_PARDISO_SCALING_OPTION 10
#define MKL_PARDISO_SCALING_ENABLE 1

#define MKL_PARDISO_SOLVE_OPTION 11
#define MKL_PARDISO_SOLVE_CONJUGATE_TRANSPOSED 1
#define MKL_PARDISO_SOLVE_TRANSPOSED 2

#define MKL_PARDISO_WEIGHTED_MATCHING_OPTION 12
#define MKL_PARDISO_WEIGHTED_MATCHING_DISABLE 0
#define MKL_PARDISO_WEIGHTED_MATCHING_ENABLE 1

#define MKL_PARDISO_N_PERTURBED_PIVOTS 13

#define MKL_PARDISO_PEAK_MEMORY_SYMBOLIC_PHASE 14
#define MKL_PARDISO_PERMANENT_MEMORY_SYMBOLIC_PHASE 15
#define MKL_PARDISO_PEAK_MEMORY_NUMERIC_PHASE 16

#define MKL_PARDISO_REPORT_NNZ_FACTORS 17
#define MKL_PARDISO_REPORT_FLOP_FACTOR_PHASE 18
#define MKL_PARDISO_REPORT_CGS_CG_DIAGNOSTIC 19

#define MKL_PARDISO_PIVOT_SYMMETRIC_INDEF_OPTION 20
#define MKL_PARDISO_PIVOT_SYMMETRIC_INDEF_DIAG 0
#define MKL_PARDISO_PIVOT_SYMMETRIC_INDEF_BUNCH_KAUFMAN 1
#define MKL_PARDISO_PIVOT_SYMMETRIC_INDEF_DIAG_NO_AUTO 2
#define MKL_PARDISO_PIVOT_SYMMETRIC_INDEF_BUNCH_KAUFMAN_NO_AUTO 3

#define MKL_PARDISO_POS_EIGENVALUES 21
#define MKL_PARDISO_NEG_EIGENVALUES 22

#define MKL_PARDISO_PARALLEL_FACTORIZATION_OPTION 23
#define MKL_PARDISO_PARALLEL_FACTORIZATION_TWO_LEVEL 1

#define MKL_PARDISO_PARALLEL_SOLVE_OPTION 24
#define MKL_PARDISO_SEQUENTIAL_SOLVE 1

#define MKL_PARDISO_MATRIX_CHECK_OPTION 26
#define MKL_PARDISO_MATRIX_CHECK_ENABLE 1

#define MKL_PARDISO_PRECISION_OPTION 27
#define MKL_PARDISO_PRECISION_DOUBLE 0
#define MKL_PARDISO_PRECISION_SINGLE 1

#define MKL_PARDISO_ZERO_NEG_PIVOTS 29

#define MKL_PARDISO_PARTIAL_SOLVE_OPTION 30
#define MKL_PARDISO_PARTIAL_SOLVE_SPARSE_RHS_PARTIAL 1
#define MKL_PARDISO_PARTIAL_SOLVE_SPARSE_RHS_FULL 2
#define MKL_PARDISO_PARTIAL_SOLVE_SELECTED 3

#define MKL_PARDISO_OPTIMAL_THREADS 33

// This is not a mistake!
#define MKL_PARDISO_INDEXING_OPTION 34
#define MKL_PARDISO_INDEXING_ZERO 1
#define MKL_PARDISO_INDEXING_ONE 0

#define MKL_PARDISO_SCHUR_COMPLEMENT_OPTION 35
#define MKL_PARDISO_SCHUR_COMPLEMENT_SOLVE 1
#define MKL_PARDISO_SCHUR_COMPLEMENT_PARTIAL 2

#define MKL_PARDISO_MATRIX_FORMAT_OPTION 36
#define MKL_PARDISO_MATRIX_FORMAT_CSR 0
#define MKL_PARDISO_MATRIX_FORMAT_BSR 1
#define MKL_PARDISO_MATRIX_FORMAT_CONVERT_BSR -1

#define MKL_PARDISO_PIVOT_OPTION 55
#define MKL_PARDISO_PIVOT_CALLBACK 1

#define MKL_PARDISO_OOC_OPTION 59
#define MKL_PARDISO_OOC_OVERFLOW 1
#define MKL_PARDISO_OOC_ALWAYS 2

#define MKL_PARDISO_OOC_MIN_RAM 62

#endif /* MKL_COMPAT_H */
