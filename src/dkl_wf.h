#include "dkl.h"

/**
 * wf_parameters - parameters of the Wright-Fisher model
 */
typedef struct wf_parameters {
  DKL_INT population_size;       // Effective population size
  double selection;              // Selective advantage of 'A' over 'a'
  double forward_mutation_rate;  // Mutation rate from 'A' to 'a'
  double backward_mutation_rate; // Mutation rate from 'a' to 'A'
  double dominance_coefficient;  // Proportion of selective advantage of 'Aa'
                                 // over 'aa'
} wf_parameters;

/**
 * wf_statistics - summary statistics of the Wright-Fisher model
 */
typedef struct wf_statistics {
  double probability_extinction; // Probability of absorption through exinction
  double probability_fixation;   // Probability of absorption through fixation
  double time_extinction; // Expected number of generations until extinction
  double time_fixation;   // Expected number of generations until fixation
  double count_before_extinction;   // Total expected count of 'A' in all
                                    // generations before extinction
  double *extinction_probabilities; // Probability of extinction vector,
                                    // conditional on the starting state
  double *fixation_probabilities;   // Probability of fixation vector
  double *generations;              // Sojourn time vector
} wf_statistics;

wf_statistics *wf_statistics_new(DKL_INT population_size) {
  wf_statistics *r = dkl_new(wf_statistics);

  r->probability_extinction = 0;
  r->probability_fixation = 0;
  r->time_extinction = 0;
  r->time_fixation = 0;
  r->count_before_extinction = 0;

  DKL_INT size = (2 * population_size) - 1;
  r->extinction_probabilities = dkl_alloc(size, double);
  r->fixation_probabilities = dkl_alloc(size, double);
  r->generations = dkl_alloc(size, double);

  return r;
}

void wf_statistics_del(wf_statistics *r) {
  dkl_del(r->extinction_probabilities);
  dkl_del(r->fixation_probabilities);
  dkl_del(r->generations);
  dkl_del(r);
}

/**
 * wf_sampling_coefficient: calculate the binomial sampling coefficient (psi)
 * for the Wright-Fisher matrix
 *
 * @param[in] N Population size (range: 2-Inf)
 * @param[in] s Selection coefficient (range: -1:Inf)
 * @param[in] u Forward mutation rate (range: 0:1/2N)
 * @param[in] v Backward mutation rate (range: 0:1/2N)
 * @param[in] h Dominance coefficient (range: 0:1)
 * @param[in] i Current number of copies of 'A' in the population
 *
 * @return The sampling coefficient (psi) for the Wright-Fisher matrix
 */
double wf_sampling_coefficient(wf_parameters *wf, DKL_INT i) {
  double a = (1 + wf->selection) * (i * i);
  double b = (1 + (wf->selection * wf->dominance_coefficient)) * i *
             ((2 * wf->population_size) - i);
  double c = ((2 * wf->population_size) - i) * ((2 * wf->population_size) - i);
  double q = (a + b) / (a + (2 * b) + c);
  return ((1 - wf->backward_mutation_rate) * q) +
         ((1 - q) * wf->forward_mutation_rate);
}

/**
 * wf_make_block: calculate 'block_size' rows of the Wright-Fisher matrix
 *
 * @param[in] N Population size (range: 2-Inf)
 * @param[in] s Selection coefficient (range: -1:Inf)
 * @param[in] u Forward mutation rate (range: 0:1/2N)
 * @param[in] v Backward mutation rate (range: 0:1/2N)
 * @param[in] h Dominance coefficient (range: 0:1)
 * @param[in] t Numeric zero threshold (range: 0:1e-10)
 * @param[in] row_offset The starting row of the block
 * @param[in] block_size The number of rows in the block
 * @param[in] buffer A memory buffer to hold the block
 * @param[in] binom Pre-calculated binomial coefficients (i=1)
 * @param[in, out] Q CSR matrix to which the block belongs
 * @param[in] current_size Current number of non-zero elements in Q
 * @param[in] global_nnz Number of non-zero elements in each row of Q
 * @param[in] old_size Number of non-zero elements in Q before the block in
 * built
 */
void wf_make_block(wf_parameters *wf, double threshold, DKL_INT row_offset,
                   DKL_INT block_size, double *buffer, double *binom,
                   csr_sparse_matrix *Q, DKL_INT *current_size,
                   DKL_INT *global_nnz, DKL_INT *old_size) {
  DKL_INT size = (2 * wf->population_size) - 1;
  DKL_INT i, j;
#pragma omp parallel for private(j)
  for (i = 0; i < block_size; i++) {

    DKL_INT local_nnz = 0;
    double lq = wf_sampling_coefficient(wf, row_offset + i + 1);
    double lnq = log(1.0 - lq);
    lq = log(lq);
    double lz = lq - lnq;
    double lp = binom[0] + lq + (lnq * size);
    double prev = lp;

    buffer[i * size] = (lp > MIN_LOG) ? lp : -INFINITY;

    for (j = 1; j < size; j++) { // Iterate columns
      lp = prev + binom[j] + lz;
      buffer[(i * size) + j] = (lp > MIN_LOG) ? lp : -INFINITY;
      prev = lp;
    } // end columns

    vdExp(size, (double *)&(buffer[i * size]), (double *)&(buffer[i * size]));
    cblas_dscal(size, -1.0, (double *)&(buffer[i * size]), 1);

    buffer[(i * size) + row_offset + i] += 1.0;

    for (j = 0; j < size; j++) {
      buffer[(i * size) + j] = (((fabs(buffer[(i * size) + j]) >= threshold) ||
                                 (j == row_offset + i))
                                    ? buffer[(i * size) + j]
                                    : 0);
    }
#pragma omp parallel for reduction(+ : local_nnz)
    for (j = 0; j < size; j++) {
      local_nnz += ((buffer[(i * size) + j] == 0) ? 0 : 1);
    }
    global_nnz[i + row_offset] = local_nnz;
  }

  DKL_INT nnz = 0;
  for (i = row_offset; i < row_offset + block_size; i++) {
    nnz += global_nnz[i];
  }

  // Re-allocate the storage for Q matrix
  *old_size = *current_size;
  *current_size += nnz;
  Q->cols = dkl_realloc(Q->cols, *current_size, DKL_INT);
  assert(Q->cols);
  Q->data = dkl_realloc(Q->data, *current_size, double);
  assert(Q->data);

  // Copy the buffer in the CSR representation into Q
  DKL_INT info;
  DKL_INT *job = dkl_alloc(6, DKL_INT);
  job[0] = 0;
  job[1] = 0;
  job[2] = 0;
  job[3] = 2;
  job[4] = nnz;
  job[5] = 1;

  mkl_ddnscsr(job, &block_size, &size, buffer, &size, &(Q->data[*old_size]),
              &(Q->cols[*old_size]), Q->row_index, &info);
  assert(info == 0);
  dkl_dealloc(job);
}
/**
 * wf_matrix_csr: build the Wright-Fisher matrix
 *
 * @param[in] N Population size (range: 2-Inf)
 * @param[in] s Selection coefficient (range: -1:Inf)
 * @param[in] u Forward mutation rate (range: 0:1/2N)
 * @param[in] v Backward mutation rate (range: 0:1/2N)
 * @param[in] h Dominance coefficient (range: 0:1)
 * @param[in] block_size Number of rows to process per cycle
 * @param[in] threshold Numeric zero threshold
 *
 * @return The CSR of the Wright-Fisher matrix
 */
csr_sparse_matrix *wf_matrix_csr(wf_parameters *wf, DKL_INT block_size,
                                 double threshold) {
  double allele_number = (double)(2 * wf->population_size);
  DKL_INT size = (2 * wf->population_size) - 1;
  DKL_INT i;
  DKL_INT current_size = 0;
  DKL_INT old_size = 0;

  csr_sparse_matrix *Q = dkl_alloc(1, csr_sparse_matrix);

  Q->nrows = size;
  Q->ncols = size;
  Q->nnz = 0;

  Q->data = dkl_new(double);
  Q->cols = dkl_new(DKL_INT);
  Q->row_index = dkl_alloc(Q->nrows + 1, DKL_INT);
  DKL_INT *global_nnz = dkl_alloc(Q->nrows, DKL_INT);

  double *buffer = dkl_alloc(block_size * size, double);
  double *binom = dkl_alloc(size, double);

#pragma omp parallel for
  for (i = 0; i < size; i++) {
    binom[i] = log(allele_number - i) - log(i + 1);
  }

  DKL_INT row_offset;
  DKL_INT rem = size % block_size;

  for (row_offset = 0; row_offset < (size - rem); row_offset += block_size) {
    wf_make_block(wf, threshold, row_offset, block_size, buffer, binom, Q,
                  &current_size, global_nnz, &old_size);
  }
  wf_make_block(wf, threshold, row_offset, rem, buffer, binom, Q, &current_size,
                global_nnz, &old_size);

  DKL_INT csum = 0;
  for (i = 0; i < size; i++) {
    Q->row_index[i] = csum;
    csum += global_nnz[i];
  }

  Q->nnz = csum;
  // assert(csum == total_nnz);
  // Last realloc to minimize the final memory usage
  Q->cols = dkl_realloc(Q->cols, Q->nnz, DKL_INT);
  Q->data = dkl_realloc(Q->data, Q->nnz, double);
  dkl_dealloc(buffer);
  dkl_dealloc(binom);
  dkl_dealloc(global_nnz);
  Q->row_index[Q->nrows] = Q->nnz;
  return (Q);
}
