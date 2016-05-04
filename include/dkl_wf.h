#include "dkl.h"

#ifndef DKL_WF_H
#define DKL_WF_H
/**
 * wf_parameters - parameters of the Wright-Fisher model
 */
typedef struct wf_parameters_t {
  DKL_INT population_size;       // Effective population size
  double selection;              // Selective advantage of 'A' over 'a'
  double forward_mutation_rate;  // Mutation rate from 'A' to 'a'
  double backward_mutation_rate; // Mutation rate from 'a' to 'A'
  double dominance_coefficient;  // Proportion of selective advantage of 'Aa'
                                 // over 'aa'

  int selection_mode;        // 0: fecundity; 1: viability, 2: haploid fecundity
  int observed_allele_count; // (If allele age requested, this is the observed
                             // count)
} wf_parameters;

wf_parameters *wf_parameters_new(void);
void wf_parameters_del(wf_parameters *wf);

/**
 * wf_statistics - summary statistics of the Wright-Fisher model
 */
typedef struct wf_statistics_t {
  double probability_extinction; // Probability of absorption through exinction
  double probability_fixation;   // Probability of absorption through fixation
  double time_extinction; // Expected number of generations until extinction
  double time_fixation;   // Expected number of generations until fixation
  double count_before_extinction; // Total expected count of 'A' in all
                                  // generations before extinction

  double phylogenetic_substitution_rate; // Rate of substitution over long time
                                         // periods
  double *extinction_probabilities;      // Probability of extinction vector,
                                         // conditional on the starting state
  double *fixation_probabilities;        // Probability of fixation vector
  double *generations;                   // Sojourn time vector

  double expected_age; // Expected allele age (DeSanctis and de Koning, 2016)
  // double age_variance; // Expected allele age-variance (DeSanctis and de
  // Koning, 2016)
} wf_statistics;

wf_statistics *wf_statistics_new(DKL_INT population_size);
void wf_statistics_del(wf_statistics *r);

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
double wf_sampling_coefficient(wf_parameters *wf, DKL_INT i);

/**
 * wf_sampling_coefficient_viability: calculate the binomial sampling
 * coefficient (psi)
 * for the Wright-Fisher matrix using viability selection (de Koning)
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
double wf_sampling_coefficient_viability(wf_parameters *wf, DKL_INT i);

/**
 * wf_sampling_coefficient_haploid: calculate the binomial sampling coefficient
 * (psi)
 * for the Wright-Fisher matrix using haploid selection+mutation and fecundity
 * selection (de Koning)
 *
 * @param[in] N Population size (range: 2-Inf)
 * @param[in] s Selection coefficient (range: -1:Inf)
 * @param[in] u Forward mutation rate (range: 0:1/2N)
 * @param[in] v Backward mutation rate (range: 0:1/2N)
 * @param[in] i Current number of copies of 'A' in the population
 *
 * @return The sampling coefficient (psi) for the Wright-Fisher matrix
 */
double wf_sampling_coefficient_haploid(wf_parameters *wf, DKL_INT i);

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
                   DKL_INT *global_nnz, DKL_INT *old_size);

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
                                 double threshold);

/**
 * wfes: solve for conditional times to absorption, sojourn times etc.
 *
 * Input parameters
 * @param[in] N Population size (range: 2-Inf)
 *   The effective population size
 * @param[in] s Selection coefficient (range: -1:Inf)
 *   The relative selective advantage of allele 'A' over allele 'a'
 * @param[in] u Forward mutation rate (range: 0:1/2N)
 *   The rate of mutation from 'A' into 'a'
 * @param[in] v Backward mutation rate (range: 0:1/2N)
 *   The rate of mutation from 'a' into 'A'
 * @param[in] h Dominance coefficient (range: 0:1)
 *   The proportion of selective advantage a heterozygote 'Aa' caries
 * @param[in] t Numeric zero threshold (range: 0:1e-10)
 *   Any number below 't' is considered zero
 * @param[in] m Selection mode (Range: 0,1)
 *   Type of selection (fecundity-default or viability)
 *
 * @param[out] extinction_probabilities Probability of extinction vector (size:
 * 2N-1)
 *   extinction_probabilities[i] is the probability of extinction, given the
 * population starts with i+1 copies of A
 * @param[out] fixation_probabilities Probability of fixation vector (size:
 * 2N-1)
 *   fixation_probabilities[i] is the probability of fixation, given the
 * population starts with i+1 copies of A
 *
 * THIS DESCRIPTION NEEDS TO BE UPDATED
 *
 * @param[out] generations Sojourn time vector (size: 2N-1)
 *  generations[i] is the number of generations the population is expected to
 * spend with i+1 copies of A, given the population starts with 1 copy of A
 *
 */
void wf_solve(wf_parameters *wf, wf_statistics *r, double zero_threshold);
#endif /* DKL_WF_H */
