#include "dkl_wf.h"

wf_statistics *wf_statistics_new(DKL_INT population_size) {
  wf_statistics *r = dkl_new(wf_statistics);

  r->probability_extinction = 0;
  r->probability_fixation = 0;
  r->time_extinction = 0;
  r->time_fixation = 0;
  r->count_before_extinction = 0;
  r->phylogenetic_substitution_rate = 0;

  DKL_INT size = (2 * population_size) - 1;
  r->extinction_probabilities = dkl_alloc(size, double);
  r->fixation_probabilities = dkl_alloc(size, double);
  r->generations = dkl_alloc(size, double);

  r->sojourn_conditional_extinction = dkl_alloc(size, double);
  r->sojourn_conditional_fixation = dkl_alloc(size, double);

  r->expected_age = 0;
  r->expected_age_stdev = 0;

  return r;
}

void wf_statistics_del(wf_statistics *r) {
  dkl_del(r->extinction_probabilities);
  dkl_del(r->fixation_probabilities);
  dkl_del(r->generations);
  dkl_del(r->sojourn_conditional_extinction);
  dkl_del(r->sojourn_conditional_fixation);
  dkl_del(r);
}

wf_parameters *wf_parameters_new(void) {
  wf_parameters *wf = dkl_new(wf_parameters);
  return wf;
}

void wf_parameters_del(wf_parameters *wf) { dkl_del(wf); }

// Ewens - Mathematical population genetics equation 3.29
double wf_sampling_coefficient(wf_parameters *wf, DKL_INT i) {
  double N = wf->population_size;
  double v = wf->forward_mutation_rate;
  double u = wf->backward_mutation_rate;
  double s = wf->selection;
  double h = wf->dominance_coefficient;

  double a = (1 + s) * (i * i);
  double b = (1 + (s * h)) * i * ((2 * N) - i);
  double c = ((2 * N) - i) * ((2 * N) - i);
  double q = (a + b) / (a + (2 * b) + c);

  return ((1 - u) * q) + ((1 - q) * v);
}

double wf_sampling_coefficient_viability(wf_parameters *wf, DKL_INT i) {
  double x = i / (2.0 * wf->population_size);
  double v = wf->forward_mutation_rate;
  double u = wf->backward_mutation_rate;
  double h = wf->dominance_coefficient;
  double s = wf->selection;

  double psi = ((-v + (-1 + u + v) * x) *
                (1 + s * (h + v - h * v + (h - 1) * (v + u - 1) * x))) /
               (-1 +
                s * (-v + (-1 + u + v) * x) *
                    (2 * h + v - 2 * h * v + (2 * h - 1) * (v + u - 1) * x));
  return psi;
}

double wf_sampling_coefficient_haploid(wf_parameters *wf, DKL_INT i) {
  double N = wf->population_size * 2.0;
  double v = wf->forward_mutation_rate;
  double u = wf->backward_mutation_rate;
  double s = wf->selection;

  double psi = (i * (s + 1) * (1 - u) + (N - i) * v) / (i * (1 + s) + N - i);
  return psi;
}

void wf_make_block(wf_parameters *wf, double threshold, DKL_INT row_offset,
                   DKL_INT block_size, double *buffer, double *binom,
                   csr_sparse_matrix *Q, DKL_INT *current_size,
                   DKL_INT *global_nnz, DKL_INT *old_size) {
  DKL_INT size = (2 * wf->population_size) - 1;
  DKL_INT i, j;
#pragma omp parallel for private(j)
  for (i = 0; i < block_size; i++) {

    double lq;
    if (wf->selection_mode == 1) {
      lq = wf_sampling_coefficient_viability(wf, row_offset + i + 1);
    } else if (wf->selection_mode == 2) {
      lq = wf_sampling_coefficient_haploid(wf, row_offset + i + 1);
    } else {
      lq = wf_sampling_coefficient(wf, row_offset + i + 1);
    }

    DKL_INT local_nnz = 0;
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

long int factorial(int n) {
  // recursive factorial for small n
  if (n == 0)
    return 1;
  else
    return(n * factorial(n-1));
}

double poisson_prob( int num, double rate ) {
  return ( pow( rate, num )  * exp(-rate)/factorial(num) );
}

int return_max_poisson( double rate, double epsilon ) {
  double current = 0;
  int index = 0;

  do {
    current = poisson_prob( index, rate );
    if ( current >= epsilon ) {
        index++;
    }
  } while ( (current >= epsilon) && ( current <= 1.0 ) );

  index--;
  return index;
}

double return_max_poisson_sum( double rate, double epsilon ) {
  double current = 0;
  double sum = 0.0;
  int index = 1;

  do {
    current = poisson_prob( index, rate );
    if ( current >= epsilon ) {
        sum += current;
        index++;
    }
  } while ( current >= epsilon );

  if (sum == 0) sum = 1.0;
  return sum;
}

void wf_solve(wf_parameters *wf, wf_statistics *r, double zero_threshold) {
  // Declaration
  DKL_INT matrix_size = (2 * wf->population_size) - 1;

  // Pardiso control parameters
  DKL_INT pardiso_matrix_type = MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC;
  DKL_INT pardiso_number_right_hand_sides = 1;
  void *pardiso_internal[MKL_IFS];
  DKL_INT pardiso_control[MKL_IFS];
  DKL_INT pardiso_maximum_factors = 1;
  DKL_INT pardiso_matrix_number = 1;
  DKL_INT pardiso_phase;
  DKL_INT pardiso_error = 0;
  DKL_INT pardiso_message_level = 0;

  // Auxiliary variables
  DKL_INT i;
  double double_dummy;   // double dummy
  DKL_INT integer_dummy; // integer dummy

  // Right hand side vector
  double *y = dkl_alloc(matrix_size, double);
  double *workspace = dkl_alloc(matrix_size, double);

  // If allele age
  double *M_2 = dkl_alloc(matrix_size, double);

  // Allele age variance
  double *M_3 = dkl_alloc(matrix_size, double);
  double *IQ_col = dkl_alloc(matrix_size, double);
  double *Ax = dkl_alloc(matrix_size, double);

  // Block size to calculate at once
  DKL_INT block_size;
  if (matrix_size >= 100) {
    block_size = (DKL_INT)(matrix_size * 0.01); // 1% per block
  } else {
    block_size = matrix_size;
  }

#ifdef DEBUG
  double start_time, end_time, global_start_time, global_end_time;
  MKL_Peak_Mem_Usage(MKL_PEAK_MEM_ENABLE);
  global_start_time = get_current_time();
  start_time = get_current_time();
#endif

  // Build sparse WF matrix
  csr_sparse_matrix *A = wf_matrix_csr(wf, block_size, zero_threshold);

#ifdef DEBUG
  assert(csr_sparse_is_correct(A) == true);
  end_time = get_current_time();
  printf("Building matrix: %gs\n", end_time - start_time);
#endif

  // RHS for solving for column of B
  for (i = 0; i < matrix_size; i++) {
    double q = wf_sampling_coefficient(wf, i + 1);
    y[i] = pow(1 - q, 2 * wf->population_size);
  }

  // Setup Pardiso control parameters
  for (i = 0; i < MKL_IFS; i++) {
    pardiso_control[i] = MKL_PARDISO_DEFAULT;
    pardiso_internal[i] = MKL_PARDISO_DEFAULT;
  }

  pardiso_control[MKL_PARDISO_DEFAULT_SETTINGS] = MKL_PARDISO_FALSE;
  pardiso_control[MKL_PARDISO_FILL_IN_REDUCING_ORDERING_OPTION] =
  MKL_PARDISO_FILL_IN_REDUCING_ORDERING_NESTED_DISSECTION_OMP;
  pardiso_control[MKL_PARDISO_RETURN_OPTION] = MKL_PARDISO_RETURN_OVERRIDE;
  pardiso_control[MKL_PARDISO_ITERATIVE_REFINEMENT_MAX] = 2;
  pardiso_control[MKL_PARDISO_PIVOTING_PERTURBATION] = 20; // Perturb the pivot elements with 1E-20
  pardiso_control[MKL_PARDISO_SCALING_OPTION] = MKL_PARDISO_SCALING_ENABLE;
  pardiso_control[MKL_PARDISO_SOLVE_OPTION] = MKL_PARDISO_DEFAULT;
  pardiso_control[MKL_PARDISO_WEIGHTED_MATCHING_OPTION] =
  MKL_PARDISO_WEIGHTED_MATCHING_ENABLE;
  pardiso_control[MKL_PARDISO_PRECISION_OPTION] = MKL_PARDISO_PRECISION_DOUBLE;
  pardiso_control[MKL_PARDISO_INDEXING_OPTION] = MKL_PARDISO_INDEXING_ZERO;
  pardiso_control[MKL_PARDISO_OOC_OPTION] = MKL_PARDISO_OOC_OVERFLOW;

#ifdef DEBUG
  pardiso_message_level = 1; // Print statistical information
  pardiso_control[MKL_PARDISO_REPORT_NNZ_FACTORS] = MKL_PARDISO_REPORT_ENABLE;
  pardiso_control[MKL_PARDISO_REPORT_FLOP_FACTOR_PHASE] =
      MKL_PARDISO_REPORT_ENABLE;
  pardiso_control[MKL_PARDISO_REPORT_CGS_CG_DIAGNOSTIC] =
      MKL_PARDISO_REPORT_ENABLE;
  pardiso_control[MKL_PARDISO_MATRIX_CHECK_OPTION] =
      MKL_PARDISO_MATRIX_CHECK_ENABLE;
#endif

  //}}}1

  // Symbolic Factorization
  pardiso_phase = MKL_PARDISO_SOLVER_PHASE_ANALYSIS;
#ifdef DEBUG
  start_time = get_current_time();
#endif
  pardiso_64(pardiso_internal, &pardiso_maximum_factors, &pardiso_matrix_number,
             &pardiso_matrix_type, &pardiso_phase, &matrix_size, A->data,
             A->row_index, A->cols, &integer_dummy,
             &pardiso_number_right_hand_sides, pardiso_control,
             &pardiso_message_level, &double_dummy, &double_dummy,
             &pardiso_error);
  if (pardiso_error != 0) {
    printf("ERROR during symbolic factorization: %" PRId64 "\n", pardiso_error);
    exit(pardiso_phase);
  }
#ifdef DEBUG
  end_time = get_current_time();
  printf("Symbolic factorization %gs\n", end_time - start_time);
#endif

  // Numeric factorization
  pardiso_phase = MKL_PARDISO_SOLVER_PHASE_NUMERICAL_FACTORIZATION;
#ifdef DEBUG
  start_time = get_current_time();
#endif
  pardiso_64(pardiso_internal, &pardiso_maximum_factors, &pardiso_matrix_number,
             &pardiso_matrix_type, &pardiso_phase, &matrix_size, A->data,
             A->row_index, A->cols, &integer_dummy,
             &pardiso_number_right_hand_sides, pardiso_control,
             &pardiso_message_level, &double_dummy, &double_dummy,
             &pardiso_error);
  if (pardiso_error != 0) {
    printf("ERROR during numeric factorization: %" PRId64 "\n", pardiso_error);
    exit(pardiso_phase);
  }
#ifdef DEBUG
  end_time = get_current_time();
  printf("Numeric factorization %gs\n", end_time - start_time);
#endif

#ifdef DEBUG
  start_time = get_current_time();
#endif

// [ LU decomposition finished ]

// For integrating over p when requested
int largest_p = return_max_poisson( 2.0 * wf->population_size * wf->forward_mutation_rate, wf->integration_cutoff );
double z = return_max_poisson_sum( 2.0 * wf->population_size * wf->forward_mutation_rate, wf->integration_cutoff );

// Default is to just calculate for requested p
int max_value = wf->initial_count;

// Otherwise use the largest_p that satisfies our criterion
if (wf->integration_cutoff > 0) {
  max_value = largest_p;
  wf->initial_count = 1;
}

// Store the initial count since this will be changing inside the loop
double stored_initial = wf->initial_count;

// Initialize our variables, since we'll be accumulating contributions over values of p
r->expected_age = 0.0;
r->expected_age_stdev = 0.0;
r->time_extinction = 0;
r->time_fixation = 0;
r->count_before_extinction = 0;
r->phylogenetic_substitution_rate = 0;

for (i = 0; i < matrix_size; i++) {
    r->extinction_probabilities[i] = 0.0;   // Pext
    r->fixation_probabilities[i] = 0.0;     // Pfix
    r->generations[i] = 0.0;                // expected number of visits given p - not integrable!
    r->sojourn_conditional_extinction[i] = 0.0;
    r->sojourn_conditional_fixation[i] = 0.0;
}

// Temporary storage for current iteration's quantities
double *fix_probs = dkl_alloc(matrix_size, double);
double *ext_probs = dkl_alloc(matrix_size, double);
double exp_age, exp_age_var, phylo;

// Iterate over all starting values
for (int pp = stored_initial; pp <= max_value; pp++) {
  // The probability used is the poisson probability of pp starting copies
  // normalized by the total probability of between 1 and largest_p copies
  double prob = poisson_prob( pp, 2.0 * wf->population_size * wf->forward_mutation_rate ) / z;
  if ( stored_initial == max_value ) prob = 1.0;

#ifdef DEBUG
  printf("Solving linear systems assuming p=%d with probability %f...\n", pp, prob );
#endif

  // Set our use variable
  wf->initial_count = pp;

  // RHS for solving for column of B
  for (i = 0; i < matrix_size; i++) {
    double q = wf_sampling_coefficient(wf, i + 1);
    y[i] = pow(1 - q, 2 * wf->population_size);
  }

  // Solution
  pardiso_phase = MKL_PARDISO_SOLVER_PHASE_SOLVE_ITERATIVE_REFINEMENT;
  pardiso_control[11] = 0;

  // Solve for the second column of B matrix (absorption probs), given y=R_0
  // (equation 20 and 8; WFES)
  pardiso_64(pardiso_internal, &pardiso_maximum_factors, &pardiso_matrix_number,
             &pardiso_matrix_type, &pardiso_phase, &matrix_size, A->data,
             A->row_index, A->cols, &integer_dummy,
             &pardiso_number_right_hand_sides, pardiso_control,
             &pardiso_message_level, y, workspace, &pardiso_error);

  if (pardiso_error != 0) {
    printf("ERROR during solution: %" PRId64 "\n", pardiso_error);
    exit(pardiso_phase);
  } else {
    memcpy(ext_probs, y, matrix_size * sizeof(double));
    for (i = 0; i < matrix_size; i++) {
      if (ext_probs[i] < 0) {
        ext_probs[i] = 0.0;
      }
      fix_probs[i] = 1.0 - ext_probs[i];
    }
  }

  // Solve for the p'th row of N -> "generations" (equation 19; WFES)

  // Set y = Ip for p=initial_count [Note, first index is 1 not 0]
  for (i = 0; i < wf->initial_count; i++) {
    y[i] = 0.0;
  }
  y[wf->initial_count - 1] = 1.0;
  for (i = wf->initial_count; i < matrix_size; i++) {
    y[i] = 0.0;
  }

  pardiso_control[MKL_PARDISO_SOLVE_OPTION] = MKL_PARDISO_SOLVE_TRANSPOSED;
  pardiso_64(pardiso_internal, &pardiso_maximum_factors, &pardiso_matrix_number,
             &pardiso_matrix_type, &pardiso_phase, &matrix_size, A->data,
             A->row_index, A->cols, &integer_dummy,
             &pardiso_number_right_hand_sides, pardiso_control,
             &pardiso_message_level, y, workspace, &pardiso_error);
  if (pardiso_error != 0) {
    printf("ERROR during solution: %" PRId64 "\n", pardiso_error);
    exit(pardiso_phase);
  } else {
    //
    memcpy(r->generations, y, matrix_size * sizeof(double));
    for (i = 0; i < matrix_size; i++) {
      if (r->generations[i] < 0) {
        r->generations[i] = 0;
      }
    }
  }

  if (wf->observed_allele_count > 0) {
    // Solve for M2 (equation 23; WFES)
    pardiso_control[MKL_PARDISO_SOLVE_OPTION] = MKL_PARDISO_SOLVE_TRANSPOSED;
    pardiso_64(pardiso_internal, &pardiso_maximum_factors,
               &pardiso_matrix_number, &pardiso_matrix_type, &pardiso_phase,
               &matrix_size, A->data, A->row_index, A->cols, &integer_dummy,
               &pardiso_number_right_hand_sides, pardiso_control,
               &pardiso_message_level, y, workspace, &pardiso_error);
    if (pardiso_error != 0) {
      printf("ERROR during solution: %" PRId64 "\n", pardiso_error);
      exit(pardiso_phase);
    } else {
      // Now store M2 as intermediate result (equation 23; WFES)
      memcpy(M_2, y, matrix_size * sizeof(double));
    }

    // Set x-th column of Q
    memset(y, 0.0, matrix_size * sizeof(double));
    memset(IQ_col, 0.0, matrix_size * sizeof(double));

    int j;
    for (i = 1; i < A->nrows; i++) {
      for (j = A->row_index[i - 1]; j < A->row_index[i]; j++) {
        if (A->cols[j] == (wf->observed_allele_count - 1)) {
          y[i - 1] = (-1 * A->data[j]);
          y[i - 1] += (i == wf->observed_allele_count) ? 1 : 0;
        }
      }
    }

    // Copy column x of Q to IQ_col
    memcpy(IQ_col, y, matrix_size * sizeof(double));

    exp_age = 0;
    for (i = 0; i < matrix_size; i++) {
      exp_age += (M_2[i] * y[i]);
    }
    exp_age /= r->generations[wf->observed_allele_count - 1];

    // Solve for M3 for variance
    memcpy(y, M_2, matrix_size * sizeof(double));
    pardiso_control[MKL_PARDISO_SOLVE_OPTION] = MKL_PARDISO_SOLVE_TRANSPOSED;
    pardiso_64(pardiso_internal, &pardiso_maximum_factors,
               &pardiso_matrix_number, &pardiso_matrix_type, &pardiso_phase,
               &matrix_size, A->data, A->row_index, A->cols, &integer_dummy,
               &pardiso_number_right_hand_sides, pardiso_control,
               &pardiso_message_level, y, workspace, &pardiso_error);
    if (pardiso_error != 0) {
      printf("ERROR during solution: %" PRId64 "\n", pardiso_error);
      exit(pardiso_phase);
    } else {
      // Now store M3 as intermediate result (equation 23; WFES)
      memcpy(M_3, y, matrix_size * sizeof(double));
    }

    // Calculate Ax = xth column of Q(I+Q)

    // Set IQ_col to be xth column of I+Q
    IQ_col[wf->observed_allele_count - 1] += 1.0;

    int k;
    // Get column Ax
    for (i = 1; i < A->nrows; i++) {
      Ax[i - 1] = 0.0;
      // in row i, iterate over all non-zero entries, k
      for (k = A->row_index[i - 1]; k < A->row_index[i]; k++) {
        if ((i - 1) == A->cols[k]) {
          // in column A->cols[k]
          Ax[i - 1] += (-1 * A->data[k] + 1.0) * IQ_col[A->cols[k]];
        } else {
          Ax[i - 1] += (-1 * A->data[k]) * IQ_col[A->cols[k]];
        }
      }
    }

    exp_age_var = 0.0;
    double second_moment = 0;
    for (i = 0; i < matrix_size; i++) {
      second_moment += (M_3[i] * Ax[i]);
    }
    second_moment /= r->generations[wf->observed_allele_count - 1];

    exp_age_var = sqrt(second_moment - pow(r->expected_age, 2.0));

  } else {
    exp_age = NAN;
    exp_age_var = NAN;
  }

  double t_ext, t_fix, count;
  t_ext = t_fix = count = 0.0;

  // Calculate the summary statistics
  for (i = 0; i < matrix_size; i++) {
    t_ext += (ext_probs[i] * r->generations[i]);
    t_fix += (fix_probs[i] * r->generations[i]);
    count +=
        (r->generations[i] * ext_probs[i]) * (i + 1);
    r->sojourn_conditional_extinction[i] = ext_probs[i] / ext_probs[wf->initial_count - 1] * r->generations[i];
    r->sojourn_conditional_fixation[i] = fix_probs[i] / fix_probs[wf->initial_count - 1] * r->generations[i];
  }
  t_ext /= ext_probs[wf->initial_count - 1];
  t_fix /= fix_probs[wf->initial_count - 1];

  count /=
      ext_probs[wf->initial_count - 1];

  double Pfix, Pext;
  Pfix = Pext = 0;

  if (ext_probs[wf->initial_count - 1] <= 0) {
    Pext = 0;
    t_ext = NAN;
  } else {
    Pext =
        ext_probs[wf->initial_count - 1];
  }
  // This could happen if p_ext is very close to zero
  if (t_ext <= 0) {
    t_ext = NAN;
  }
  if (fix_probs[wf->initial_count - 1] <= 0) {
    Pfix = 0;
    t_fix = NAN;
  } else {
    Pfix = fix_probs[wf->initial_count - 1];
  }
  // This could happen if p_fix is very close to zero
  if (t_fix <= 0) {
    t_fix = NAN;
  }
  // Kimura's substitution rate given exact Pfix
  phylo = 2.0 * wf->population_size * wf->forward_mutation_rate * Pfix;

  // Accumulate integrals
  r->probability_fixation += (prob * Pfix);
  r->probability_extinction += (prob * Pext);
  r->time_fixation   += (prob * t_fix);
  r->time_extinction += (prob * t_ext);
  r->count_before_extinction += (prob * count);
  r->phylogenetic_substitution_rate += (prob * phylo);
  r->expected_age_stdev += (prob * exp_age_var);
  r->expected_age += (prob * exp_age);

  // Diagnostic output
  // printf("Expected age for %d is %f (prob %f)\n", pp, exp_age, prob);
}

wf->initial_count = stored_initial;

#ifdef DEBUG
  end_time = get_current_time();
  printf("Solution %gs\n", end_time - start_time);
#endif

  // Memory release
  pardiso_phase = MKL_PARDISO_SOLVER_PHASE_RELEASE_MEMORY_ALL;
  pardiso_64(pardiso_internal, &pardiso_maximum_factors, &pardiso_matrix_number,
             &pardiso_matrix_type, &pardiso_phase, &matrix_size, &double_dummy,
             A->row_index, A->cols, &integer_dummy,
             &pardiso_number_right_hand_sides, pardiso_control,
             &pardiso_message_level, &double_dummy, &double_dummy,
             &pardiso_error);

  dkl_dealloc( fix_probs );
  dkl_dealloc( ext_probs );

  dkl_dealloc(y);
  dkl_dealloc(workspace);

  dkl_dealloc(A->data);
  dkl_dealloc(A->cols);
  dkl_dealloc(A->row_index);

#ifdef DEBUG
  printf("Memory used: %.3g GB\n",
         (double)MKL_Peak_Mem_Usage(MKL_PEAK_MEM) / (GB_CONV));
  global_end_time = get_current_time();
  printf("Total runtime %gs\n", global_end_time - global_start_time);
#endif
}
