#include "dkl_wf.h"


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

wf_parameters *wf_parameters_new(void) {
  wf_parameters *wf = dkl_new(wf_parameters);
  return wf;
}

void wf_parameters_del(wf_parameters *wf) {
  dkl_del(wf);
}

double wf_sampling_coefficient(wf_parameters *wf, DKL_INT i) {
  double a = (1 + wf->selection) * (i * i);
  double b = (1 + (wf->selection * wf->dominance_coefficient)) * i *
             ((2 * wf->population_size) - i);
  double c = ((2 * wf->population_size) - i) * ((2 * wf->population_size) - i);
  double q = (a + b) / (a + (2 * b) + c);
  return ((1 - wf->forward_mutation_rate) * q) +
         ((1 - q) * wf->backward_mutation_rate);
}

double wf_sampling_coefficient_viability(wf_parameters *wf, DKL_INT i) {
  double x = i/(2.0*wf->population_size);
  double v = wf->forward_mutation_rate;
  double u = wf->backward_mutation_rate;
  double h = wf->dominance_coefficient;
  double s = wf->selection;

  double psi = ( ( -v + (-1 + u + v)*x ) * (1 + s * (h + v - h*v + (h-1)*(v+u-1)*x) ) ) / ( -1+s * (-v+(-1+u+v)*x ) * (2*h+v-2*h*v+( 2*h-1 ) * (v+u-1) * x ) );
  return psi;
}

double wf_sampling_coefficient_haploid(wf_parameters *wf, DKL_INT i) {
  double N = wf->population_size * 2.0;
  double v = wf->forward_mutation_rate;
  double u = wf->backward_mutation_rate;
  double s = wf->selection;

  double psi = ( i * (s+1) * (1-u) + (N-i)*v ) / ( i * (1+s) + N-i );
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

void wf_solve(wf_parameters *wf, wf_statistics *r, double zero_threshold) {
  // Declaration
  DKL_INT matrix_size = (2 * wf->population_size) - 1;
  DKL_INT mtype = 11;
  DKL_INT nrhs = 1;
  void *pt[64];

  // Pardiso control parameters
  DKL_INT iparm[64];
  DKL_INT maxfct, mnum, phase, error, msglvl;

  // Auxiliary variables
  DKL_INT i;
  double ddum;  // double dummy
  DKL_INT idum; // integer dummy

  double *rhs = dkl_alloc(matrix_size, double);
  double *workspace = dkl_alloc(matrix_size, double);

  // If allele age
  double *_M2 = dkl_alloc(matrix_size, double);
  // double *_M3 = dkl_alloc(matrix_size, double);

  DKL_INT block_size;
  if (matrix_size >= 100) {
    block_size = matrix_size * 0.01;
  } else {
    block_size = matrix_size;
  }

  double start_time, end_time, global_start_time, global_end_time;

#ifdef DEBUG
  MKL_Peak_Mem_Usage(MKL_PEAK_MEM_ENABLE);
  global_start_time = get_current_time();
  start_time = get_current_time();
#endif

  csr_sparse_matrix *A = wf_matrix_csr(wf, block_size, zero_threshold);

#ifdef DEBUG
  assert(csr_sparse_is_correct(A) == true);
  end_time = get_current_time();
  printf("Building matrix: %gs\n", end_time - start_time);
#endif

  for (i = 0; i < matrix_size; i++) {
        rhs[i]=0;
  }

  // RHS for solving for column of B
  for (i = 0; i < matrix_size; i++) {
    double q = wf_sampling_coefficient(wf, i + 1);
    rhs[i] = pow(1 - q, 2 * wf->population_size);
  }

  // Setup Pardiso control parameters
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }

  iparm[0] = 1; // No solver default
  iparm[1] = 3; // Fill-in reordering from METIS
  iparm[2] = 1;
  iparm[3] = 0;   // No iterative-direct algorithm
  iparm[4] = 0;   // No user fill-in reducing permutation
  iparm[5] = 1;   // Don't write solution into x
  iparm[6] = 0;   // Not in use
  iparm[7] = 2;   // Max numbers of iterative refinement steps
  iparm[8] = 0;   // Not in use
  iparm[9] = 20;  // Perturb the pivot elements with 1E-20
  iparm[10] = 1;  // Use nonsymmetric permutation and scaling MPS
  iparm[11] = 0;  // Conjugate transposed/transpose solve
  iparm[12] = 1;  // Maximum weighted matching algorithm is switched-on
  iparm[13] = 0;  // Output: Number of perturbed pivots
  iparm[14] = 0;  // Not in use
  iparm[15] = 0;  // Not in use
  iparm[16] = 0;  // Not in use
  iparm[17] = -1; // Output: Number of nonzeros in the factor LU
  iparm[18] = -1; // Output: Mflops for LU factorization
  iparm[19] = 0;  // Output: Numbers of CG Iterations
  iparm[26] = 0;  // Double precision
  iparm[34] = 1;
  iparm[59] = 1;

#ifdef DEBUG
  msglvl = 1; // Print statistical information
#else
  msglvl = 0; // Print statistical information in file
#endif

  maxfct = 1; // Maximum number of numerical factorizations.
  mnum = 1;   // Which factorization to use.
  error = 0;  // Initialize error flag

  for (i = 0; i < 64; i++) {
    pt[i] = 0;
  }
  //}}}1

  // Symbolic Factorization
  phase = 11;
#ifdef DEBUG
  start_time = get_current_time();
#endif
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &matrix_size, A->data,
             A->row_index, A->cols, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum,
             &error);
  if (error != 0) {
    printf("ERROR during symbolic factorization: %" PRId64 "\n", error);
    exit(phase);
  }
#ifdef DEBUG
  end_time = get_current_time();
  printf("Symbolic factorization %gs\n", end_time - start_time);
#endif

  // Numeric factorization
  phase = 22;
#ifdef DEBUG
  start_time = get_current_time();
#endif
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &matrix_size, A->data,
             A->row_index, A->cols, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum,
             &error);
  if (error != 0) {
    printf("ERROR during numeric factorization: %" PRId64 "\n", error);
    exit(phase);
  }
#ifdef DEBUG
  end_time = get_current_time();
  printf("Numeric factorization %gs\n", end_time - start_time);
#endif

  // Solution
  phase = 33;

  iparm[11] = 0;
#ifdef DEBUG
  start_time = get_current_time();
#endif

  // Solve for the second column of B matrix (absorption probs), given rhs=R_2N (equation 20 and 8; WFES)
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &matrix_size, A->data,
             A->row_index, A->cols, &idum, &nrhs, iparm, &msglvl, rhs,
             workspace, &error);
  if (error != 0) {
    printf("ERROR during solution: %" PRId64 "\n", error);
    exit(phase);
  } else {
    // NOTE: This seems backwards; should be fixation_probabilities!
    memcpy(r->extinction_probabilities, rhs, matrix_size * sizeof(double));
    for (i = 0; i < matrix_size; i++) {
      if (r->extinction_probabilities[i] < 0) {
        r->extinction_probabilities[i] = 0.0;
      }
      r->fixation_probabilities[i] = 1.0 - r->extinction_probabilities[i];
    }
  }

  // Solve for the first row of N -> "generations" (equation 19; WFES)
  // Set rhs = Ip for p=1 [Note, first index is 1 not 0]
  rhs[0] = 1.0;
  for (i = 1; i < matrix_size; i++) {
    rhs[i] = 0.0;
  }
  iparm[11] = 2;
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &matrix_size, A->data,
             A->row_index, A->cols, &idum, &nrhs, iparm, &msglvl, rhs,
             workspace, &error);
  if (error != 0) {
    printf("ERROR during solution: %" PRId64 "\n", error);
    exit(phase);
  } else {
    //
    memcpy(r->generations, rhs, matrix_size * sizeof(double));
    for (i = 0; i < matrix_size; i++) {
      if (r->generations[i] < 0) {
        r->generations[i] = 0;
      }
    }
  }

  if ( wf->observed_allele_count > 0) {

  for (i = 0; i < matrix_size; i++) {
	rhs[i] = r->generations[i];
  }

  // Solve for M2 (equation 23; WFES)
  iparm[11] = 2;
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &matrix_size, A->data,
             A->row_index, A->cols, &idum, &nrhs, iparm, &msglvl, rhs,
             workspace, &error);
  if (error != 0) {
    printf("ERROR during solution: %" PRId64 "\n", error);
    exit(phase);
  } else {
    // Now store M2 as intermediate result (equation 23; WFES)
    memcpy( _M2, rhs, matrix_size * sizeof(double));
  }

  /*
  // Solve for M3 (equation 26; WFES)
  iparm[11] = 2;
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &matrix_size, A->data,
             A->row_index, A->cols, &idum, &nrhs, iparm, &msglvl, rhs,
             workspace, &error);
  if (error != 0) {
    printf("ERROR during solution: %" PRId64 "\n", error);
    exit(phase);
  } else {
    // Now store M3 as intermediate result (equation 26; WFES)
    memcpy(_M3, rhs, matrix_size * sizeof(double));
    }
  }
  */

  // Set x-th column of Q
  for (i = 0; i < matrix_size; i++) {
	rhs[i]=0;
  }
/*
  int count=0;
  int start = A->row_index[wf->observed_allele_count-1];
  int end = A->row_index[wf->observed_allele_count];
  printf("Check of first rowindex %d and col %d\n", A->row_index[0], A->cols[0]);

  for (i=start; i<end; i++) {
	if (count == (wf->observed_allele_count-1) ) {
		rhs[count++] = (A->data[i] * -1) + 1;
	} else {
		rhs[count++] = (A->data[i] * -1);
	}
	printf("Got A[%d][%d] = %f\n", A->cols[i], wf->observed_allele_count-1, rhs[count-1]);
  }  */

  int j;
  for (i=1; i<A->nrows; i++) {
	for (j=A->row_index[i-1]; j<A->row_index[i]; j++) {
		if ( A->cols[j] == (wf->observed_allele_count-1) ) {
			if ( i == wf->observed_allele_count ) {
				rhs[i-1] = ( A->data[j] * -1) + 1;
			} else {
				rhs[i-1] = A->data[j] * -1;
			}
			continue;
		}
	}
  }
  r->expectedAge = 0;
  for (i=0; i < matrix_size; i++) {
	r->expectedAge += ( _M2[i] * rhs[i] );
  }
  r->expectedAge /= r->generations[ wf->observed_allele_count -1 ];

  // printf("Expected age given %d obs is %f\n", wf->observed_allele_count, r->expectedAge);

  // The variance requires a column of Q(I+Q), which is O((2N-1)^2) to compute. Let's skip this.
  // r->ageVariance = xx / r->generations[ wf->observed_allele_count ]; //

}
#ifdef DEBUG
  end_time = get_current_time();
  printf("Solution %gs\n", end_time - start_time);
#endif


  // Calculate the summary statistics
  for (i = 0; i < matrix_size; i++) {
    r->time_extinction += (r->extinction_probabilities[i] * r->generations[i]);
    r->time_fixation += (r->fixation_probabilities[i] * r->generations[i]);
    r->count_before_extinction +=
        (r->generations[i] * r->extinction_probabilities[i]) * (i + 1);
  }
  r->time_extinction /= r->extinction_probabilities[0];
  r->time_fixation /= r->fixation_probabilities[0];

  r->count_before_extinction /= r->extinction_probabilities[0];

  if (r->extinction_probabilities[0] <= 0) {
    r->probability_extinction = 0;
    r->time_extinction = NAN;
  } else {
    r->probability_extinction = r->extinction_probabilities[0];
  }
  if (r->fixation_probabilities[0] <= 0) {
    r->probability_fixation = 0;
    r->time_fixation = NAN;
  } else {
    r->probability_fixation = r->fixation_probabilities[0];
  }

  // Memory release
  phase = -1; // Release internal memory
  pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &matrix_size, &ddum,
             A->row_index, A->cols, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum,
             &error);

  dkl_dealloc(rhs);
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