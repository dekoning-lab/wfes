#define _POSIX_C_SOURCE 200809L
#include "dkl.h"

void print_help() {
  printf("WFES: Wright-Fisher model solver\n"
         "USAGE:\n"
         " -N, --population_size:        Population size\n"
         " -s, --selection_coefficient:  Selection coefficient\n"
         " -u, --forward_mutation_rate:  Mutation rate from a to A\n"
         " -v, --backward_mutation_rate: Mutation rate from A to a\n"
         " -d, --dominance_coefficient:  Proportion of selection Aa recieves\n"
         "[-m, --selection_mode]:        Selection mode (1: viability; 2: haploid)\n"
         "[-x, --observed_allele_count]: Observed count in the population (for allele age)\n"
         "[-z, --zero_threshold]:        Any number below this is considered "
         "0. Default 1e-25\n"
         "[-g, --generations_file]:      Generations spent with a given number "
         "of copies\n"
         "[-e, --extinction_file]:       Probability of extinction, given the "
         "starting number of copies\n"
         "[-f, --fixation_file]:         Probability of fixation, given the "
         "starting number of copies\n"
         "[--force]:                     Do not preform any parameter validity "
         "checks\n"
         "[--help]:                      Print this message and exit\n");
}

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
void wfes(wf_parameters *wf, wf_statistics *r, double zero_threshold) {
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
  double *_M3 = dkl_alloc(matrix_size, double);

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

/**
* main
*/
int main(int argc, char **argv) {
  // Parse the command line
  bool help = dkl_args_parse_flag(argc, argv, false, "h", "-h", "--h", "help",
                                  "-help", "--help", NULL);
  if (help) {
    print_help();
    exit(DKL_HELP_EXIT);
  }

  // Initialize WF input parameters
  wf_parameters *wf = dkl_new(wf_parameters);

  // Parse the required arguments
  wf->population_size = dkl_args_parse_int(
      argc, argv, true, "n", "-n", "--n", "N", "-N", "--N", "population_size",
      "-population_size", "--population_size", NULL);
  wf->selection = dkl_args_parse_double(
      argc, argv, true, "s", "-s", "--s", "selection_coefficient",
      "-selection_coefficient", "--selection_coefficient", NULL);
  wf->forward_mutation_rate = dkl_args_parse_double(
      argc, argv, true, "u", "-u", "--u", "forward_mutation_rate",
      "-forward_mutation_rate", "--forward_mutation_rate", NULL);
  wf->backward_mutation_rate = dkl_args_parse_double(
      argc, argv, true, "v", "-v", "--v", "backward_mutation_rate",
      "-backward_mutation_rate", "--backward_mutation_rate", NULL);
  wf->dominance_coefficient = dkl_args_parse_double(
      argc, argv, true, "d", "-d", "--d", "dominance_coefficient",
      "-dominance_coefficient", "--dominance_coefficient", NULL);

  // Parse optional parameters
  double zero_threshold = dkl_args_parse_double(
      argc, argv, false, "z", "-z", "--z", "zero_threshold", "-zero_threshold",
      "--zero_threshold", NULL);
  if (dkl_errno == DKL_OPTION_NOT_FOUND) {
#ifdef DEBUG
    println("Using default threshold: 1e-30");
#endif
    zero_threshold = 1e-30;
    dkl_clear_errno();
  }
  
  wf->selection_mode = dkl_args_parse_int(
      argc, argv, false, "m", "-m", "--m", "M", "-M", "--M", "selection_mode",
      "-selection_mode", "--selection_mode", NULL);
  if (dkl_errno == DKL_OPTION_NOT_FOUND) {
#ifdef DEBUG
    println("Using default selection mode: fecundity");
#endif
    wf->selection_mode = 0;
    dkl_clear_errno();
  }
  
  wf->observed_allele_count = dkl_args_parse_int(
      argc, argv, false, "x", "-x", "--x", "X", "-X", "--X", "observed_allele_count",
      "-observed_allele_count", "--observed_allele_count", NULL);
  if (dkl_errno == DKL_OPTION_NOT_FOUND) {
#ifdef DEBUG
    println("Defaulting to not calculating allele age");
#endif
    wf->observed_allele_count = -1;
    dkl_clear_errno();
  }

  // Haploid model check
  if (wf->selection_mode == 2 ) {
	// This is needed because we multiply N by 2 throughout for diploid population
	// (We fix this here and account for it in the WF calculation and the final output)
	wf->population_size /= 2.0;
  }

  char *generations_file = dkl_args_parse_string(
      argc, argv, false, "g", "-g", "--g", "generations_file",
      "-generations_file", "--generations_file", "sojourn_time_file",
      "-sojourn_time_file", "--sojourn_time_file", NULL);
  char *extinction_file = dkl_args_parse_string(
      argc, argv, false, "e", "-e", "--e", "extinction_file",
      "-extinction_file", "--extinction_file", NULL);
  char *fixation_file = dkl_args_parse_string(
      argc, argv, false, "f", "-f", "--f", "fixation_file", "-fixation_file",
      "--fixation_file", NULL);

  bool force = dkl_args_parse_flag(argc, argv, false, "force", "-force",
                                   "--force", NULL);
  if (!force) {
    if (wf->population_size > 500000) {
      error_print("The population_size parameter is too large - the "
                  "computation might take a very long time");
      println("Use `--force` to override");
      exit(DKL_PARAM_ERROR);
    }
    double max_mutation_rate = 1.0 / (2.0 * wf->population_size);
    if (wf->forward_mutation_rate > max_mutation_rate ||
        wf->backward_mutation_rate > max_mutation_rate) {
      error_print(
          "The mutation rate might violate the Wright-Fisher assumptions");
      println("Use `--force` to override");
      exit(DKL_PARAM_ERROR);
    }
    if (zero_threshold > 1e-10) {
      error_print("The zero threshold is too high - this will produce "
                  "inaccurate results");
      println("Use `--force` to override");
      exit(DKL_PARAM_ERROR);
    }
  }

  DKL_INT matrix_size = (2 * wf->population_size) - 1;

  wf_statistics *results = wf_statistics_new(wf->population_size);

  wfes(wf, results, zero_threshold);

  double gensubRate = 1.0 / ( ( 1.0/(2*wf->population_size * wf->forward_mutation_rate) + results->time_extinction ) * ( 1.0/results->probability_fixation - 1) + ( 1.0/(2*wf->population_size * wf->forward_mutation_rate) ) + results->time_fixation );

  // Correct for halpoid size if necessary
  if (wf->selection_mode == 2) {
	wf->population_size *= 2.0;
  }
  // Output the results

  printf("%" PRId64 ",%g,%g,%g,%g,%g,%g,%g,%g,%g,%g", wf->population_size,
         wf->selection, wf->forward_mutation_rate, wf->backward_mutation_rate,
         wf->dominance_coefficient, results->probability_extinction,
         results->probability_fixation, results->time_extinction,
         results->time_fixation, results->count_before_extinction, gensubRate);

  if ( wf->observed_allele_count > 0) {
    printf(",%g", results->expectedAge);
  }
  printf("\n");

  if (generations_file) {
    FILE *f = fopen(generations_file, "w");
    if (f != NULL) {
      for (DKL_INT i = 0; i < matrix_size - 1; i++) {
        fprintf(f, "%g,", results->generations[i]);
      }
      fprintf(f, "%g\n", results->generations[matrix_size - 1]);
      fclose(f);
    }
  }
  if (extinction_file) {
    FILE *f = fopen(extinction_file, "w");
    if (f != NULL) {
      for (DKL_INT i = 0; i < matrix_size - 1; i++) {
        fprintf(f, "%g,", results->extinction_probabilities[i]);
      }
      fprintf(f, "%g\n", results->extinction_probabilities[matrix_size - 1]);
      fclose(f);
    }
  }
  if (fixation_file) {
    FILE *f = fopen(fixation_file, "w");
    if (f != NULL) {
      for (DKL_INT i = 0; i < matrix_size - 1; i++) {
        fprintf(f, "%g,", results->fixation_probabilities[i]);
      }
      fprintf(f, "%g\n", results->fixation_probabilities[matrix_size - 1]);
      fclose(f);
    }
  }

  // Deallocation
  dkl_del(wf);
  wf_statistics_del(results);

  return 0;
}
