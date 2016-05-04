#include "dkl.h"

void print_help(void) {
  printf("WFES: Wright-Fisher model solver\n"
         "USAGE:\n"
         " -n, --population_size:        Population size\n"
         " -s, --selection_coefficient:  Selection coefficient\n"
         " -v, --forward_mutation_rate:  Mutation rate from a to A\n"
         " -u, --backward_mutation_rate: Mutation rate from A to a\n"
         " -h, --dominance_coefficient:  Proportion of selection Aa recieves\n"
         "[-m, --selection_mode]:        Selection mode (1: viability; 2: "
         "haploid)\n"
         "[-x, --observed_allele_count]: Observed count in the population (for "
         "allele age)\n"
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
* main
*/
int main(int argc, char **argv) {
  // Parse the command line
  bool help = dkl_args_parse_flag(argc, argv, false, "--help", NULL);
  if (help) {
    print_help();
    exit(DKL_HELP_EXIT);
  }

  // Initialize WF input parameters
  wf_parameters *wf = wf_parameters_new();

  // Parse the required arguments
  wf->population_size =
      dkl_args_parse_int(argc, argv, true, "-n", "--population_size", NULL);
  wf->selection =
      dkl_args_parse_double(argc, argv, true, "-s", "--selection", NULL);
  wf->forward_mutation_rate =
      dkl_args_parse_double(argc, argv, true, "-u", "--forward_mutation", NULL);
  wf->backward_mutation_rate = dkl_args_parse_double(
      argc, argv, true, "-v", "--backward_mutation", NULL);
  wf->dominance_coefficient =
      dkl_args_parse_double(argc, argv, true, "-h", "--dominance", NULL);

  // Parse optional parameters
  double zero_threshold =
      dkl_args_parse_double(argc, argv, false, "-z", "--zero_threshold", NULL);
  if (dkl_errno == DKL_OPTION_NOT_FOUND) {
#ifdef DEBUG
    println("Using default threshold: 1e-30");
#endif
    zero_threshold = 1e-30;
    dkl_clear_errno();
  }

  wf->selection_mode =
      dkl_args_parse_int(argc, argv, false, "-m", "--selection_mode", NULL);
  if (dkl_errno == DKL_OPTION_NOT_FOUND) {
#ifdef DEBUG
    println("Using default selection mode: fecundity");
#endif
    wf->selection_mode = 0;
    dkl_clear_errno();
  }

  wf->observed_allele_count = dkl_args_parse_int(
      argc, argv, false, "-x", "--observed_allele_count", NULL);
  if (dkl_errno == DKL_OPTION_NOT_FOUND) {
#ifdef DEBUG
    println("Defaulting to not calculating allele age");
#endif
    wf->observed_allele_count = NAN;
    dkl_clear_errno();
  }

  // Haploid model check
  if (wf->selection_mode == 2) {
    // This is needed because we multiply N by 2 throughout for diploid
    // population
    // (We fix this here and account for it in the WF calculation and the final
    // output)
    wf->population_size /= 2.0;
  }

  char *generations_file =
      dkl_args_parse_string(argc, argv, false, "-g", "--generations_file",
                            "--sojourn_time_file", NULL);
  char *extinction_file =
      dkl_args_parse_string(argc, argv, false, "-e", "--extinction_file", NULL);
  char *fixation_file =
      dkl_args_parse_string(argc, argv, false, "-f", "--fixation_file", NULL);

  bool force = dkl_args_parse_flag(argc, argv, false, "--force", NULL);
  if (!force) {
    if (wf->population_size > 500000) {
      error_print("The population_size parameter is too large - the "
                  "computation might take a very long time");
      error_print("Use `--force` to override");
      exit(DKL_PARAM_ERROR);
    }
    double max_mutation_rate = 1.0 / (2.0 * wf->population_size);
    if (wf->forward_mutation_rate > max_mutation_rate ||
        wf->backward_mutation_rate > max_mutation_rate) {
      error_print(
          "The mutation rate might violate the Wright-Fisher assumptions");
      error_print("Use `--force` to override");
      exit(DKL_PARAM_ERROR);
    }
    if (zero_threshold > 1e-10) {
      error_print("The zero threshold is too high - this will produce "
                  "inaccurate results");
      error_print("Use `--force` to override");
      exit(DKL_PARAM_ERROR);
    }
  }

  DKL_INT matrix_size = (2 * wf->population_size) - 1;

  wf_statistics *results = wf_statistics_new(wf->population_size);

  wf_solve(wf, results, zero_threshold);

  // Correct for halpoid size if necessary
  if (wf->selection_mode == 2) {
    wf->population_size *= 2.0;
  }
  // Output the results

  printf("%" PRId64 ",%g,%g,%g,%g,%g,%g,%g,%g,%g,%g", wf->population_size,
         wf->selection, wf->forward_mutation_rate, wf->backward_mutation_rate,
         wf->dominance_coefficient, results->probability_extinction,
         results->probability_fixation, results->time_extinction,
         results->time_fixation, results->count_before_extinction,
         results->phylogenetic_substitution_rate);

  if (wf->observed_allele_count > 0) {
    printf(",%g", results->expected_age);
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
  wf_parameters_del(wf);
  wf_statistics_del(results);

  return 0;
}