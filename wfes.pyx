from libc.stdint cimport int64_t
ctypedef int64_t dkl_int


cdef extern from "include/dkl_wf.h":
    struct wf_parameters_t:
        dkl_int population_size
        double selection
        double forward_mutation_rate
        double backward_mutation_rate
        double dominance_coefficient
        int selection_mode
        int observed_allele_count

    struct wf_statistics_t:
        double probability_extinction
        double probability_fixation
        double time_extinction
        double time_fixation
        double count_before_extinction

        double *extinction_probabilities
        double *fixation_probabilities
        double *generations

        double expected_age

    ctypedef wf_parameters_t wf_parameters
    ctypedef wf_statistics_t wf_statistics

    wf_statistics *wf_statistics_new(dkl_int population_size)
    void wf_statistics_del(wf_statistics *r)

    wf_parameters *wf_parameters_new()
    void wf_parameters_del(wf_parameters *wf)

    void wf_solve(wf_parameters *wf, wf_statistics *result, double zero_threshold)


def solve(population_size, selection_coefficient, forward_mutation_rate, backward_mutation_rate, dominance_coefficient, zero_threshold=1e-25):
    cdef:
        wf_parameters *wf
        wf_statistics *results

    results = wf_statistics_new(population_size)
    wf = wf_parameters_new()

    wf.population_size = population_size
    wf.selection = selection_coefficient
    wf.forward_mutation_rate = forward_mutation_rate
    wf.backward_mutation_rate = backward_mutation_rate
    wf.dominance_coefficient = dominance_coefficient

    wf_solve(wf, results, zero_threshold)
    print(results.probability_extinction, results.probability_fixation,
          results.time_extinction, results.time_fixation,
          results.count_before_extinction)
    wf_parameters_del(wf)
    wf_statistics_del(results)
