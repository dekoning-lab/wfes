import numpy as np
from collections import namedtuple
from libc.stdint cimport int64_t
ctypedef int64_t dkl_int

# Define python interface
WF_Statistics = namedtuple("WF_Statistics", ["probability_extinction",
                                             "probability_fixation",
                                             "time_extinction",
                                             "time_fixation",
                                             "count_before_extinction",
                                             "phylogenetic_substitution_rate",
                                             "extinction_probabilities",
                                             "fixation_probabilities",
                                             "generations",
                                             "expected_age"])

# Define C interface
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
        double phylogenetic_substitution_rate

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


def solve(population_size, selection_coefficient, forward_mutation_rate, backward_mutation_rate, dominance_coefficient, zero_threshold=1e-25, observed_allele_count = 0):
    cdef:
        dkl_int matrix_size = (2 * population_size) - 1
        wf_parameters *wf
        wf_statistics *results
        double[::1] extinction_probabilities
        double[::1] fixation_probabilities
        double[::1] generations

    results = wf_statistics_new(population_size)
    wf = wf_parameters_new()

    wf.population_size = population_size
    wf.selection = selection_coefficient
    wf.forward_mutation_rate = forward_mutation_rate
    wf.backward_mutation_rate = backward_mutation_rate
    wf.dominance_coefficient = dominance_coefficient
    wf.observed_allele_count = observed_allele_count

    wf_solve(wf, results, zero_threshold)

    extinction_probabilities = (<double[:matrix_size:1]> results.extinction_probabilities).copy()
    fixation_probabilities = (<double[:matrix_size:1]> results.fixation_probabilities).copy()
    generations = (<double[:matrix_size:1]> results.generations).copy()

    wf_parameters_del(wf)

    #Should this be deferred?
    wf_statistics_del(results)

    return WF_Statistics(results.probability_extinction,
                         results.probability_fixation,
                         results.time_extinction,
                         results.time_fixation,
                         results.count_before_extinction,
                         results.phylogenetic_substitution_rate,
                         np.asarray(extinction_probabilities),
                         np.asarray(fixation_probabilities),
                         np.asarray(generations),
                         results.expected_age)
