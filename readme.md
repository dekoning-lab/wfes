# The Wright Fisher Exact Solver

The Wright Fisher Exact Solver ,`WFES` ['double-u fez'] is a toolbox for making fast, scalable computations in population genetics. Given parameters for a Wright Fisher model, `WFES` solves for the probability of fixation or extinction of an allele, the mean time to absorption, and the conditional mean time to fixation or extinction.

## Building
The default target is the shared executable and the `python`-compatible shared library
```
https://github.com/dekoning-lab/wfes
cd wfes
make
```

## Usage

Short command line options:
```
wfes N 1000 s 0.001 u 1e-8 v 1e-8 h 0.5
```

Full command line options:
```
wfes --population_size 1000
     --selection 0.001
     --forward_mutation_rate 1e-8
     --backward_mutation_rate 1e-8
     --dominance 0.5
```

## Output

`WFES` echoes back the model parameters it was invoked with, followed by comma-separated calculated statistics:

|Field|
|---|
|`population_size`|
|`selection`|
|`forward_mutation_rate`|
|`backward_mutation_rate`|
|`dominance`|
|`probability_of_extinction`|
|`probability_of_fixation`|
|`time_to_extinction`|
|`time_to_fixation`|
|`total_count_before_extinction`|

For example:
```
wfes N 1000 s 0.001 u 1e-8 v 1e-8 h 0.5
1000,0.001,1e-8,1e-8,0.5,0.995,0.005,10.0038,396.511,200
```

## Summary

​Given the effective population size, the selection coefficient, the forward and backward mutation rates (into and out of the mutant allele, respectively), and the dominance coefficient, WFES first builds the appropriate Wright-Fisher probability transition matrix. Since this matrix will have a lot of near-zero values in it, we may ask it to truncate these to zero for a reduction in memory and runtime, but a loss in accuracy (see Table 1 of the manuscript). Given a truncation threshold, WFES will consider all matrix entries under the threshold to be zero. We recommend a maximum truncation value of 1e-10.

The selection coefficient "s" and dominance coefficient "h" are related to the fitnesses of the alleles as follows.

Genotype | Fitness
---- | ----
AA | 1+s
Aa | 1+sh
aa | 1

Once the Wright Fisher transition matrix has been computed, WFES will compute the 1) first row of the corresponding fundamental matrix and 2) the entire absorption probability matrix. Using these, properties of the Markov Chain are computed (for further details, see Online Methods of the manuscript). These properties are all dependent on the assumption that the allele entered the population as a single copy. They are:

- The probability of extinction of the allele
- The probability of fixation of the allele
- The mean time to extinction, given that the allele eventually goes to extinction
- The mean time to fixation, given that the allele eventually goes to fixation
- The mean total count before extinction, given that the allele eventually goes to extinction (the sum of all copies of the allele, over all generations before extinction)

The first row of the fundamental matrix is also written to a user-specified file. The i-th entry represents the mean amount of steps that the Markov Chain will spend in state i before absorption (see Online Methods of the manuscript for details).

WFES includes a check for input parameters that are either 1) not within acceptable parameter ranges (see table below) or 2) technically within parameter ranges but unreasonable for the specific model (for example, for a population size of 1 million the computer will likely not have the needed memory requirements). For the latter case, WFES also includes a force option to override these checks. In some specific cases when the selection coefficient is too low, the mean time to fixation given that the allele eventually fixes will be reported as NA, but this is unavoidable.


​#### Valid Input Parameter Values

Parameter | Lower Bound | Upper Bound |
---- | ---- | ----
Population Size | 2 | Inf
Selection Coefficient | -1 | Inf
Forward Mutation Rate | 0 | Inf
Backward Mutation Rate | 0 | Inf
Dominance Coefficient | 0 | 1
Threshold Value | 0 | 1e-10

Note that the first four parameters will give an error if they are too high, even though they have a theoretically infinite upper bound. These errors can be overrode with the force option.
​
