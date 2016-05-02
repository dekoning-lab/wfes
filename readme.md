# The Wright Fisher Exact Solver

The Wright Fisher Exact Solver, `WFES` ['double-u fez'] is a toolbox for making **fast, scalable matrix computations** in population genetics and molecular evolution **without diffusion theory approximations or simulation**. Given parameters for a Wright Fisher model, `WFES` exactly calculates a variety of transient and long-term behaviours using efficient sparse parallel linear algebra techniques. On most computers, `WFES` supports population sizes up to about `Ne=100,000`. Truncation of very small values in the transition matrix can be performed (`-z`) to ameliorate the large space-complexity of such models (see Kryukov, DeSanctis and de Koning, 2016). An experimental "out of core" option is also available to support even larger population sizes, however, this option is unstable in the current release.

### Supported exact computations

WFES' exact statistics currently include: 
* the probabilities of fixation and extinction of an allele; 
* the mean time to absorption; 
* the conditional mean time to fixation or extinction; 
* the mean sojourn times in each frequency class (optional); 
* the mean number of copies of an allele on its way to extinction (the "window of opportunity" for secondary mutations in our stochastic tunnelling codon models; in prep.);
* the exact expected age of an allele (optional, if an observed frequency is provided, `-x`; DeSanctis and de Koning, 2016);
* the expected rate of phylogenetic substitution accounting for fast recurrent mutation (in prep.). 

### Model variations

Several variations of the general Wright-Fisher model are supported by default that incorporate two-way mutation, selection, and dominance. These include the standard model of fecundity selection (diploid) `-m 0`, an alternative viabiity selection model (diploid) `-m 1`, and a haploid model accounting for mutation and selection `-m 2`.

*Please cite:* **Kryukov I, DeSanctis B, and APJ de Koning (2016). Efficient techniques for direct analysis of discrete-time population genetic models. Submitted. (BioArXiv link to be added.)**

---
## Building
The default target is the shared executable and the `python`-compatible shared library
```
git clone https://github.com/dekoning-lab/wfes
cd wfes
make
```
Currently, only `linux` systems are supported.

## Usage

Short command line options:
```lang=bash
wfes N 1000 s 0.001 u 1e-8 v 1e-8 d 0.5
```

Full command line options:
```lang=bash
$ ./wfes --help
WFES: Wright-Fisher exact solver
USAGE:
 -N, --population_size:        Population size
 -s, --selection_coefficient:  Selection coefficient
 -u, --forward_mutation_rate:  Mutation rate from a to A (per locus per generation)
 -v, --backward_mutation_rate: Mutation rate from A to a (per locus per generation)
 -d, --dominance_coefficient:  Proportion of selection Aa recieves
[-m, --selection_mode]:        Selection mode (0: fecundity, default; 1: viability; 2: haploid)
[-x, --observed_allele_count]: Observed count in the population (for allele age)
[-z, --zero_threshold]:        Any number below this is considered 0. Default 1e-25
[-g, --generations_file]:      Generations spent with a given number of copies
[-e, --extinction_file]:       Probability of extinction, given the starting number of copies
[-f, --fixation_file]:         Probability of fixation, given the starting number of copies
[--force]:                     Do not preform any parameter validity checks
[--help]:                      Print this message and exit
```

Python interface:
```lang=python
import wfes
mu = 1e-8
wfes.solve(population_size = 1000,
    selection_coefficient = 0.001,
    forward_mutation_rate = mu,
    backward_mutation_rate = mu,
    dominance_coefficient = 0.5)
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
|`expected_allele_age`|
|`phylogenetic_substitution_rate`|

For example:
```
wfes N 1000 s 0.001 u 1e-8 v 1e-8 d 0.5
1000,0.001,1e-8,1e-8,0.5,0.995,0.005,10.0038,396.511,200
```

The output format is dictated by the convenience of producing tables.

## Optional parameters

###Sparsity

Since the Wright-Fisher matrices are generally sparse, there is an optional `zero-threshold` parameter, below which all numbers are considered zero. Judicious use of this cutoff can increase the sparsity of system, while still producing an accurate result. The default cutoff is `1e-25`.

###Vector output

`WFES` can write the full vector solutions that it solves for to disk via the executable. The `python` interface returns these as vectors in the result struct.

There are three optional parameters, which denote file names of the output. Note that currently all the output vectors are indexed from `1`.

- `--generations_file/--sojourn_time_file`: output the expected number of generations the population spends with a given number of copies to file.
- `--extinction_file`: output the probability of extinction, conditioned on starting with a given number of copies.
- `--fixation_file`: output the probability of fixation, conditioned on starting with a given number of copies.

##Python interface

There is a convenience wrapper for `python`, which has been tested against version `3.4` and `3.5`. The module assumes that the `libwfes.so` shared library (build by default) is located in the script directory. The module depends on `numpy`.

## Diploid fitness model

​Given the effective population size, the selection coefficient, the forward and backward mutation rates, and the dominance coefficient, `WFES` first builds the appropriate Wright-Fisher probability transition matrix. It then solves for exact long-term behaviors of the model using sparse direct linear solver.

For a wildtype allele `a` and a mutant allele `A`, the selection coefficient `s` and dominance coefficient `h` are related to the fitnesses of the alleles as follows:

Genotype | Fitness
---- | ----
AA | 1+s
Aa | 1+sh
aa | 1

​
##Disk offload

`WFES` uses `MKL PARDISO` linear system solver, which has out-of-core capabilities. Please refer to `MKL` [documentation](https://software.intel.com/en-us/articles/how-to-use-ooc-pardiso). The `pardiso_ooc.cfg` contains the relevant configuration. *Warning:* this feature is currently experimental.
