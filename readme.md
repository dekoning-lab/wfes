# The Wright Fisher Exact Solver

The Wright Fisher Exact Solver, `WFES` ['double-u fez'] is a toolbox for making **fast, scalable matrix computations** in population genetics and molecular evolution **without diffusion theory approximations or simulation**. Given parameters for a Wright Fisher model, `WFES` exactly calculates a variety of transient and long-term behaviours using efficient sparse parallel linear algebra techniques. On most computers, `WFES` supports population sizes up to about `Ne=100,000`. Truncation of very small values in the transition matrix can be performed (`-z`) to ameliorate the large space-complexity of such models (see Kryukov, DeSanctis and de Koning, 2016). An experimental "out of core" option is also available to support even larger population sizes, however, this option is unstable in the current release.

### Supported exact computations

WFES exact statistics currently include:
* the probabilities of fixation and extinction of an allele;
* the mean time to absorption (extinction or fixation);
* the conditional mean time to fixation or extinction;
* the mean sojourn times in each frequency class;
* the mean number of copies of an allele on its way to extinction (the "window of opportunity" for secondary mutations in our stochastic tunnelling codon models; in prep.);
* the exact expected age of an allele (optional, if an observed frequency is provided, `-x`; DeSanctis and de Koning, 2016).

### Model variations

Several variations of the general Wright-Fisher model are supported by default that incorporate two-way mutation, selection, and dominance. These include the standard model of fecundity selection (diploid) `-m 0`, an alternative viability selection model (diploid) `-m 1`, and a haploid model accounting for mutation and selection `-m 2`.

### Notes

**This software is in beta release.** Please report bugs to the authors by [opening an issue](https://github.com/dekoning-lab/wfes/issues/new). Thanks! - [Ivan](mailto:ikryukov@ucalgary.ca) and [Jason](mailto:jason.dekoning@ucalgary.ca).

*Citation:* **Kryukov I, DeSanctis B, and APJ de Koning (2016). Efficient techniques for direct analysis of discrete-time population genetic models. Submitted. (BioArXiv link to be added.)**

---
## Building

The default target is the shared executable and the `python`-compatible shared library
```
git clone https://github.com/dekoning-lab/wfes
cd wfes
make
```

Both `linux` (with `gcc`) and `OSX` (with `clang`) are supported

### Dependencies

Binary:
* `MKL` libraries

Python:
* `MKL` libraries
* `cython`
* `numpy`

To perform fast numerical analysis, `MKL PARDISO` direct sparse solver is used. The `MKL` libraries are required at compile and run time. These are available directly from [`intel`](https://software.intel.com/en-us/intel-mkl), and through the [`anaconda`](https://www.continuum.io/downloads) `python` distribution. By default, we assume that `anaconda` is installed under `/anaconda/`, and the libraries are compiled against the provided `MKL`.

### Note on OS X

Currently, `clang` does not support `openmp`, so this functionality is restricted to `gcc` at the moment. `OS X` will still use `MKL` threads.

## Usage

Short command line options:
```lang=bash
wfes -n 1000 -s 0.001 -v 1e-8 -u 1e-8 -h 0.5
```

Full command line options:
```lang=bash
$ ./wfes --help
WFES: Wright-Fisher exact solver
USAGE:
 -n, --population_size:        Population size
 -s, --selection_coefficient:  Selection coefficient
 -u, --backward_mutation_rate: Mutation rate from A to a (per locus per generation)
 -v, --forward_mutation_rate:  Mutation rate from a to A (per locus per generation)
 -h, --dominance_coefficient:  Proportion of selection Aa receives
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
|`phylogenetic_substitution_rate`|
|`expected_allele_age`|

For example:
```
./wfes -n 1000 -s 0.01 -u 1e-8 -v 1e-8 -h 0.5 -x 10
1000,0.01,1e-08,1e-08,0.5,0.990068,0.00993225,9.58584,1409.82,212.18,1.98552e-07,56.9915
```

The output format is dictated by the convenience of producing tables.

## Batch runs

We include two scripts to perform array jobs with `slurm` job management system. Use `generate_params.py` to generate a file for the input parameters. Then, submit `array_job_wfes.sh`, which will read the parameter file and launch a separate job for each.

```lang=bash
python generate_params.py > params.txt
```


###Sparsity

Since the Wright-Fisher matrices are generally sparse, there is an optional `zero-threshold` parameter, below which all numbers are considered zero. Judicious use of this cutoff can increase the sparsity of system, while still producing an accurate result. The default cutoff is `1e-25`.

###Vector output

`WFES` can write the full vector solutions that it solves for to disk via the executable. The `python` interface returns these as vectors in the result struct.

There are three optional parameters, which denote file names of the output. Note that currently all the output vectors are indexed from `1`.

- `--generations_file/--sojourn_time_file`: output the expected number of generations the population spends with a given number of copies to file.
- `--extinction_file`: output the probability of extinction, conditioned on starting with a given number of copies.
- `--fixation_file`: output the probability of fixation, conditioned on starting with a given number of copies.

##Python interface

There is a convenience wrapper for `python`. The module uses `cython` to build a shared library module, located in the script directory. The module depends on `numpy`.

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
