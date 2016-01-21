# The Wright Fisher Exact Solver

The Wright Fisher Exact Solver, `WFES` ['double-u fez'] is a toolbox for making fast, scalable computations in population genetics. Given parameters for a Wright Fisher model, `WFES` solves for the probability of fixation or extinction of an allele, the mean time to absorption, and the conditional mean time to fixation or extinction.

## Building
The default target is the shared executable and the `python`-compatible shared library
```
git clone https://github.com/dekoning-lab/wfes
cd wfes
make
```

## Usage

Short command line options:
```lang=bash
wfes N 1000 s 0.001 u 1e-8 v 1e-8 d 0.5
```

Full command line options:
```lang=bash
wfes --population_size 1000
     --selection_coefficient 0.001
     --forward_mutation_rate 1e-8
     --backward_mutation_rate 1e-8
     --dominance 0.5
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

# Method

​Given the effective population size, the selection coefficient, the forward and backward mutation rates, and the dominance coefficient, `WFES` first builds the appropriate Wright-Fisher probability transition matrix. It then solves for exact long-term behaviors of the model using sparse direct linear solver.

The selection coefficient `s` and dominance coefficient `d` are related to the fitnesses of the alleles as follows:

Genotype | Fitness
---- | ----
AA | 1+s
Aa | 1+sd
aa | 1

`WFES` solves for the following statistics of the model:

- The probability of extinction of the allele
- The probability of fixation of the allele
- The mean time to extinction, given that the allele eventually goes to extinction
- The mean time to fixation, given that the allele eventually goes to fixation
- The mean total count before extinction, given that the allele eventually goes to extinction (the sum of all copies of the allele, over all generations before extinction)

​
##Disk offload

`WFES` uses `MKL PARDISO` linear system solver, which has out-of-core capabilities. Please refer to `MKL` [documentation](https://software.intel.com/en-us/articles/how-to-use-ooc-pardiso). The `pardiso_ooc.cfg` contains the relevant configuration.
