from __future__ import print_function
from numpy import linspace
population_sizes = [50000]
selection = linspace(-1, 1, 101)
mutation_rates = [0, 1e-20, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5]
dominance = linspace(0, 1, 11)
truncation = [1e-20]

for N in population_sizes:
	for s in selection:
		for m in mutation_rates:
			for d in dominance:
				for t in truncation:
					print("-n {0} -s {1} -v {2} -u {3} -h {4} -z {5}".format(N,s,m,m,d,t))
