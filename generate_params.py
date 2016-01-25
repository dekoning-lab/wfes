from __future__ import print_function
from numpy import linspace
population_sizes = [500000]
selection = linspace(-1, 1, 2001)
mutation_rates = [0]
dominance = linspace(0, 1, 11)
truncation = [1e-30]

for N in population_sizes:
	for s in selection:
		for m in mutation_rates:
			for d in dominance:
				for t in truncation:
					print("-N {0} -s {1} -u {2} -v {3} -d {4} -z {5}".format(N,s,m,m,d,t))
