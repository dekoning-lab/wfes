import wfes

mu = 1e-8

wf = wfes.solve(1000, 0.01, mu, mu, 0.5)

print(wf.probability_extinction, wf.probability_fixation, wf.time_extinction, wf.time_fixation)
