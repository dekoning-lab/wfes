from ctypes import cdll, c_double, c_longlong, pointer, POINTER, Structure, byref
import numpy as np

libwfes = cdll.LoadLibrary("./libwfes.so")

class wf_parameters(Structure):
    _fields_ = [("population_size", c_longlong),
                ("selection", c_double),
                ("forward_mutation_rate", c_double),
                ("backward_mutation_rate", c_double),
                ("dominance_coefficient", c_double)]

class wf_statistics(Structure):
    _fields_ = [("probability_extinction", c_double),
                ("probability_fixation", c_double),
                ("time_extinction", c_double),
                ("time_fixation", c_double),
                ("count_before_extinction", c_double),
                ("extinction_probabilities", POINTER(c_double)),
                ("fixation_probabilities", POINTER(c_double)),
                ("N", POINTER(c_double))]

libwfes.wfes.argtypes = [POINTER(wf_parameters), POINTER(wf_statistics), c_double]
libwfes.wfes.restype = None

def solve(population_size, selection_coefficient, forward_mutation_rate, backward_mutation_rate, dominance_coefficient, zero_threshold = 1e-25):
    matrix_size = (2 * population_size) - 1
    wf = wf_parameters(population_size, selection_coefficient, forward_mutation_rate, backward_mutation_rate, dominance_coefficient)

    extinction_probabilities = np.zeros(matrix_size)
    fixation_probabilities = np.zeros(matrix_size)
    n = np.zeros(matrix_size)

    p_ext = c_double()
    p_fix = c_double()
    t_ext = c_double()
    t_fix = c_double()
    c_ext = c_double()

    result = wf_statistics(p_ext, p_fix, t_ext, t_fix, c_ext, np.ctypeslib.as_ctypes(extinction_probabilities), np.ctypeslib.as_ctypes(fixation_probabilities), np.ctypeslib.as_ctypes(n))

    libwfes.wfes(byref(wf), byref(result), zero_threshold)

    return(result)

if __name__ == '__main__':
    r = solve(100, 0, 0, 0, 0.5, 0)
    print(r.probability_extinction, r.probability_fixation, r.time_extinction, r.time_fixation, r.count_before_extinction)
