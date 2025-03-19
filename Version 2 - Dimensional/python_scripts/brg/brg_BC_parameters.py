# brg_BC_parameters.py
# This python script stores the boundary condition parameters for use when generating
# a series of batch runs for CO2GraVISim

import numpy as np
from itertools import product

# h_BCs = [1, 1, 1, 1]
# P_BCs = [1, 1, 1, 1]

# h: N,E,S,W , P: N,E,S,W
# NOTE: This needs to be an array, not a list!
BC_options = [[1,1,1,1,1,1,1,1]] 
# BC_options = [[1,1,1,1,1,1,1,1] , [1,1,1,2,1,1,1,2]]


# Number of combinations generated
N_BC = len(BC_options)

BC_variations = {"N_variations": N_BC}

for j in range(0, N_BC):
    p = BC_options[j]
    D = {}
    D["current_thickness"] = {"north": p[0], "east": p[1], "south": p[2], "west": p[3]}
    D["ambient_pressure"]  = {"north": p[4], "east": p[5], "south": p[6], "west": p[7]}
    BC_variations[j] = D


# BC_dict = {
#     "current_thickness": {"north": 1, "east": 1, "south": 1, "west": 1},
#     "ambient_pressure": {"north": 1, "east": 1, "south": 1, "west": 1},
# }
