# brg_initial_profiles.py
# This python script stores the initial profile parameters for use when generating
# a series of batch runs for CO2GraVISim

import numpy as np
from itertools import product

initial_profile_options = {
    'x_c' : [0.0], 
    'y_c' : [0.0], 
    'r_x' : [30.0], 
    'r_y' : [50.0], 
    'h_max' : [0.0],
    # 'h_max' : np.linspace(1.0,5.0,10),
}

# Extract values from the dictionary
param_values = list(initial_profile_options.values())

# Generate a list of all of the possible parameter sets
initial_profile_params = list(product(*param_values))

# Number of combinations generated
N_initial = len(initial_profile_params)


initial_profile_variations = {"N_variations": N_initial}

# key_list = [
#     "x_c",
#     "y_c",
#     "r_x",
#     "r_y",
#     "h_max",
# ]

for j in range(0, N_initial):
    p = initial_profile_params[j]
    D = {}
    for i,k in enumerate(initial_profile_options.keys()):     #(key_list):
        D[k] = p[i]
    initial_profile_variations[j] = D