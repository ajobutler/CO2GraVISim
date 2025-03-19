# brg_topography_parameters.py
# This python script stores the topography parameters for use when generating
# a series of batch runs for CO2GraVISim

import numpy as np
from itertools import product

# Slope parameters
topo_params_options = {
    "ceil_angle_x": [0.1],
    "ceil_angle_y": [-0.2],
    "base_angle_x": [0.1],
    "base_angle_y": [-0.2],

    "ceil_bump_amp": [0.0], #[2.1], #np.linspace(0,0.1,20), #[0.01, 0.05, 0.1, 0.2],
    # "ceil_bump_amp": np.linspace(0,0.1,20), #[0.01, 0.05, 0.1, 0.2],
    "ceil_lambda_x": [19.0],
    "ceil_lambda_y": [15.0],
    "ceil_decay_x": [50.0],
    "ceil_decay_y": [80.0],

    "base_bump_amp": [0.0],
    "base_lambda_x": [5.0],
    "base_lambda_y": [3.0],
    "base_decay_x": [10.0],
    "base_decay_y": [9.0],
}

# Extract values from the dictionary
param_values = list(topo_params_options.values())

# Generate a list of all of the possible parameter sets
topo_params = list(product(*param_values))

#Number of combinations generated
N_topo = len(topo_params)

topo_variations = { 'N_variations' : N_topo }

for k, p in enumerate(topo_params):
    D = { 'slope' : p[0:4] , 'bump' : p[4:] }
    topo_variations[k] = D

