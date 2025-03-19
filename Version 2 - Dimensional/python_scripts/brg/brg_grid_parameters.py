# brg_grid_parameters.py
# This python script stores the grid parameters for use when generating
# a series of batch runs for CO2GraVISim

import numpy as np
from itertools import product

grid_params_options = {
    "nx": [64],  # Number of grid points in x
    "ny": [118],  # Number of grid points in y
    "nz": [58],  # Number of grid points in z
    "dx": [50.0],  # Grid spacing in x
    "dy": [50.0],  # Grid spacing in y
}


# Extract values from the dictionary
param_values = list(grid_params_options.values())

# Generate a list of all of the possible parameter sets
grid_params = list(product(*param_values))

# Number of combinations generated
N_grid = len(grid_params)


grid_parameter_variations = {"N_variations": N_grid}

key_list = [
    "nx",
    "ny",
    "nz",
    "dx",
    "dy",
]

for j in range(0, N_grid):
    p = grid_params[j]
    D = {}
    for i,k in enumerate(key_list):
        D[k] = p[i]
    grid_parameter_variations[j] = D


### Need to adjust scripts that use this, as it used to be nx, ny, etc. without the dictionary
