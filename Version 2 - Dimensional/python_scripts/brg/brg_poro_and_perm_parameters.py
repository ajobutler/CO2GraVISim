# brg_poro_and_perm_parameters.py
# This python script stores the porosity and permeability parameters for use when generating
# a series of batch runs for CO2GraVISim


from itertools import product

# Parameters for generating porosity and permeability distributions
poro_and_perm_options = {
    "ptype_h": ["Uniform"],
    "pparams_h": [[0.75, 0.3, 5, 5]],
    # "pparams_h": [[0.75, 0.3, 5, 5],[0.75, 0.3, 10, 10],[0.75, 0.3, 15, 15]],
    "ptype_v": ["Uniform"],
    "pparams_v": [[0.6, 5.0, 1.2]],
}

# Extract values from the dictionary
param_values = list(poro_and_perm_options.values())

# Generate a list of all of the possible parameter sets
poro_and_perm_params = list(product(*param_values))

# Number of combinations generated
N_poro_and_perm = len(poro_and_perm_params)

poro_and_perm_variations = { 'N_variations' : N_poro_and_perm }

key_list = [
    "ptype_h",
    "pparams_h",
    "ptype_v",
    "pparams_v",
]

for j in range(0, N_poro_and_perm):
    p = poro_and_perm_params[j]
    D = {}
    for i,k in enumerate(key_list):
        D[k] = p[i]
    poro_and_perm_variations[j] = D
