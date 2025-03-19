# brg_flow_parameters.py
# This python script stores the flow parameters for use when generating
# a series of batch runs for CO2GraVISim

import numpy as np
from itertools import product

flow_params_options = {
    "rho_c": [810.0],  # CO2 density [kg m^-3]
    # "rho_c": np.linspace(630.,810.,15),  # CO2 density [kg m^-3]
    "rho_a_unsat": [1030.0],  # Ambient density [kg m^-3]
    "rho_a_sat": [1042.0],  # Saturated ambient density [kg m^-3]
    "mu_c": [7e-5],  # CO2 density [Pa s]
    "mu_a": [9e-4],  # Ambient density [Pa s]
    # "s_c_r": np.linspace(0.,0.36,20).tolist(),  # Residual CO2 saturation in trapping region [-]
    "s_c_r": [0.36], # Residual CO2 saturation in trapping region [-]
    "s_a_i": [0.2],  # Irreducible ambient saturation in mobile CO2 region [-]
    "krn_mobile": [1.0],  # Relative permeability of the non-wetting phase (CO2) in the mobile region [-]
    "krw_residual": [1.0],  # Relative permeability of the wetting phase (ambient) in the trapping region [-]
    "C_sat": [0.04],  # Volume fraction of CO2 dissolved in ambient [-]
    "g": [9.81],  # gravitational acceleration [m s^-2]
    "D_mol": [2e-9],  # Molecular diffusivity of CO2 in ambient [m^2 s^-1]
    "perm_ratio": [0.1],  # Ratio of vertical to horizontal absolute permeability [-]
    "omega_conv": [0.0],  # Control prefactor for convective dissolution term [-]
}

# Extract values from the dictionary
param_values = list(flow_params_options.values())

# Generate a list of all of the possible parameter sets
flow_params = list(product(*param_values))

# Number of combinations generated
N_flow = len(flow_params)


flow_parameter_variations = {"N_variations": N_flow}

key_list = [
    "rho_c",
    "rho_a_unsat",
    "rho_a_sat",
    "mu_c",
    "mu_a",
    "s_c_r",
    "s_a_i",
    "krn_mobile",
    "krw_residual",
    "C_sat",
    "g",
    "D_mol",
    "perm_ratio",
    "omega_conv",
]

for j in range(0, N_flow):
    p = flow_params[j]
    D = {}
    for i,k in enumerate(key_list):
        D[k] = p[i]
    flow_parameter_variations[j] = D
