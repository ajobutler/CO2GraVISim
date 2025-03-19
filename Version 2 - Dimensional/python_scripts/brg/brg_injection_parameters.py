# brg_injection_parameters.py
# This python script stores the injection parameters for use when generating
# a series of batch runs for CO2GraVISim

import numpy as np

# n_inj_locs = 2
# inj_loc_idxs = [[120, 175], [80, 175]]

inj_loc_options = {
    # 0 : [[120, 175], [80, 175]] ,
    # 1 : [[110, 175], [80, 175]]
    # 0 : [[31,43]] ,
    0 : [[31,60]] ,
}

### !!!! I need to modify the below to deal with generating flux values for multiple wells !!!

N_loc = len(inj_loc_options)
# print(f"{N_loc = }")
# print(f"{len(inj_loc_options[0]) = }")


base_injection_time = 1096.0 #365.0  # [days]
base_flux = 8e-3 # 0.0 # 0.002  # [m^3 s^-1]

# n_flux_vals = 2
# flux_times = [ 0 , 1 ]
# flux_vals = [ 1 , 0 ]

if (np.abs(base_flux) <= 1e-10):
    N_flux_variations = 1
    flux_times = [[ 0.0, 365.0 ]]
    flux_vals  = [[ 0.0,   0.0 ]]

else:
    ## Fixed total volume #########################################################################
    N_flux_variations = 30 #1 #11
    # inj_flux_vals = base_flux * np.logspace(-1,1,N_variations,endpoint=True)
    inj_flux_vals = base_flux * np.logspace(
        -1.0, 1.0, N_flux_variations, endpoint=True
        # np.log10(2.0 / 5.0), np.log10(5.0 / 2.0), N_flux_variations, endpoint=True
    )

    fixed_total_volume = base_injection_time * base_flux
    t_stop_vals = fixed_total_volume / inj_flux_vals

    flux_times = []
    flux_vals = []

    for k in range(0, N_flux_variations):
        flux_times.append([0.0, t_stop_vals[k]])
        flux_vals.append([inj_flux_vals[k], 0.0])

N_flux = N_flux_variations


N_injection = N_loc * N_flux

injection_variations = { 'N_variations' : N_injection}

for n_loc, (k_loc,v_loc) in enumerate(inj_loc_options.items()):
    for n_flux in range(0,N_flux):
        n = n_loc*N_flux + n_flux
        injection_variations[n] = {
            "locations" : v_loc ,
            "Flux_times" : flux_times[n_flux] ,
            "Flux_vals" : [flux_vals[n_flux]] ,
        } 



# # # Injection_dict = {
# # #     "z_idx_ceil": 1,  # z index of the upper sealing layer (caprock) of the storage interval
# # #     "z_idx_base": 20,  # z index of the lower sealing layer (basement) of the storage interval
# # #     "locations": [
# # #         [120, 175],
# # #         [80, 175],
# # #     ],  # Array of (i,j) pairs for each of the injection wells.
# # #     # NOTE: This has to be an array rather than a list!
# # #     "Flux_times": [
# # #         0.0,
# # #         365.0,
# # #         730.0,
# # #     ],  # List of times [days] at which the injection flux changes
# # #     "Flux_vals": [
# # #         [0.002, 0.00, 0.00],
# # #         [0.00, 0.001, 0.0],
# # #     ],  # Array of the correspondig flux values [m^3 s^-1] at each of these flux times for each of the injection locations.
# # #     # Each row corresponds to an injection location.
# # #     # NOTE: This has to be an array rather than a list!
# # # }
