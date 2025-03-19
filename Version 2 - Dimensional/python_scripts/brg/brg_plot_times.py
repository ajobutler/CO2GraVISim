# brg_plot_times.py
# This python script stores the plot time parameters for use when generating
# a series of batch runs for CO2GraVISim

## In order to generate the plot times here, I need to add the parent folder for this script
## to the path so that I can import the relevant function.
import sys
import os
# Get the parent directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# Add the parent directory to sys.path
sys.path.append(parent_dir)

from Generate_Plot_times import pt_generate


pt_waypoints = [[0.0, 98615.0]] #[[0.0, 365.0, 1096.0]] #[days]
# pt_waypoints = [[0.0, 365.0, 1096.0], [0.0, 730.0, 1461.0]]
pt_intervals = [[20]] #[[11, 10]]
# pt_intervals = [[11, 10], [6, 20]]

# # pt_waypoints = [[0.0, 5113.0, 10226.0]] #[days]
# # pt_intervals = [[11, 10]]

N_waypoints = len(pt_waypoints)
N_intervals = len(pt_intervals)
N_plot = N_waypoints * N_intervals

plot_times_variations = { 'N_variations' : N_plot }

for k_w in range(0,N_waypoints):
    for k_i in range(0,N_intervals):
        n = k_i*N_waypoints + k_w
        
        plot_times,_,_ = pt_generate(pt_waypoints[k_w], pt_intervals[k_i])
        
        D = { "times" : plot_times }
        plot_times_variations[n] = D


# Plot_times_dict = {
#     "times": [
#         0.00000,
#         365.00000,
#         730.00000,
#         1096.00000,
#         1461.00000,
#         1826.00000,
#         2191.00000,
#         2557.00000,
#         2922.00000,
#         3287.00000,
#         3652.00000,
#         4018.00000,
#         4383.00000,
#         4748.00000,
#         5113.00000,
#     ],
# }
