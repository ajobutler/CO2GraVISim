# batch_probability.py (18/03/24)
# This script takes the output from a batch run of CO2GraVISim
# and plots relevent statistics for plume extent, etc.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec
import warnings

import argparse
from pathlib import Path
import re  # Regular expressions

###################################################################################################
## Functions ######################################################################################
###################################################################################################


def is_run_folder(name):
    # Function to check if the folder name matches the 'run_N' pattern
    return re.match(r"^run_\d+$", name) is not None


def get_run_folders(batch_folder):

    # Ensure that the folders exist
    if not batch_folder.is_dir():
        print(f"Error: The input folder {batch_folder} does not exist.")
        exit(1)

    print(f" Batch run folder : {batch_folder}")

    # List all subdirectories in the batch_folder
    subdirs = [name.name for name in batch_folder.iterdir() if name.is_dir()]

    # Filter subdirectories that match the 'run_N' pattern
    run_folders = list(filter(is_run_folder, subdirs))

    # Count the run folders
    num_run_folders = len(run_folders)

    print(f"Number of 'run' subfolders: {num_run_folders}")

    return run_folders, num_run_folders


def generate_plot_tile(
    Array,
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_title,
    colmap,
    norm,
    cont_lvls,
    ax,
    i,
):

    # Formatting string for the contour labels
    cbar_fmt = lambda x, pos: "{:.3g}".format(x)

    im_c = ax.imshow(
        np.flipud(Array),
        extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        cmap=colmap,
        norm=norm,
    )

    if i > 0:
        cont_p = ax.contour(
            X_grid,
            Y_grid,
            Array,
            colors="gray",
            alpha=1.0,
            levels=cont_lvls,
            linewidths=1.0,
        )

        ax.clabel(cont_p, inline=True, fontsize=10)

    # contours of ceiling topography
    cont_H0 = ax.contour(X_grid, Y_grid, H0, colors="gray", alpha=0.25)

    # Plot injection locations
    for k in range(0, n_inj_locs):
        ax.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

    ax.set_title(str_title, fontsize=16)

    # Add a colourbar to the right of the plot, and make it the same height as the plot
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar_mob_m = plt.colorbar(
        im_c,
        ax=ax,
        cmap=colmap,
        norm=norm,
        cax=cax,
        format=FuncFormatter(cbar_fmt),
    )

    return im_c


def load_common_files(filepath_A, filepath_B, filename, skiprows_val_A, skiprows_val_B):

    file_A = Path(f"{filepath_A}/{filename}.txt")

    file_B = Path(f"{filepath_B}/{filename}.txt")

    try:
        data = np.loadtxt(file_A, skiprows=skiprows_val_A)
    except FileNotFoundError:
        # File isn't in the filepath_A folder - look in filepath_B instead
        try:
            data = np.loadtxt(file_B, skiprows=skiprows_val_B)
        except FileNotFoundError:
            print(f"{filename} was not found in either {filepath_A} or {filepath_B}.")

    # H0 = np.loadtxt(f"{filepath_B}/Input/ceil_topo.txt")

    return data


###################################################################################################
## Command-line arguments and file management #####################################################
###################################################################################################

# Deal with inputs for the batch run folder (built with help from ChatGPT 3.5)

# Initialize the argument parser
parser = argparse.ArgumentParser(description="Batch Run Data Folders")

# Add arguments for input and output folder with default values
parser.add_argument(
    "--batch", type=str, default="./batch_runs/", help="Batch run data folder path"
)

# Parse the arguments
args = parser.parse_args()

# Extract the main batch run folder path, and remove any trailing quotation marks and slashes
batch_folder = Path(args.batch)

# Get the run folders within
run_folders, num_run_folders = get_run_folders(batch_folder)


###################################################################################################
## Common data ####################################################################################
###################################################################################################

print("\n--Loading data--\n")

## Load common values from the run_1 folder
run_1_folder = batch_folder / "run_1"
common_files_folder = batch_folder / "Common_files"

# load plot times
plot_times = np.loadtxt(run_1_folder / "Output/Other/plot_times.txt")
n_plot = len(plot_times)


# # # #!!!
# # # #  Might I want these values to change between runs??
# # # #!!!
# # # # load parameter values
# # # parameters = np.loadtxt(
# # #     f"{run_1_folder}/Other/parameter_values.txt"
# # # )  # nx ny dx dy M Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve
# # # nx, ny, nz, dx, dy, M, Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve = parameters
# # # nx, ny, nz = [ int(i) for i in [nx,ny,nz] ] #Convert these to integers

# load parameter values
parameters = np.loadtxt(
    run_1_folder / "Output/Other/parameter_values.txt"
)  # nx ny nz dx dy M Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve
nx, ny, nz, dx, dy, M, Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve = parameters
nx, ny, nz = [int(i) for i in [nx, ny, nz]]


# Build the spatial grids
X_grid = np.arange(-(nx - 1) / 2.0, (nx - 1) / 2.0 + 1.0) * dx
Y_grid = np.arange(-(ny - 1) / 2.0, (ny - 1) / 2.0 + 1.0) * dy

X, Y = np.meshgrid(X_grid, Y_grid)


# load injection locations
inj_locs_data = load_common_files(
    common_files_folder,
    run_1_folder / "Output/Other",
    "injection_locations",
    3,
    0,
)
# inj_locs_data = np.loadtxt(f"{run_1_folder}/Output/Other/injection_locations.txt")

if np.ndim(inj_locs_data) == 1:
    # If there's only one injection point, this needs to be done slightly differently
    n_inj_locs = 1
    inj_grid_vals = np.zeros((1, 2))
    inj_grid_vals[0, 0] = X_grid[int(inj_locs_data[0])]
    inj_grid_vals[0, 1] = Y_grid[int(inj_locs_data[1])]
else:
    shape_inj_locs = np.shape(inj_locs_data)
    n_inj_locs = shape_inj_locs[0]
    inj_grid_vals = np.zeros((n_inj_locs, 2))
    for k in range(0, n_inj_locs):
        inj_grid_vals[k, 0] = X_grid[int(inj_locs_data[k, 0])]
        inj_grid_vals[k, 1] = Y_grid[int(inj_locs_data[k, 1])]


H0 = load_common_files(common_files_folder, run_1_folder / "Input", "ceil_topo", 0, 0)
B0 = load_common_files(common_files_folder, run_1_folder / "Input", "base_topo", 0, 0)
Porosity = load_common_files(
    common_files_folder, run_1_folder / "Input", "porosity", 0, 0
)
Permeability = load_common_files(
    common_files_folder, run_1_folder / "Input", "permeability", 0, 0
)

# # Reservoir Ceiling and Basement topography
# H0 = np.loadtxt(f"{run_1_folder}/Input/ceil_topo.txt")
# B0 = np.loadtxt(f"{run_1_folder}/Input/base_topo.txt")

# # Porosity and Permeability fields
# Porosity = np.loadtxt(f"{run_1_folder}/Input/porosity.txt")
# Permeability = np.loadtxt(f"{run_1_folder}/Input/permeability.txt")


print(np.shape(H0), np.shape(X))

###################################################################################################
## Remaining data #################################################################################
###################################################################################################

# Prior weightings for probabilistic combinations
# Uniform weighting
weights = (1.0 / num_run_folders) * np.ones(num_run_folders)

# Threshold value for value to be included/detected
V_threshold = 1e-5
h_presence_threshold = 1e-7
P_max_threshold = 1e0
P_min_threshold = -1e0


# Preallocate arrays

# CO2 Volumes
V_mobile_presence_array = np.zeros([ny, nx, n_plot])
V_total_presence_array = np.zeros([ny, nx, n_plot])
V_mobile_mean_array = np.zeros([ny, nx, n_plot])
V_total_mean_array = np.zeros([ny, nx, n_plot])
V_mobile_sd_array = np.zeros([ny, nx, n_plot])
V_total_sd_array = np.zeros([ny, nx, n_plot])
V_max_array = np.zeros([ny, nx])

# Mobile current thickness
h_presence_likelihood_array = np.zeros([ny, nx])

# Dynamic ambient pressure
P_max_array = np.zeros([ny, nx, n_plot])
P_min_array = np.zeros([ny, nx, n_plot])
P_mean_array = np.zeros([ny, nx, n_plot])
P_sd_array = np.zeros([ny, nx, n_plot])
P_overall_max_array = np.zeros([ny, nx])
P_overall_min_array = np.zeros([ny, nx])
P_max_likelihood_array = np.zeros([ny, nx])
P_min_likelihood_array = np.zeros([ny, nx])


# Fill in these arrays from the batch-run data
for k in range(0, num_run_folders):

    print(f"- {k+1} of {num_run_folders}")

    for l in range(1, n_plot):

        filename = f"run_{k+1}/Output/Current_Volume/V" + "{0:02d}".format(l) + ".txt"

        V = np.loadtxt(
            batch_folder / filename
        )

        filename = f"run_{k+1}/Output/Current_Volume/V_res" + "{0:02d}".format(l) + ".txt"

        V_res = np.loadtxt(
            batch_folder / filename
        )

        filename = f"run_{k+1}/Output/Current_Pressure/P" + "{0:02d}".format(l) + ".txt"

        P = np.loadtxt(
            batch_folder / filename
        )

        # Probability of mobile CO2 above threshold
        V_mobile_presence_array[:, :, l] += (V > V_threshold) * weights[k]
        # Probability of mobile+trapped CO2 above threshold
        V_total_presence_array[:, :, l] += ((V + V_res) > V_threshold) * weights[k]
        # Expected volume of mobile CO2
        V_mobile_mean_array[:, :, l] += V * weights[k]
        # Expected volume of mobile+trapped CO2
        V_total_mean_array[:, :, l] += (V + V_res) * weights[k]
        # Standard deviation of mobile CO2 (first part of calculation - Doing the E[X^2] part first)
        V_mobile_sd_array[:, :, l] += (V**2) * weights[k]
        # Standard deviation of mobile+trapped CO2 (first part of calculation)
        V_total_sd_array[:, :, l] += ((V + V_res) ** 2) * weights[k]

        # Maximum dynamic ambient pressure at this time across runs
        P_max_array[:, :, l] = np.maximum(P_max_array[:, :, l], P)
        # Minimum dynamic ambient pressure at this time across runs
        P_min_array[:, :, l] = np.minimum(P_min_array[:, :, l], P)
        # Expected dynamic ambient pressure
        P_mean_array[:, :, l] += P * weights[k]
        # Standard deviation of dynamic ambient pressure (first part of calculation - Doing the E[X^2] part first)
        P_sd_array[:, :, l] += (P**2) * weights[k]

    # Maximum volume of mobile+trapped CO2 at any time across all runs
    V_max_array[:, :] = np.maximum((V + V_res), V_max_array)

    # Probability of (mobile) CO2 reaching a given location
    Max_mobile_thickness = np.loadtxt(
        batch_folder / f"run_{k+1}/Output/Other/Max_mobile_thickness.txt"
    )

    h_presence_likelihood_array[:, :] = h_presence_likelihood_array[:, :] + weights[
        k
    ] * (Max_mobile_thickness > h_presence_threshold)

    # Max/min Pressure at a given location
    P_overall_max = np.loadtxt(
        batch_folder / f"run_{k+1}/Output/Other/Max_pressure.txt"
    )
    P_overall_min = np.loadtxt(
        batch_folder / f"run_{k+1}/Output/Other/Min_pressure.txt"
    )

    P_overall_max_array[:, :] = np.maximum(P_overall_max_array[:, :], P_overall_max)
    P_overall_min_array[:, :] = np.minimum(P_overall_min_array[:, :], P_overall_min)

    # Probability of dynamic pressure being above the threshold
    P_max_likelihood_array[:, :] = P_max_likelihood_array[:, :] + weights[k] * (
        P_overall_max > P_max_threshold
    )
    P_min_likelihood_array[:, :] = P_min_likelihood_array[:, :] + weights[k] * (
        P_overall_min < P_min_threshold
    )


# Second part of the standard deviation - subtracting the (E[X])^2 part and square-rooting
# Should there be a un-biasing prefactor in front of these? Are these sample statistics?

V_mobile_sd_array = np.sqrt(
    np.maximum(0.0, V_mobile_sd_array - V_mobile_mean_array**2)
)  # np.sqrt(num_run_folders/(num_run_folders-1)) *

V_total_sd_array = np.sqrt(
    np.maximum(0.0, V_total_sd_array - V_total_mean_array**2)
)  # np.sqrt(num_run_folders/(num_run_folders-1)) *


P_sd_array = np.sqrt(np.maximum(0.0, P_sd_array - P_mean_array**2))


# Time scale and units
if plot_times[-1] > 730.0:
    # Express in terms of standard years
    Time_scale = 365.25
    Time_unit = "years"
    plot_times = plot_times / Time_scale
else:
    Time_scale = 1.0
    Time_unit = "days"


###################################################################################################
## Plot preparation ###############################################################################
###################################################################################################

# Extreme values
V_value_max = np.max(V_max_array)
V_sd_max = np.max([np.max(V_mobile_sd_array), np.max(V_total_sd_array)])

P_value_max = np.max(P_overall_max_array)
P_value_min = np.min(P_overall_min_array)
P_sd_max = np.max(P_sd_array)

print(f"Volume: \n\tMax value: {V_value_max:4.3e} \n\tMax s.d.: {V_sd_max:4.3e}")
print(
    f"Pressure: \n\tMax value: {P_value_max:4.3e} \n\tMin value: {P_value_min:4.3e} \n\tMax s.d.: {P_sd_max:4.3e}"
)


# Colour map for values
colmap_Values = cm.Blues
norm_Values_V = colors.Normalize(0.0, V_value_max)
norm_Values_Pmax = colors.Normalize(0.0, P_value_max)
norm_Values_Pmin = colors.Normalize(P_value_min, 0.0)
norm_Values_P = colors.Normalize(P_value_min, P_value_max)

# Colour map for presence
colmap_Presence = cm.Reds
norm_Presence = colors.Normalize(0.0, 1.0)

# Colour map for standard deviation
colmap_SD = cm.Greens
norm_SD_V = colors.Normalize(0.0, V_sd_max)
norm_SD_P = colors.Normalize(0.0, P_sd_max)


# Formatting string for the contour labels
cbar_fmt = lambda x, pos: "{:.3g}".format(x)


###################################################################################################
## Volume Plots ###################################################################################
###################################################################################################

## Main loop over each of the output plots
# Start with the final plot time, so it can be viewed first

print("\n--Volume plots--\n")

for i in range(n_plot - 1, -1, -1):
    t = plot_times[i]
    # If there are no contours (because the function is zero everywhere),
    # supress the warning message
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="No contour levels were found within the data range."
        )

        print("Plot " + str(i))

        ##### Figure setup ########################################################################

        fig = plt.figure(figsize=(15, 10))

        gs = GridSpec(2, 4, figure=fig)

        # str_title = "[{:02d} of {:02d}]: t = {:3.2e}".format(i, n_plot - 1, t)
        str_title = f"[{i:02d} of {n_plot - 1:02d}]: t = {t:.2f} {Time_unit}"

        fig.suptitle(
            str_title,
            fontsize=16,
            x=0.01,
            y=0.99,
            horizontalalignment="left",
            verticalalignment="top",
            bbox=dict(facecolor="none", edgecolor="black"),
        )

        ##### Mobile Presence #####################################################################

        ax_c = fig.add_subplot(gs[0, 0])

        cont_lvls = [1e-1, 5e-1, 9e-1, 1e0]
        str_title = "Current Mobile CO$_{2}$ Presence"

        generate_plot_tile(
            V_mobile_presence_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Presence,
            norm_Presence,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Total Presence ######################################################################

        ax_c = fig.add_subplot(gs[1, 0])

        cont_lvls = [1e-1, 5e-1, 9e-1, 1e0]
        str_title = "Current Total CO$_{2}$ Presence"

        generate_plot_tile(
            V_total_presence_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Presence,
            norm_Presence,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Mobile Mean #########################################################################

        ax_c = fig.add_subplot(gs[0, 1])

        V_max_i = np.amax(V_mobile_mean_array[:, :, i])

        cont_lvls = [0.1 * V_max_i, 0.5 * V_max_i, 1.0 * V_max_i]
        str_title = "Mean Mobile Plume"

        generate_plot_tile(
            V_mobile_mean_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_V,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Total Mean ##########################################################################

        ax_c = fig.add_subplot(gs[1, 1])

        V_max_i = np.amax(V_total_mean_array[:, :, i])

        cont_lvls = [0.1 * V_max_i, 0.5 * V_max_i, 1.0 * V_max_i]
        str_title = "Mean Total Plume"

        generate_plot_tile(
            V_total_mean_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_V,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Mobile S.D. #########################################################################

        ax_c = fig.add_subplot(gs[0, 2])

        sd_max_i = np.max(V_mobile_sd_array[:, :, i])

        cont_lvls = [0.1 * sd_max_i, 0.5 * sd_max_i, 1.0 * sd_max_i]
        str_title = "S.D. Mobile Plume"

        generate_plot_tile(
            V_mobile_sd_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_SD,
            norm_SD_V,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Total S.D. ##########################################################################

        ax_c = fig.add_subplot(gs[1, 2])

        sd_max_i = np.max(V_total_sd_array[:, :, i])

        cont_lvls = [0.1 * sd_max_i, 0.5 * sd_max_i, 1.0 * sd_max_i]
        str_title = "S.D. Total Plume"

        generate_plot_tile(
            V_total_sd_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_SD,
            norm_SD_V,
            cont_lvls,
            ax_c,
            i,
        )

        ##### High Water Mark #####################################################################

        ax_c = fig.add_subplot(gs[0, 3])

        cont_lvls = [
            0.01 * V_value_max,
            0.1 * V_value_max,
            0.5 * V_value_max,
            1.0 * V_value_max,
        ]
        str_title = "High Water Mark"

        generate_plot_tile(
            V_max_array[:, :],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_V,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Presence Likelihood #################################################################

        ax_c = fig.add_subplot(gs[1, 3])

        cont_lvls = [1e-1, 5e-1, 9e-1, 1e0]
        str_title = r"CO$_{2}$ Presence Likelihood"

        generate_plot_tile(
            h_presence_likelihood_array[:, :],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Presence,
            norm_Presence,
            cont_lvls,
            ax_c,
            i,
        )

        # Add a contour showing where any CO2 has been observed across all of the runs
        cont_any = ax_c.contour(
            X_grid,
            Y_grid,
            (h_presence_likelihood_array[:, :] > 0.0).astype("float"),
            colors="k",
            alpha=1.0,
            levels=1,
            linewidths=0.5,
            linestyles="-",
        )

        ##### Save plot ###########################################################################

        # Tidy up the layout of subplots
        # plt.tight_layout(rect=[0, 0, 1, 0.95])
        gs.tight_layout(fig)

        # plt.show()

        # # Save this frame
        plt.savefig(
            f"./plots/batch_probability/temp/Volume/t{i:02d}.png"
        )
        plt.close()


###################################################################################################
## CO2 Likelihood Plot ############################################################################
###################################################################################################

clicked_points = []

fig = plt.figure()

gs = GridSpec(1, 1, figure=fig)

ax_c = fig.add_subplot(gs[0, 0])

cont_lvls = [1e-1, 5e-1, 9e-1, 1e0]
str_title = r"CO$_{2}$ Presence Likelihood"

im_c = generate_plot_tile(
    h_presence_likelihood_array[:, :],
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_title,
    colmap_Presence,
    norm_Presence,
    cont_lvls,
    ax_c,
    0,
)

# Add a contour showing where any CO2 has been observed across all of the runs
cont_any = ax_c.contour(
    X_grid,
    Y_grid,
    (h_presence_likelihood_array[:, :] > 0.0).astype("float"),
    colors="k",
    alpha=1.0,
    levels=1,
    linewidths=0.5,
    linestyles="-",
)

gs.tight_layout(fig)
# plt.show()

# # Save this frame
plt.savefig("./plots/batch_probability/temp/Volume/Volume_Likelihood_plot.png")
plt.close()


###################################################################################################
## Pressure Plots #################################################################################
###################################################################################################

## Main loop over each of the output plots
# Start with the final plot time, so it can be viewed first

print("\n--Pressure plots--\n")


for i in range(n_plot - 1, -1, -1):
    t = plot_times[i]
    # If there are no contours (because the function is zero everywhere),
    # supress the warning message
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="No contour levels were found within the data range."
        )

        print("Plot " + str(i))

        ##### Figure setup ########################################################################

        fig = plt.figure(figsize=(15, 10))

        gs = GridSpec(2, 4, figure=fig)

        str_title = "[{:02d} of {:02d}]: t = {:3.2e}".format(i, n_plot - 1, t)

        fig.suptitle(
            str_title,
            fontsize=16,
            x=0.01,
            y=0.98,
            horizontalalignment="left",
            verticalalignment="top",
            bbox=dict(facecolor="none", edgecolor="black"),
        )

        ##### Current Pressure Maximum ############################################################

        ax_c = fig.add_subplot(gs[0, 0])

        cont_lvls = [val * P_value_max for val in [0.0, 0.1, 0.5, 1.0]]
        str_title = "Current Max Pressure"

        generate_plot_tile(
            P_max_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_Pmax,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Mean Pressure #######################################################################

        ax_c = fig.add_subplot(gs[0, 1])

        P_max_i = np.amax(P_mean_array[:, :, i])

        cont_lvls = [val * P_max_i for val in [-1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0]]
        str_title = "Mean Pressure"

        generate_plot_tile(
            P_mean_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_P,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Pressure S.D. ##########################################################################

        ax_c = fig.add_subplot(gs[0, 2])

        sd_max_i = np.max(P_sd_array[:, :, i])

        cont_lvls = [val * sd_max_i for val in [-1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0]]
        str_title = "S.D. Pressure"

        generate_plot_tile(
            P_sd_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_SD,
            norm_SD_P,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Overall pressure maximum ############################################################

        ax_c = fig.add_subplot(gs[0, 3])

        cont_lvls = [val * P_value_max for val in [1e-1, 5e-1, 9e-1, 1e0]]
        str_title = r"Overall pressure maximum"

        generate_plot_tile(
            P_overall_max_array[:, :],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_Pmax,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Current Pressure Minimum ############################################################

        ax_c = fig.add_subplot(gs[1, 0])

        cont_lvls = [val * P_value_min for val in [1.0, 0.5, 0.1, 0.0]]
        str_title = "Current Min Pressure"

        generate_plot_tile(
            P_min_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_Pmin,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Mean Pressure #######################################################################

        ax_c = fig.add_subplot(gs[1, 1])

        P_max_i = np.amax(P_mean_array[:, :, i])

        cont_lvls = [val * P_max_i for val in [-1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0]]
        str_title = "Mean Pressure"

        generate_plot_tile(
            P_mean_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_P,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Pressure S.D. ##########################################################################

        ax_c = fig.add_subplot(gs[1, 2])

        sd_max_i = np.max(P_sd_array[:, :, i])

        cont_lvls = [val * sd_max_i for val in [-1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0]]
        str_title = "S.D. Pressure"

        generate_plot_tile(
            P_sd_array[:, :, i],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_SD,
            norm_SD_P,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Overall pressure minimum ############################################################

        ax_c = fig.add_subplot(gs[1, 3])

        cont_lvls = [val * P_value_min for val in [1e0, 9e-1, 5e-1, 1e-1]]
        str_title = r"Overall pressure minimum"

        generate_plot_tile(
            P_overall_min_array[:, :],
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            colmap_Values,
            norm_Values_Pmin,
            cont_lvls,
            ax_c,
            i,
        )

        ##### Save plot ###########################################################################

        # Tidy up the layout of subplots
        # plt.tight_layout(rect=[0, 0, 1, 0.95])
        gs.tight_layout(fig)

        # plt.show()

        # # Save this frame
        plt.savefig(
            "./plots/batch_probability/temp/Pressure/t" + "{0:02d}".format(i) + ".png"
        )
        plt.close()


###################################################################################################
## Pressure Max and Min Plot ######################################################################
###################################################################################################

fig = plt.figure(figsize=(8, 6))

gs = GridSpec(2, 2, figure=fig)

ax_c = fig.add_subplot(gs[0, 0])

cont_lvls = [1e-1, 5e-1, 9e-1, 1e0]
str_title = r"Pressure Max Exceedance Likelihood"

generate_plot_tile(
    P_max_likelihood_array[:, :],
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_title,
    colmap_Presence,
    norm_Presence,
    cont_lvls,
    ax_c,
    0,
)


ax_c = fig.add_subplot(gs[0, 1])

cont_lvls = [val * P_value_max for val in [1e-1, 5e-1, 9e-1, 1e0]]
str_title = r"Overall pressure maximum"

generate_plot_tile(
    P_overall_max_array[:, :],
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_title,
    colmap_Values,
    norm_Values_Pmax,
    cont_lvls,
    ax_c,
    1,
)


ax_c = fig.add_subplot(gs[1, 0])

cont_lvls = [1e-1, 5e-1, 9e-1, 1e0]
str_title = r"Pressure Min Exceedance Likelihood"

generate_plot_tile(
    P_min_likelihood_array[:, :],
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_title,
    colmap_Presence,
    norm_Presence,
    cont_lvls,
    ax_c,
    0,
)


ax_c = fig.add_subplot(gs[1, 1])

cont_lvls = [val * P_value_min for val in [1e0, 9e-1, 5e-1, 1e-1]]
str_title = r"Overall pressure minimum"

generate_plot_tile(
    P_overall_min_array[:, :],
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_title,
    colmap_Values,
    norm_Values_Pmin,
    cont_lvls,
    ax_c,
    1,
)


gs.tight_layout(fig)
# plt.show()

# # Save this frame
plt.savefig("./plots/batch_probability/temp/Pressure/Pressure_Likelihood_plot.png")
plt.close()
