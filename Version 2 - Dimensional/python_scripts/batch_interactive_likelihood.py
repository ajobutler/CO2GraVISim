# batch_interactive_likelihood.py (01/07/24)
# This script takes the output from a batch run of CO2GraVISim
# and produces an interactive plot of where the CO2 can reach, so that
# the user can focus in on which runs see CO2 at specific locations

# Built with the help of
# https://matplotlib.org/stable/users/explain/figure/event_handling.html
# https://stackoverflow.com/questions/48066779/get-the-indices-of-the-clicked-point-relative-to-a-subplot


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec

import argparse
from pathlib import Path
import re  # Regular expressions

import functools
from matplotlib.widgets import RadioButtons

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

    print(f"Batch run folder : {batch_folder}")

    # List all subdirectories in the batch_folder
    subdirs = [name.name for name in batch_folder.iterdir() if name.is_dir()]

    # Filter subdirectories that match the 'run_N' pattern
    run_folders = list(filter(is_run_folder, subdirs))

    # Count the run folders
    num_run_folders = len(run_folders)

    # Sort the folders to avoid lexicographical issues, i.e. so that
    # run_2 comes before run_10
    run_folders.sort(key=lambda x: int(re.search(r"\d+", x).group()))

    print(f"Number of 'run' subfolders: {num_run_folders}")

    return run_folders, num_run_folders


def onclick_function(fig, axs, event):
    # Function to determine information based on where the user has clicked.

    global h_point_L, h_point_R  # global so that these values are remembered between calls
    global specific_runs_L, specific_runs_R

    if axs["main"].in_axes(event):

        if event.button == 1:
            print("\n[Left click]")
            h_point_L, specific_runs_L = plot_click_update(
                event,
                h_point_L,
                fig,
                "Likelihood_L",
                "Detail_L",
                "gs",
            )

        elif event.button == 3:
            print("\n[Right click]")
            h_point_R, specific_runs_R = plot_click_update(
                event,
                h_point_R,
                fig,
                "Likelihood_R",
                "Detail_R",
                "bo",
            )

    axs["iterations"].clear()
    run_indentifier_plot(axs["iterations"])
    fig.canvas.draw()
    fig.canvas.flush_events()

    return True


def plot_click_update(event, h_point, fig, str_target_p, str_target_d, marker_style):
    # This function calculates the indices of the point chosen on the main plot,
    # and updates the appropriate other plots depending on the mouse button used
    # (i.e. left click or right click)

    # Find nearest point on x grid
    # searchsorted (from the left) returns i such that a[i-1] < v <= a[i]
    idx_x = np.searchsorted(X_grid, event.xdata, side="left")
    idx_y = np.searchsorted(Y_grid, event.ydata, side="left")

    x_val = X_grid[idx_x]
    y_val = Y_grid[idx_y]

    # Add a marker of the chosen point to the plot
    if h_point is None:
        (h_point,) = axs["main"].plot(x_val, y_val, marker_style)
        fig.canvas.draw()
        fig.canvas.flush_events()
    else:
        h_point.set_data([x_val], [y_val])
        fig.canvas.draw()
        fig.canvas.flush_events()

    # Determine which runs have sufficient CO2 at this location
    pass_runs = np.where(
        Max_mobile_thickness_array[idx_y, idx_x, :] > h_presence_threshold
    )
    pass_runs = pass_runs[0].tolist()
    # pass_runs[:] = [x+1 for x in pass_runs]

    # A boolean list of 'pass'/'fail' for each run
    specific_runs = np.zeros(num_run_folders, dtype=bool)
    specific_runs[pass_runs] = True

    detail_plot(specific_runs, str_target_p, str_target_d, marker_style)

    print(f"Indices (i,j) = ({idx_x:d} , {idx_y:d})")
    print(f"Coordinates (x,y) = ({x_val:.3f} , {y_val:.3f})")
    print(f"Runs:")
    print([x + 1 for x in pass_runs])  # Switch from 0-indexed to 1-indexed
    print("\n")

    return h_point, specific_runs


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
    str_ax,
    i,
):
    # This function produces a 2D plot coloured by a given array, with a specified
    # colourscheme, ceiling topography contours, injection locations, and a colourbar

    global cbaxs

    # Formatting string for the contour labels
    cbar_fmt = lambda x, pos: "{:.3g}".format(x)

    im_c = axs[str_ax].imshow(
        np.flipud(Array),
        extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        cmap=colmap,
        norm=norm,
    )

    if i > 0:
        cont_p = axs[str_ax].contour(
            X_grid,
            Y_grid,
            Array,
            colors="gray",
            alpha=1.0,
            levels=cont_lvls,
            linewidths=2,
        )

        axs[str_ax].clabel(cont_p, inline=True, fontsize=10)

    # contours of ceiling topography
    cont_H0 = axs[str_ax].contour(X_grid, Y_grid, H0, colors="gray", alpha=0.25)

    # Plot injection locations
    for k in range(0, n_inj_locs):
        axs[str_ax].plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

    axs[str_ax].set_title(str_title, fontsize=16)

    # Add a colourbar to the right of the plot, and make it the same height as the plot
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="5%", pad=0.05)

    if cbaxs[str_ax] is not None:
        # Update existing colourbar
        cbaxs[str_ax].update_normal(im_c)
    else:
        # Create colourbar
        cbaxs[str_ax] = plt.colorbar(
            im_c,
            ax=axs[str_ax],
            cmap=colmap,
            norm=norm,
            aspect=50,
            format=FuncFormatter(cbar_fmt),
        )

    return im_c


def detail_plot(specific_runs, str_target_p, str_target_d, marker_style):
    # The conditional probability and mean topography arrays just for the specific runs chosen
    Likelihood_array, h_max_array, Topo_array, kappa_c_array = specific_run_properties(
        specific_runs
    )

    str_title = f"{{{marker_style}}}: {np.sum(specific_runs)} runs"

    # Update the relevant axes
    axs[str_target_p].clear()
    axs[str_target_d].clear()

    likelihood_plot(Likelihood_array, str_title, str_target_p)

    if plot_mode == str_topo:
        topography_plot(Likelihood_array, Topo_array, str_title, str_target_d)
    elif plot_mode == str_perm:
        avg_perm_plot(Likelihood_array, kappa_c_array, str_title, str_target_d)
    elif plot_mode == str_hmax:
        hmax_plot(Likelihood_array, h_max_array, str_title, str_target_d)
    elif plot_mode == str_sweep:
        sweep_plot(specific_runs, str_title, str_target_d, marker_style)
    elif plot_mode == str_areal:
        areal_extent_plot(specific_runs, str_title, str_target_d, marker_style)
    elif plot_mode == str_dist:
        distance_plot(specific_runs, str_title, str_target_d, marker_style)

    return True


def specific_run_properties(specific_runs):
    # This function calculates the likelihood array and mean topography array
    # just for the specific runs that pass the threshold at a given point

    # Calculate appropriate weighted average
    N_specific_runs = np.sum(specific_runs)

    if N_specific_runs > 0:

        # Prior weightings for probabilistic combinations
        # Uniform weighting
        weights = (1.0 / N_specific_runs) * specific_runs

        Likelihood_array = np.sum(
            (Max_mobile_thickness_array > h_presence_threshold)
            * weights.reshape([1, 1, -1]),
            axis=2,
        )

        h_max_array = np.sum(
            Max_mobile_thickness_array * weights.reshape([1, 1, -1]),
            axis=2,
        )

        Topo_array = np.sum(
            Topography_variation * weights.reshape([1, 1, -1]),
            axis=2,
        )

        kappa_c_array = np.sum(
            kappa_c_variation * weights.reshape([1, 1, -1]),
            axis=2,
        )

    else:
        # No runs passed the threshold at this point
        Likelihood_array = np.zeros([ny, nx])
        h_max_array = np.zeros([ny, nx])
        Topo_array = np.zeros([ny, nx])
        kappa_c_array = np.zeros([ny, nx])

    return Likelihood_array, h_max_array, Topo_array, kappa_c_array


def likelihood_plot(Likelihood_array, str_title, str_ax):
    # This function creates a likelihood plot for a specified set of data,
    # and plots it to the appropriate axes.
    # specific_runs should be a 1D array of boolean values indicating whether a specific
    # run should be included

    str_title = str_title + " - cond. prob."

    # Plot
    im_c = generate_plot_tile(
        Likelihood_array[:, :],
        X_grid,
        Y_grid,
        H0,
        n_inj_locs,
        inj_grid_vals,
        str_title,
        colmap_Presence,
        norm_Presence,
        cont_lvls,
        str_ax,
        0,
    )

    # Add a contour showing where any CO2 has been observed across the specied runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (Likelihood_array[:, :] > 0.0).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles="-",
    )

    # Add a contour showing where any CO2 has been observed across all runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (
            np.sum(Max_mobile_thickness_array > h_presence_threshold, axis=2) > 0.0
        ).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles=":",
    )

    return im_c


def hmax_plot(Likelihood_array, h_max_array, str_title, str_ax):
    # This function creates a plot of the average maximum mobile CO2 thickness
    # for a specified set of data, and plots it to the appropriate axes.
    # specific_runs should be a 1D array of boolean values indicating whether a specific
    # run should be included

    str_title = str_title + " - h_max"

    hmax_max = np.max(h_max_array)
    norm_hmax = colors.Normalize(0.0, hmax_max)

    hmax_lvls = np.array([0.1, 0.25, 0.5, 0.75, 0.9]) * hmax_max

    # Plot
    im_c = generate_plot_tile(
        h_max_array[:, :],
        X_grid,
        Y_grid,
        H0,
        n_inj_locs,
        inj_grid_vals,
        str_title,
        colmap_hmax,
        norm_hmax,
        [],  # cont_lvls,
        str_ax,
        0,
    )

    # Add a contour showing where any CO2 has been observed across the specied runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (Likelihood_array[:, :] > 0.0).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles="-",
    )

    # Add percentage contours for this mean of the specified runs
    cont_hmax = axs[str_ax].contour(
        X_grid,
        Y_grid,
        h_max_array[:, :],
        colors="b",
        alpha=0.5,
        levels=hmax_lvls,
        linewidths=0.5,
        linestyles="--",
    )

    # Add contour labels for this specific-run-average hmax plot
    axs[str_ax].clabel(cont_hmax, inline=True, fmt=f"%3.2f", fontsize=10)

    # # # Add horizontal marker lines to the colourbar corresponding to the current contours
    # # # cbaxs[str_ax].ax.axhline(max_val, c="r")
    # # for _, L in enumerate(hmax_lvls):
    # #     cbaxs[str_ax].ax.axhline(L, c="b", ls="--")

    ## Currently these horizontal lines cause an issue, because they persist when switching
    ## to different runs. I need to work out how to appropriately update these / clear the
    ## colourbar

    # Add a contour showing where any CO2 has been observed across all runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (
            np.sum(Max_mobile_thickness_array > h_presence_threshold, axis=2) > 0.0
        ).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles=":",
    )

    return im_c


def avg_perm_plot(Likelihood_array, Avg_perm_array, str_title, str_ax):
    # This function creates a plot of the mean average permeability variation
    # for a specified set of data, and plots it to the appropriate axes.
    # specific_runs should be a 1D array of boolean values indicating whether a specific
    # run should be included

    str_title = str_title + " - avg. perm."

    # Plot
    im_c = generate_plot_tile(
        Avg_perm_array[:, :] * (Likelihood_array[:, :] > 0.0),
        X_grid,
        Y_grid,
        H0,
        n_inj_locs,
        inj_grid_vals,
        str_title,
        colmap_Permeability,
        norm_Permeability,
        cont_lvls_Avg_Perm,
        str_ax,
        0,
    )

    # Add a contour showing where any CO2 has been observed across the specified runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (Likelihood_array[:, :] > 0.0).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles="-",
    )

    # Add a contour showing where any CO2 has been observed across all runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (
            np.sum(Max_mobile_thickness_array > h_presence_threshold, axis=2) > 0.0
        ).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles=":",
    )

    return im_c


def topography_plot(Likelihood_array, Topo_array, str_title, str_ax):
    # This function creates a plot of the mean topography variation
    # for a specified set of data and plots it to the appropriate axes.
    # specific_runs should be a 1D array of boolean values indicating
    # whether a specific run should be included

    str_title = str_title + " - mean topo."

    # Plot
    im_c = generate_plot_tile(
        Topo_array[:, :],
        X_grid,
        Y_grid,
        H0,
        n_inj_locs,
        inj_grid_vals,
        str_title,
        colmap_Topography,
        norm_Topography,
        cont_lvls_Topography,
        str_ax,
        0,
    )

    # Add a contour showing where any CO2 has been observed across the specified runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (Likelihood_array[:, :] > 0.0).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles="-",
    )

    # Add a contour showing where any CO2 has been observed across all runs
    axs[str_ax].contour(
        X_grid,
        Y_grid,
        (
            np.sum(Max_mobile_thickness_array > h_presence_threshold, axis=2) > 0.0
        ).astype("float"),
        colors="k",
        alpha=1.0,
        levels=1,
        linewidths=0.5,
        linestyles=":",
    )

    return im_c


def sweep_plot(specific_runs, str_title, str_ax, marker_style):

    I = np.arange(0, num_run_folders)

    axs[str_ax].plot(I, Sweep_array, ls="-", c="b")

    axs[str_ax].scatter(
        I[specific_runs],
        Sweep_array[specific_runs],
        marker=marker_style[1],
        c=marker_style[0],
    )

    axs[str_ax].set_xlim([-1, num_run_folders])
    axs[str_ax].set_ylim([0, 1.05 * Sweep_scale])

    axs[str_ax].set_aspect("auto")

    axs[str_ax].minorticks_on()
    axs[str_ax].grid(which="both")

    axs[str_ax].set_xlabel("Run")

    axs[str_ax].set_ylabel("Sweep (%)")

    axs[str_ax].set_title(str_title, fontsize=16)


def areal_extent_plot(specific_runs, str_title, str_ax, marker_style):

    I = np.arange(1, num_run_folders + 1)

    axs[str_ax].plot(I, Areal_extent_array, ls="-", c="b")

    axs[str_ax].scatter(
        I[specific_runs],
        Areal_extent_array[specific_runs],
        marker=marker_style[1],
        c=marker_style[0],
    )

    axs[str_ax].set_xlim([0, num_run_folders])
    axs[str_ax].set_ylim([0, 1.05 * Areal_extent_scale])

    axs[str_ax].set_aspect("auto")

    axs[str_ax].minorticks_on()
    axs[str_ax].grid(which="both")

    axs[str_ax].set_xlabel("Run")

    axs[str_ax].set_ylabel("Areal Extent (%)")

    axs[str_ax].set_title(str_title, fontsize=16)


def distance_plot(specific_runs, str_title, str_ax, marker_style):

    I = np.arange(1, num_run_folders + 1)

    axs[str_ax].plot(I, Distance_array, ls="-", c="b")

    axs[str_ax].scatter(
        I[specific_runs],
        Distance_array[specific_runs],
        marker=marker_style[1],
        c=marker_style[0],
    )

    axs[str_ax].set_xlim([0, num_run_folders])
    axs[str_ax].set_ylim([0, 1.05 * Distance_scale])

    axs[str_ax].set_aspect("auto")

    axs[str_ax].minorticks_on()
    axs[str_ax].grid(which="both")

    axs[str_ax].set_xlabel("Run")

    axs[str_ax].set_ylabel("Distance")

    axs[str_ax].set_title(str_title, fontsize=16)


def read_in_3D_array(folder_A, folder_B, filename, nx, ny, nz):
    # This function reads in a 3D array. These are stored in files as
    # nz consecutive chunks of (nx*ny) data

    Array = np.zeros((ny, nx, nz))

    filepath_A = folder_A / f"{filename}.txt"
    filepath_B = folder_B / f"{filename}.txt"

    if filepath_A.is_file():
        filepath = filepath_A
    elif filepath_B.is_file():
        filepath = filepath_B
    else:
        raise FileNotFoundError(f"{filename} not found in {folder_A} or {folder_B}.")

    with open(filepath) as f:
        Lines = f.readlines()

        l_count = 0
        for k in range(0, nz):
            for j in range(0, ny):
                L = Lines[l_count].strip()
                L = L.split()
                Array[j, :, k] = [float(x) for x in L]
                l_count += 1

    return Array


def calculate_cumulative_integral(F, D0):

    F_Cint = np.zeros((ny, nx, nz))

    Z_array = np.zeros((ny, nx, nz))

    for k in range(0, nz):
        Z_array[:, :, k] = (k / (nz - 1)) * (D0[:, :])

    Sum = np.zeros((ny, nx))

    for k in range(1, nz):
        dz = Z_array[:, :, k] - Z_array[:, :, k - 1]
        Sum = Sum + 0.5 * dz[:, :] * (F[:, :, k - 1] + F[:, :, k])
        F_Cint[:, :, k] = Sum

    return F_Cint, Z_array


def vertical_interpolation(F, h, D0, Z_array):

    F_interp = np.zeros((ny, nx))

    xi = np.divide(h, D0)

    for j in range(0, ny):
        for i in range(0, nx):
            z = np.squeeze(Z_array[j, i, :])
            f = np.squeeze(F[j, i, :])
            F_interp[j, i] = np.interp(xi[j, i], z, f)

    return F_interp


def perm_averages(Permeability, h, D0):

    perm_Cint, Z_array = calculate_cumulative_integral(Permeability, D0)

    CO2_perm_int = vertical_interpolation(perm_Cint, h, D0, Z_array)

    Amb_perm_int = np.zeros([ny, nx])
    Amb_perm_int[:, :] = perm_Cint[:, :, -1] - CO2_perm_int[:, :]

    kappa_c = integral_average(CO2_perm_int, h)
    kappa_a = integral_average(Amb_perm_int, (D0 - h))

    return kappa_c, kappa_a


def integral_average(Integral, z):
    # This function calculates the integral average
    # Avg = (1/z) * Integral, where Integral =  \int_{0}^{z} F(u) du
    # We need to make sure we deal with z=0 appropriately here.
    # Integral and z should be the same shape.

    Avg = np.zeros_like(Integral)

    mask = z > 1e-5
    Avg[mask] = Integral[mask] / z[mask]

    return Avg


def Array_variation(Array):
    # This function calculates the variation of a 3D array about its
    # mean over its 3rd axis

    Array_variation = Array - np.reshape(np.mean(Array, axis=2), [ny, nx, 1])

    return Array_variation


def run_indentifier_plot(ax):
    # This function produces a scatter plot indicating which runs are in the
    # two groups identified by clicking on the main plot

    Indices = np.arange(0, num_run_folders)
    Indices_ticks = [x + 1 for x in Indices]

    mrkr_style_L = ["s" if s else "x" for s in specific_runs_L]
    mrkr_colour_L = ["g" if s else "k" for s in specific_runs_L]

    mrkr_style_R = ["o" if s else "x" for s in specific_runs_R]
    mrkr_colour_R = ["b" if s else "k" for s in specific_runs_R]

    for k in Indices:
        ax.scatter(
            -1.0,
            Indices[k] + 1,  # 1-indexed
            marker=mrkr_style_L[k],
            c=mrkr_colour_L[k],
        )
        ax.scatter(
            +1.0,
            Indices[k] + 1,  # 1-indexed
            marker=mrkr_style_R[k],
            c=mrkr_colour_R[k],
        )

    ax.set_xlim([-2, 2])

    ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax.set_yticks(Indices_ticks)


def update_plot_mode(label):

    global plot_mode

    # Update plot mode string
    plot_mode = label

    # #Refresh detail plots, if they exist
    if h_point_L is not None:
        detail_plot(specific_runs_L, "Likelihood_L", "Detail_L", "gs")
    if h_point_R is not None:
        detail_plot(specific_runs_R, "Likelihood_R", "Detail_R", "bo")

    fig.canvas.draw()
    fig.canvas.flush_events()


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
plot_times = np.loadtxt( run_1_folder / "Output/Other/plot_times.txt")
n_plot = len(plot_times)

# load parameter values
parameters = np.loadtxt(
    run_1_folder / "Output/Other/parameter_values.txt"
)  # nx ny nz dx dy M Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve
# Assign each to their own variables
nx, ny, nz, dx, dy, M, Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve = parameters
nx, ny, nz = [int(i) for i in [nx, ny, nz]]  # Convert these ones to integers


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


###################################################################################################
## Remaining data #################################################################################
###################################################################################################

# Threshold value for value to be included/detected
h_presence_threshold = 1e-7

# Preallocate arrays
Max_mobile_thickness_array = np.zeros([ny, nx, num_run_folders])
h_presence_likelihood_array = np.zeros([ny, nx])
Topography_variation = np.zeros([ny, nx, num_run_folders])
kappa_c_variation = np.zeros([ny, nx, num_run_folders])
kappa_a_variation = np.zeros([ny, nx, num_run_folders])

D0_array = np.zeros([ny, nx, num_run_folders])

Average_permeability_array = np.zeros([ny, nx, nz])

Sweep_array = np.zeros([num_run_folders])

Areal_extent_array = np.zeros([num_run_folders])

Distance_array = np.zeros([num_run_folders])

# Fill in these arrays from the batch-run data
for k in range(0, num_run_folders):

    print(f"- {k+1} of {num_run_folders}")

    run_folder = batch_folder / run_folders[k]

    # Maximum mobile CO2 thickness at each location over the run
    h_max = np.loadtxt(run_folder / "Output/Other/Max_mobile_thickness.txt")

    # Reservoir topography
    H0 = load_common_files(
        common_files_folder, run_folder / "Input", "ceil_topo", 0, 0
    )
    B0 = load_common_files(
        common_files_folder, run_folder / "Input", "base_topo", 0, 0
    )
    D0 = B0 - H0

    # Store relevant values
    Max_mobile_thickness_array[:, :, k] = h_max
    Topography_variation[:, :, k] = H0
    D0_array[:, :, k] = D0

    # Permeability for this run
    Permeability = read_in_3D_array(
        common_files_folder, run_folder / "Input", "permeability", nx, ny, nz
    )

    Average_permeability_array[:, :, :] += (1.0 / num_run_folders) * Permeability[
        :, :, :
    ]

    # Vertical averages for the CO2 and ambient based on the max mobile CO2 thickness
    kappa_c_variation[:, :, k], kappa_a_variation[:, :, k] = perm_averages(
        Permeability, h_max, D0
    )

    # Percentage sweep of the entire storage interval
    Sweep_array[k] = 100.0 * np.sum(np.divide(h_max, D0)) / (nx * ny)

    # Percentage coverage of the 2D areal footprint of the storage interval
    Areal_extent_array[k] = 100.0 * np.sum(h_max > h_presence_threshold) / (nx * ny)

    # Maximum distance of significant CO2 presence from the(!) injection location
    Distance_array[k] = np.max(
        (h_max > h_presence_threshold)
        * np.sqrt((X - inj_grid_vals[0, 0]) ** 2 + (Y - inj_grid_vals[0, 1]) ** 2)
    )
    # !!! This is currently only built for a single injection well


# Convert these to the variation around their respective means
Topography_variation = Array_variation(Topography_variation)


# Convert kappa values to variation compared to integral done over the run-average permeability
for k in range(0, num_run_folders):

    kappa_c_avgperm, kappa_a_avgperm = perm_averages(
        Average_permeability_array,
        Max_mobile_thickness_array[:, :, k],
        D0_array[:, :, k],
    )

    kappa_c_variation[:, :, k] -= kappa_c_avgperm
    kappa_a_variation[:, :, k] -= kappa_a_avgperm


# Scales for colourbars and colourschemes
Topography_scale = np.max(np.abs(Topography_variation))

kappa_c_variation_scale = np.max(kappa_c_variation)
kappa_a_variation_scale = np.max(kappa_a_variation)

avg_perm_variation_scale = np.max([kappa_c_variation_scale, kappa_a_variation_scale])

Sweep_scale = np.max(Sweep_array)

Areal_extent_scale = np.max(Areal_extent_array)

Distance_scale = np.max(Distance_array)

###################################################################################################
## Plot preparation ###############################################################################
###################################################################################################

# Colour map for presence
colmap_Presence = cm.Reds
norm_Presence = colors.Normalize(0.0, 1.0)

# Colour map for hmax
colmap_hmax = cm.Reds
# Norm is set for each specific run, so that features are clear

# Colour map for topography
colmap_Topography = cm.seismic
norm_Topography = colors.Normalize(-Topography_scale, Topography_scale)

# Colour map for permeability
colmap_Permeability = cm.seismic
norm_Permeability = colors.Normalize(
    -avg_perm_variation_scale, avg_perm_variation_scale
)

# Formatting string for the contour labels
cbar_fmt = lambda x, pos: "{:.3g}".format(x)

cont_lvls = [1e-1, 5e-1, 9e-1, 1e0]
# cont_lvls_Topography = [ i*np.max(Topography_variation) for i in [1e-1, 5e-1, 9e-1, 1e0] ]
cont_lvls_Topography = np.linspace(-Topography_scale, Topography_scale, 10)
cont_lvls_Avg_Perm = np.linspace(
    -avg_perm_variation_scale, avg_perm_variation_scale, 10
)

# plot marker handles
h_point_L = None
h_point_R = None

# Initial values for the specific_runs arrays - these are boolean arrays
# identifying which runs belong to which plot sets
specific_runs = np.ones(num_run_folders, dtype=bool)
specific_runs_L = np.zeros(num_run_folders, dtype=bool)
specific_runs_R = np.zeros(num_run_folders, dtype=bool)

# Plotting mode label
str_topo = "Topo."
str_perm = "Perm."
str_hmax = "hmax"
str_sweep = "Sweep"
str_areal = "Area"
str_dist = "Distance"

plot_mode = str_sweep

radio_strings = (str_topo, str_hmax, str_perm, str_sweep, str_areal, str_dist)
radio_clrs = ["k" for s in radio_strings]
radio_fnts = [8 for s in radio_strings]
radio_sizes = [15 for s in radio_strings]

radio_active = 3  # Areal

# Radio buttons
radio_bg_clr = "gainsboro"

###################################################################################################
## Plot ###########################################################################################
###################################################################################################

fig, axs = plt.subplot_mosaic(
    [
        ["buttons", "iterations", "main", "Likelihood_L", "Likelihood_R"],
        [".", "iterations", "main", "Likelihood_L", "Likelihood_R"],
        [".", "iterations", "main", "Detail_L", "Detail_R"],
        [".", "iterations", "main", "Detail_L", "Detail_R"],
    ],
    width_ratios=[0.2, 0.05, 1, 0.75, 0.75],
    layout="constrained",
    figsize=(15, 8),
)

# Associated colourbars
cbaxs = {
    "main": None,
    "Likelihood_L": None,
    "Likelihood_R": None,
    "Detail_L": None,
    "Detail_R": None,
}

str_title = r"CO$_{2}$ Presence Likelihood"
fig.suptitle(str_title, fontsize=18)

# Button
axs["buttons"].set_facecolor(radio_bg_clr)
radio = RadioButtons(
    axs["buttons"],
    radio_strings,
    active=radio_active,
    label_props={"color": radio_clrs, "fontsize": radio_fnts},
    radio_props={"s": radio_sizes},
)
radio.on_clicked(update_plot_mode)


# Run identifier plot
run_indentifier_plot(axs["iterations"])


# Main plot, for interaction
str_title = f"All {num_run_folders} runs"
specific_runs = np.ones(num_run_folders, dtype=bool)
Likelihood_array, _, _, _ = specific_run_properties(specific_runs)
im_c = likelihood_plot(Likelihood_array, str_title, "main")

# Response to a mouse click
onclick = functools.partial(onclick_function, fig, axs)
cid = fig.canvas.mpl_connect("button_press_event", onclick)

# fig.set_tight_layout(True)
# fig.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

plt.show()
