# Overview_plot.py
# This script produces plots from the output of CO2GraVISim,
# showing
# - the Permeability and Porosity fields,
# - the ceiling topography,
# - the local thickness and volumes of the mobile and trapped CO2
# - the pressure fields for the CO2 and the ambient fluid
# - the total volumes of free and trapped CO2 over time

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec
import warnings

from scipy import integrate, interpolate


import argparse
from pathlib import Path

###################################################################################################
## Functions ######################################################################################
###################################################################################################


def read_in_3D_array(filepath, nx, ny, nz):
    # This function reads in a 3D array. These are stored in files as nz consecutive chunks
    # of (nx*ny) data

    Array = np.zeros((ny, nx, nz))

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


def calculate_cumulative_integral(F):

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


def vertical_interpolation(F, h):

    F_interp = np.zeros((ny, nx))

    xi = np.divide(h, D0)

    for j in range(0, ny):
        for i in range(0, nx):
            z = np.squeeze(Z_array[j, i, :])
            f = np.squeeze(F[j, i, :])
            # interp_fn = interpolate.interp1d(z, f, kind="linear")
            # F_interp[j, i] = interp_fn(xi[j, i])
            F_interp[j, i] = np.interp(xi[j, i], z, f)

    return F_interp


def colourmap_arrays(Array, colmap, minn, maxx):
    # This function produces various colourmap arrays used for the 2D and
    # 3D plots generated here.

    norm = colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap=colmap)
    m.set_array([])
    fcolors = m.to_rgba(Array)

    return fcolors, norm, m


def generate_pressure_plot(
    P,
    Ux,
    Uy,
    V,
    V_res,
    colmap,
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_type,
    ax,
):
    # This function produces a 2D plot of a pressure array and the corresponding Darcy velocity vectors field,
    # along with relevant contours and the ceiling topography

    min_val, max_val = [-np.max(np.abs(P)), np.max(np.abs(P))]

    # Formatting string for the contour labels
    cbar_fmt = lambda x, pos: "{:.3e}".format(x)

    # Masked array
    # masked_array = np.ma.masked_where(P < threshold, P)
    _, norm, _ = colourmap_arrays(P, colmap, min_val, max_val)

    # surface plot
    im_c = ax.imshow(
        np.flipud(P),
        cmap=colmap_P,
        norm=norm,
        extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
    )
    # ambient Darcy velocity field
    ax.quiver(
        X[0:-1:flow_field_step, 0:-1:flow_field_step],
        Y[0:-1:flow_field_step, 0:-1:flow_field_step],
        Uy[0:-1:flow_field_step, 0:-1:flow_field_step],
        Ux[0:-1:flow_field_step, 0:-1:flow_field_step],
    )
    # contours of ceiling topography
    cont_H0 = ax.contour(X_grid, Y_grid, H0, colors="grey", alpha=0.25)
    ax.clabel(cont_H0, inline=True, fmt="H0 = %.2g", fontsize=10)

    # edge contour of mobile current
    cont_Vb = ax.contour(
        X_grid,
        Y_grid,
        V * (V >= V_threshold),
        colors=clr_V,
        alpha=1.0,
        levels=[V_threshold],
        linewidths=3,
    )
    #
    # edge contour of trapped current
    cont_Vresb = ax.contour(
        X_grid,
        Y_grid,
        V_res * (V_res >= V_res_threshold),
        colors=clr_V_res,
        alpha=1.0,
        levels=[V_res_threshold],
        linewidths=3,
    )

    # Plot injection locations
    for k in range(0, n_inj_locs):
        ax.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "kx")

    ax.set_title(str_type + " Pressure and Darcy velocity")
    ax.set_xlabel(r"$x$ $[m]$", fontsize=20)
    ax.set_ylabel(r"$y$ $[m]$", fontsize=20, rotation=0)

    # Add a colourbar to the right of the plot, and make it the same height as the plot
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    cbar = plt.colorbar(
        im_c,
        ax=ax,
        cmap=colmap_P,
        norm=norm,
        cax=cax,
        format=FuncFormatter(cbar_fmt),
    )

    # # # # Add horizontal marker lines to the colourbar corresponding to the current contours
    # # # cbar.ax.axhline(max_val, c="r")
    # # # for _, L in enumerate(contour_lvls):
    # # #     cbar.ax.axhline(L, c="k", ls=":")


###################################################################################################
## Command-line arguments and file management #####################################################
###################################################################################################


# Initialize the argument parser
parser = argparse.ArgumentParser(description="Input and Output data folders")

# Add arguments for input and output folder with default values
parser.add_argument(
    "--input", type=str, default="./Input/", help="Input data folder path"
)
parser.add_argument(
    "--output", type=str, default="./Output/", help="Output data folder path"
)
parser.add_argument(
    "--batch",
    type=str,
    default=None,
    help="Path to a particular run from a batch. This overrides --input and --output",
)

# Parse the arguments
args = parser.parse_args()

# If the --batch option has been specifed, that should override
# the --input and --output options
if args.batch != None:
    # Extract the batch folder path, and remove any trailing quotation marks and slashes
    batch_folder = Path(args.batch)
    # Set the input and output paths relative to this
    input_folder = batch_folder / "Input"
    output_folder = batch_folder / "Output"
else:
    # Extract the input and output folder paths, and remove any trailing quotation marks and slashes
    input_folder = Path(args.input)
    output_folder = Path(args.output)


# Ensure that the folders exist
if not input_folder.is_dir():
    print(f"Error: The input folder {input_folder} does not exist.")
    exit(1)
if not output_folder.is_dir():
    print(f"Warning: The output folder {output_folder} does not exist.")
    exit(1)


print(f" Input folder : {input_folder}")
print(f" Output folder: {output_folder}")

###################################################################################################
## Load data ######################################################################################
###################################################################################################

# load plot times
plot_times = np.loadtxt(f"{output_folder}/Other/plot_times.txt")
n_plot = len(plot_times)

# load parameter values
parameters = np.loadtxt(
    f"{output_folder}/Other/parameter_values.txt"
)  # nx ny dx dy M Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve
nx, ny, nz, dx, dy, M, Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve = parameters
nx, ny, nz = [int(i) for i in [nx, ny, nz]]  # Convert these to integers


# Spatial grids
X_grid = np.arange(-(nx - 1) / 2.0, (nx - 1) / 2.0 + 1.0) * dx
Y_grid = np.arange(-(ny - 1) / 2.0, (ny - 1) / 2.0 + 1.0) * dy

X, Y = np.meshgrid(X_grid, Y_grid)


# load injection locations
inj_locs_data = np.loadtxt(f"{output_folder}/Other/injection_locations.txt")

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


# Mesh coarsening for vector flow field plots
flow_field_step = 5


# Reservoir Ceiling and Basement topography
H0 = np.loadtxt(f"{input_folder}/ceil_topo.txt")
B0 = np.loadtxt(f"{input_folder}/base_topo.txt")

D0 = B0 - H0

# # Porosity and Permeability fields
# Porosity = np.loadtxt(f"{input_folder}/porosity.txt")
# Permeability = np.loadtxt(f"{input_folder}/permeability.txt")

# Porosity and Permeability fields
Porosity_full = read_in_3D_array(f"{input_folder}/porosity.txt", nx, ny, nz)
Permeability_full = read_in_3D_array(f"{input_folder}/permeability.txt", nx, ny, nz)

Porosity = np.squeeze(Porosity_full[:, :, 0])
Permeability = np.squeeze(Permeability_full[:, :, 0])

Perm_int, Z_array = calculate_cumulative_integral(Permeability_full)


## Load plot data into arrays
# Free and trapped CO2 thicknesses
h_array = np.zeros([ny, nx, n_plot])
h_res_array = np.zeros([ny, nx, n_plot])

# Free and trapped CO2 local volumes
V_mobile_array = np.zeros([ny, nx, n_plot])
V_trapped_array = np.zeros([ny, nx, n_plot])

# Pressure fields for the CO2 and ambient
P_amb_array = np.zeros([ny, nx, n_plot])
P_cur_array = np.zeros([ny, nx, n_plot])

# Pressure gradients for the ambient fluid
dPdx_amb_array = np.zeros([ny, nx, n_plot])
dPdy_amb_array = np.zeros([ny, nx, n_plot])

# Pressure gradients for the CO2
dPdx_cur_array = np.zeros([ny, nx, n_plot])
dPdy_cur_array = np.zeros([ny, nx, n_plot])

# Darcy velocities for the ambient fluid
Ux_amb_array = np.zeros([ny, nx, n_plot])
Uy_amb_array = np.zeros([ny, nx, n_plot])


# Darcy velocity for the CO2
Ux_cur_array = np.zeros([ny, nx, n_plot])
Uy_cur_array = np.zeros([ny, nx, n_plot])

# Permeability integrals
perm_int_c = np.zeros([ny, nx])
perm_int_a = np.zeros([ny, nx])


# Load data from each of the output times and fill in the corresponding arrays
for i, t in enumerate(plot_times):

    h = np.loadtxt(
        f"{output_folder}/Current_Thickness/h" + "{0:02d}".format(i) + ".txt"
    )

    h_array[:, :, i] = h

    h_res_array[:, :, i] = np.loadtxt(
        f"{output_folder}/Current_Thickness/h_res" + "{0:02d}".format(i) + ".txt"
    )

    V_mobile_array[:, :, i] = np.loadtxt(
        f"{output_folder}/Current_Volume/V" + "{0:02d}".format(i) + ".txt"
    )

    V_trapped_array[:, :, i] = np.loadtxt(
        f"{output_folder}/Current_Volume/V_res" + "{0:02d}".format(i) + ".txt"
    )

    P_amb_array[:, :, i] = np.loadtxt(
        f"{output_folder}/Current_Pressure/P" + "{0:02d}".format(i) + ".txt"
    )

    P_cur_array[:, :, i] = P_amb_array[:, :, i] + Gamma_val * (
        H0[:, :] + h_array[:, :, i] - np.average(H0[:, :])
    )

    dPdx_amb_array[:, :, i], dPdy_amb_array[:, :, i] = np.gradient(
        P_amb_array[:, :, i], dx, dy
    )

    dPdx_cur_array[:, :, i], dPdy_cur_array[:, :, i] = np.gradient(
        P_cur_array[:, :, i], dx, dy
    )

    perm_int_c[:, :] = vertical_interpolation(Perm_int, h)
    perm_int_a[:, :] = Perm_int[:, :, -1] - perm_int_c[:, :]

    Ux_amb_array[:, :, i] = -M * perm_int_a[:, :] * dPdx_amb_array[:, :, i]
    Uy_amb_array[:, :, i] = -M * perm_int_a[:, :] * dPdy_amb_array[:, :, i]

    Ux_cur_array[:, :, i] = -perm_int_c[:, :] * dPdx_cur_array[:, :, i]
    Uy_cur_array[:, :, i] = -perm_int_c[:, :] * dPdy_cur_array[:, :, i]


# Load the volume data recorded at each step, rather than just at the output times
# Volume_data = np.loadtxt(f"{output_folder}/Other/Volumes.txt",skiprows=1)
# Iterations = Volume_data[:, 0]
# Times = Volume_data[:, 1]
# Volumes_mobile = Volume_data[:, 2]
# Volumes_trapped = Volume_data[:, 3]
# Volumes_injected = Volume_data[:, 4]
# # max_Vol = np.max(Volumes_mobile)
# max_Vol = np.max(Volumes_injected)

Volume_data = np.loadtxt(f"{output_folder}/Other/Volumes.txt", skiprows=1)
# Transpose Volume_data and assign the variables to each of the rows
(
    Iterations,
    Times,
    Volumes_injected,
    Volumes_total,
    Volumes_mobile,
    Volumes_trapped,
    Volumes_dissolution,
    Volumes_boundary,
    Volumes_mobile_local_max,
    Volumes_trapped_local_max,
) = Volume_data.T

Volume_profile_mobile = np.zeros([n_plot])
Volume_profile_trapped = np.zeros([n_plot])

for i, t in enumerate(plot_times):
    Volume_profile_mobile[i] = np.sum(V_mobile_array[:, :, i])
    Volume_profile_trapped[i] = np.sum(V_trapped_array[:, :, i])


# Time scale and units
if Times[-1] > 730.0:
    # Express in terms of standard years
    Time_scale = 365.25
    Time_unit = "years"
    Times = Times / Time_scale
    plot_times = plot_times / Time_scale
else:
    Time_scale = 1.0
    Time_unit = "days"


# Print appropriate maximum and minimum values
max_current_thickness = np.amax(h_array)
print("Maximum mobile current thickness is  " + str(max_current_thickness))

max_h_res = np.max(h_res_array)
print("Maximum trapped current thickness is " + str(max_h_res))

max_V_mobile = np.amax(Volumes_mobile)
print("Maximum mobile volume is             " + str(max_V_mobile))

max_V_trapped = np.amax(Volumes_trapped)
print("Maximum trapped volume is            " + str(max_V_trapped))

max_amb_pressure = np.max(P_amb_array)
print("Maximum ambient pressure is          " + str(max_amb_pressure))
min_amb_pressure = np.min(P_amb_array)
print("Minimum ambient pressure is          " + str(min_amb_pressure))

max_cur_pressure = np.max(P_cur_array)
print("Maximum current pressure is          " + str(max_cur_pressure))
min_cur_pressure = np.min(P_cur_array)
print("Minimum current pressure is          " + str(min_cur_pressure))

max_abs_pressure = np.max(np.abs(P_amb_array))


# Colour maps for the respective plots below
colmap_h = cm.BuPu
colmap_h_res = cm.Reds
colmap_V_mobile = cm.BuPu  # cm.BuGn
colmap_V_trapped = cm.Reds  # cm.OrRd
colmap_P = cm.seismic

colmap_perm = cm.gray
colmap_poro = cm.gray

# norms for the permeability and porosity (the rest are calculated within the loop below
# as they change between plots)
norm_perm = colors.Normalize(0.0, np.max(Permeability))
norm_poro = colors.Normalize(0.0, np.max(Porosity))


# How detailed to make the 3D plot of the current - number of steps between points to plot
stride_val = 1

# Threshold thickness value for masking and bounding contour - used in determining the
# edge of the current
h_threshold = 1e-4
h_res_threshold = 1e-4
V_threshold = 1e-4
V_res_threshold = 1e-4
P_threshold = 1e-4

# Remaining percentage contours to plot
h_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_current_thickness
cont_h_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]

h_res_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_h_res
cont_h_res_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]

V_mobile_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_V_mobile
cont_V_mobile_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]

V_trapped_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_V_trapped
cont_V_trapped_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]

# more levels for the Pressure as it can be negative
P_levels = (
    np.array([-0.95, -0.8, -0.6, -0.4, -0.2, -0.1, 0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95])
    * max_abs_pressure
)
cont_P_strs = [
    "-95%",
    "-80%",
    "-60%",
    "-40%",
    "-20%",
    "-10%",
    "0",
    "10%",
    "20%",
    "40%",
    "60%",
    "80%",
    "95%",
]


# Formatting string for the contour labels
cbar_fmt = lambda x, pos: "{:.3f}".format(x)

# Colours for the bounding contours of the free and trapped CO2 regions
clr_V = "green"
clr_V_res = "tab:brown"


## Main loop over each of the output plots
# Plot in reverse order, so we can see the final state sooner.
for i in range(n_plot - 1, -1, -1):
    t = plot_times[i]
    # If there are no contours (because the function is zero everywhere), supress the warning message
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="No contour levels were found within the data range."
        )

        print("Plot " + str(i))

        # load specific data for this frame
        h = h_array[:, :, i]
        h_res = h_res_array[:, :, i]
        V_mobile = V_mobile_array[:, :, i]
        V_trapped = V_trapped_array[:, :, i]
        P_amb = P_amb_array[:, :, i]
        P_cur = P_cur_array[:, :, i]
        Ux_amb = Ux_amb_array[:, :, i]
        Uy_amb = Uy_amb_array[:, :, i]
        Ux_cur = Ux_cur_array[:, :, i]
        Uy_cur = Uy_cur_array[:, :, i]

        fig = plt.figure(figsize=(15, 10))

        gs = GridSpec(3, 6, figure=fig)

        # str_title = "[{:02d} of {:02d}]: t = {:3.2e}".format(i, n_plot - 1, t)

        # fig.suptitle(
        #     str_title,
        #     fontsize=16,
        #     x=0.01,
        #     y=0.98,
        #     horizontalalignment="left",
        #     verticalalignment="top",
        #     bbox=dict(facecolor="none", edgecolor="black"),
        # )

        str_title = f"[{i:02d} of {n_plot - 1:02d}]: t = {t:.2f} {Time_unit}"

        fig.suptitle(
            str_title,
            fontsize=14,
            x=0.99,
            y=0.01,
            horizontalalignment="right",
            verticalalignment="bottom",
            bbox=dict(facecolor="none", edgecolor="black"),
        )


        ### Ceiling Topography ####################################################################
        ax_c = fig.add_subplot(gs[0, 0])

        cont_H0 = ax_c.contour(X_grid, Y_grid, H0, colors="grey", alpha=0.25)
        ax_c.clabel(cont_H0, inline=True, fmt="H0 = %.2g", fontsize=10)

        ax_c.set_title("Ceil. Topography")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

        ### Basement Topography ###################################################################
        ax_c = fig.add_subplot(gs[1, 0])

        cont_B0 = ax_c.contour(X_grid, Y_grid, B0, colors="grey", alpha=0.25)
        ax_c.clabel(cont_B0, inline=True, fmt="B0 = %.2g", fontsize=10)

        ax_c.set_title("Base. Topography")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

        ### Porosity ##############################################################################
        ax_c = fig.add_subplot(gs[0, 1])

        # surface plot of Porosity
        im_c = ax_c.imshow(
            np.flipud(Porosity),
            cmap=colmap_poro,
            norm=norm_poro,
            extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        )

        # Plot injection locations
        for k in range(0, n_inj_locs):
            ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

        ax_c.set_title("Porosity")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

        # Add a colourbar to the right of the plot, and make it the same height as the plot
        divider = make_axes_locatable(ax_c)
        cax_c = divider.append_axes("right", size="5%", pad=0.05)

        cbar_poro = plt.colorbar(
            im_c,
            ax=ax_c,
            cmap=colmap_poro,
            norm=norm_poro,
            cax=cax_c,
            format=FuncFormatter(cbar_fmt),
        )

        ### Permeability ########################################################################
        ax_c = fig.add_subplot(gs[1, 1])

        # surface plot of Permeability
        im_c = ax_c.imshow(
            np.flipud(Permeability),
            cmap=colmap_perm,
            norm=norm_perm,
            extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        )

        # Plot injection locations
        for k in range(0, n_inj_locs):
            ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

        ax_c.set_title("Permeability")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

        # Add a colourbar to the right of the plot, and make it the same height as the plot
        divider = make_axes_locatable(ax_c)
        cax_c = divider.append_axes("right", size="5%", pad=0.05)

        cbar_perm = plt.colorbar(
            im_c,
            ax=ax_c,
            cmap=colmap_perm,
            norm=norm_perm,
            cax=cax_c,
            format=FuncFormatter(cbar_fmt),
        )

        ### Current Pressure #####################################################################
        ax_c = fig.add_subplot(gs[0:, 2:4])

        generate_pressure_plot(
            P_cur,
            Ux_cur,
            Uy_cur,
            V_mobile,
            V_trapped,
            colmap_P,
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            "Current",
            ax_c,
        )

        ### Ambient Pressure #####################################################################
        ax_c = fig.add_subplot(gs[0:, 4:])

        generate_pressure_plot(
            P_amb,
            Ux_amb,
            Uy_amb,
            V_mobile,
            V_trapped,
            colmap_P,
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            "Ambient",
            ax_c,
        )

        ### Total Volume ##############################################################################

        ax_c = fig.add_subplot(gs[2, 0:2])

        # # Total mobile volume - output times
        # ax_c.plot(plot_times, Volume_profile_mobile, ".", markersize=8, color=clr_V)
        # # Total trapped volume - output times
        # ax_c.plot(
        #     plot_times, Volume_profile_trapped, ".", markersize=8, color=clr_V_res
        # )
        # # Total mobile volume - every step
        # ax_c.plot(Times, Volumes_mobile, "-.", color=clr_V, label="Mobile CO$_{2}$")
        # # Total trapped volume - every step
        # ax_c.plot(
        #     Times, Volumes_trapped, ":", color=clr_V_res, label="Trapped CO$_{2}$"
        # )
        # # Total CO2 phase - every step
        # ax_c.plot(
        #     Times, Volumes_mobile + Volumes_trapped, "k--", label="Total CO$_{2}$ phase"
        # )
        # # Total CO2 injected - every step
        # ax_c.plot(Times, Volumes_injected, "b-", label="Total CO$_{2}$ injected")

        # # Current output time for this plot
        # ax_c.axvline(x=t, color="gray")

        # ax_c.set_title("Volume")
        # ax_c.set_xlabel(r"$t$", fontsize=20)
        # ax_c.set_ylabel(r"$V$", fontsize=20, labelpad=20, rotation=0)

        # ax_c.legend(loc="upper left", fontsize=9)

        # # Add a grid to the volume plot
        # ax_c.grid(True, linestyle="-", alpha=0.5)


        # Total CO2 injected - every step
        ax_c.plot(Times, Volumes_injected, "b-", label="CO$_{2}$ injected")

        # Total active volume - output times
        ax_c.plot(plot_times, Volume_profile_mobile, ".", markersize=8, color=clr_V)
        # Total mobile volume - every step
        ax_c.plot(Times, Volumes_mobile, "-.", color=clr_V, label="Mobile")

        # Total trapped volume - output times
        ax_c.plot(
            plot_times, Volume_profile_trapped, ".", markersize=8, color=clr_V_res
        )
        # Total trapped volume - every step
        ax_c.plot(Times, Volumes_trapped, ":", color=clr_V_res, label="Trapped")

        # Total dissolved volume - every step
        ax_c.plot(Times, Volumes_dissolution, "m-", label="Dissolved")

        # Total volume lost through the boundaries - every step
        ax_c.plot(
            Times,
            Volumes_boundary,
            "c--",
            label="Boundaries",
        )

        # Total free-phase CO2 - every step
        ax_c.plot(
            Times,
            Volumes_mobile + Volumes_trapped,
            "k--",
            label="(Mo. + Tr.)",
        )

        # Total CO2 phase (including dissolved) - every step
        ax_c.plot(
            Times,
            Volumes_total,
            "r:",
            label="Total",
        )

        # Current output time for this plot
        ax_c.axvline(x=t, color="gray")

        # ax_c.set_title("Volume")
        ax_c.set_xlabel(f"$t$ [{Time_unit}]", fontsize=16)
        ax_c.set_ylabel(r"Volume [$m^3$]", fontsize=16, labelpad=20)

        # ax_c.legend(loc="upper left", fontsize=6)
        ax_c.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=6)

        # Add a grid to the volume plot
        ax_c.grid(True, linestyle="-", alpha=0.5)






        # Tidy up the layout of subplots
        # plt.tight_layout(rect=[0, 0, 1, 0.95])
        gs.tight_layout(figure=fig)

        # Save this frame
        plt.savefig("./plots/Overview_plot/temp/t" + "{0:02d}".format(i) + ".png")
        plt.close()
