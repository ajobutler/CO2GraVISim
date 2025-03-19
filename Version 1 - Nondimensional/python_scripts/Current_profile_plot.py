# Current_profile_plot.py
# This script produces plots from the output of CO2GraVISim,
# showing
# - the ceiling topography,
# - the local volumes of the mobile and trapped CO2
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

import argparse
import os

### Deal with inputs for input and output data folders ######################################

# Initialize the argument parser
parser = argparse.ArgumentParser(description='Input and Output data folders')

# Add arguments for input and output folder with default values
parser.add_argument('-input', type=str, default='./Input/',
                    help='Input data folder path')
parser.add_argument('-output', type=str, default='./Output/',
                    help='Output data folder path')

# Parse the arguments
args = parser.parse_args()

# Extract the input and output folder paths, and remove any trailing slashes
input_folder = args.input.rstrip('/\\')
output_folder = args.output.rstrip('/\\')


# Ensure that the folders exist
if not os.path.isdir(input_folder):
    print(f"Error: The input folder {input_folder} does not exist.")
    exit(1)
if not os.path.isdir(output_folder):
    print(f"Warning: The output folder {output_folder} does not exist.")
    exit(1)


print(f" Input folder : {input_folder}")
print(f" Output folder: {output_folder}")

#############################################################################################


# load plot times
plot_times = np.loadtxt(f"{output_folder}/Other/plot_times.txt")


# load parameter values
parameters = np.loadtxt(
    f"{output_folder}/Other/parameter_values.txt"
)  # nx ny dx dy M Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve
nx = int(parameters[0])
ny = int(parameters[1])
dx = parameters[2]
dy = parameters[3]

M           = parameters[4]
Gamma_val   = parameters[5]
s_c_r       = parameters[6]
s_a_i       = parameters[7]
C_sat       = parameters[8]
q_dissolve  = parameters[9]


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

# Porosity and Permeability fields
Porosity = np.loadtxt(f"{input_folder}/porosity.txt")
Permeability = np.loadtxt(f"{input_folder}/permeability.txt")


## Load plot data into arrays
# Free and trapped CO2 thicknesses
h_array     = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])
h_res_array = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])

# Free and trapped CO2 local volumes
V_active_array  = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])
V_trapped_array = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])



# Load data from each of the output times and fill in the corresponding arrays
for i, t in enumerate(plot_times):
    # load height data for this frame
    if i == 0:
        h_array[:, :, i] = np.zeros([np.shape(H0)[0], np.shape(H0)[1]])
    else:
        h_array[:, :, i] = np.loadtxt(
            f"{output_folder}/Current_Thickness/h" + "{0:02d}".format(i) + ".txt"
        )

    h_res_array[:, :, i] = np.loadtxt(
        f"{output_folder}/Current_Thickness/h_res" + "{0:02d}".format(i) + ".txt"
    )


    V_active_array[:, :, i] = Porosity * (1.0 - s_a_i) * h_array[:, :, i]
    V_trapped_array[:, :, i] = Porosity * s_c_r * h_res_array[:, :, i]


# Print appropriate maximum and minimum values
max_current_thickness = np.amax(h_array)
print("Maximum active current thickness is  " + str(max_current_thickness))

max_h_res = np.max(h_res_array)
print("Maximum trapped current thickness is " + str(max_h_res))

max_V_active = np.amax(V_active_array)
print("Maximum active volume is             " + str(max_V_active))

max_V_trapped = np.amax(V_trapped_array)
print("Maximum trapped volume is            " + str(max_V_trapped))



# Load the volume data recorded at each step, rather than just at the output times
Volume_data = np.loadtxt(f"{output_folder}/Other/Volumes.txt")
Times = Volume_data[:, 0]
Volumes_active = Volume_data[:, 1]
Volumes_trapped = Volume_data[:, 2]
Volumes_injected = Volume_data[:, 3]

# max_Vol = np.max(Volumes_active)
max_Vol = np.max(Volumes_injected)

Volume_profile_active = np.zeros([len(plot_times)])
Volume_profile_trapped = np.zeros([len(plot_times)])

for i, t in enumerate(plot_times):
    Volume_profile_active[i] = np.sum(V_active_array[:, :, i]) * dx * dy
    Volume_profile_trapped[i] = np.sum(V_trapped_array[:, :, i]) * dx * dy


# Colour maps for the respective plots below
colmap_h = cm.BuPu
colmap_h_res = cm.Reds
colmap_V_active = cm.BuPu  # cm.BuGn
colmap_V_trapped = cm.Reds  # cm.OrRd
colmap_P = cm.seismic

colmap_perm = cm.gray
colmap_poro = cm.gray

# norms for the permeability and porosity (the rest are calculated within the loop below
# as they change between plots)
norm_perm = colors.Normalize(0., np.max(Permeability))
norm_poro = colors.Normalize(0., np.max(Porosity))


# How detailed to make the 3D plot of the current - number of steps between points to plot
stride_val = 1

# Threshold thickness value for masking and bounding contour - used in determining the
# edge of the current
h_threshold = 1e-4
h_res_threshold = 1e-4
V_active_threshold = 1e-4
V_trapped_threshold = 1e-4
P_threshold = 1e-4

# Remaining percentage contours to plot
h_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_current_thickness
cont_h_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]

h_res_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_h_res
cont_h_res_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]

V_active_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_V_active
cont_V_active_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]

V_trapped_levels = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.95]) * max_V_trapped
cont_V_trapped_strs = ["10%", "20%", "40%", "60%", "80%", "95%"]



# Formatting string for the contour labels
cbar_fmt = lambda x, pos: "{:.3f}".format(x)

# Colours for the bounding contours of the free and trapped CO2 regions
clr_h = "green"
clr_h_res = "tab:brown"


## Main loop over each of the output plots
# Plot in reverse order, so we can see the final state sooner.
for i in range(len(plot_times)-1,-1,-1):
    t = plot_times[i]
    # If there are no contours (because the function is zero everywhere), supress the warning message
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="No contour levels were found within the data range."
        )

        print("Plot " + str(i))

        # load height data for this frame
        h = h_array[:, :, i]
        h_res = h_res_array[:, :, i]
        V_active = V_active_array[:, :, i]
        V_trapped = V_trapped_array[:, :, i]

        # create colormap according to h value, masked with a threshold value so that it doesn't colour the h=0 regions
        color_dimension_h = np.ma.masked_where(h < h_threshold, h)
        minn, maxx = 0.0, max_current_thickness
        norm_h = colors.Normalize(minn, maxx)
        m = plt.cm.ScalarMappable(norm=norm_h, cmap=colmap_h)
        m.set_array([])
        fcolors_h = m.to_rgba(color_dimension_h)

        # h_res colormap
        color_dimension_h_res = np.ma.masked_where(h_res < h_res_threshold, h_res)
        minn, maxx = 0.0, max_h_res  # color_dimension.max()
        norm_h_res = colors.Normalize(minn, maxx)
        m = plt.cm.ScalarMappable(norm=norm_h_res, cmap=colmap_h_res)
        m.set_array([])
        fcolors_h_res = m.to_rgba(color_dimension_h_res)

        # V_active colormap
        color_dimension_V_active = np.ma.masked_where(
            V_active < V_active_threshold, V_active
        )
        minn, maxx = 0.0, max_V_active
        norm_V_active = colors.Normalize(minn, maxx)
        m = plt.cm.ScalarMappable(norm=norm_V_active, cmap=colmap_V_active)
        m.set_array([])
        fcolors_V_active = m.to_rgba(color_dimension_V_active)

        # V_trapped colormap
        color_dimension_V_trapped = np.ma.masked_where(
            V_trapped < V_trapped_threshold, V_trapped
        )
        minn, maxx = 0.0, max_V_trapped
        norm_V_trapped = colors.Normalize(minn, maxx)
        m = plt.cm.ScalarMappable(norm=norm_V_trapped, cmap=colmap_V_trapped)
        m.set_array([])
        fcolors_V_trapped = m.to_rgba(color_dimension_V_trapped)


        fig = plt.figure(figsize=(15, 10))

        gs = GridSpec(6, 3, figure=fig)

        str_title = "[{:02d} of {:02d}]: t = {:3.2e}".format(i, len(plot_times) - 1, t)


        fig.suptitle(
            str_title,
            fontsize=16,
            x=0.01,
            y=0.98,
            horizontalalignment="left",
            verticalalignment="top",
            bbox=dict(facecolor="none", edgecolor="black"),
        )

        ### Porosity ############################################################################
        ax_c = fig.add_subplot(gs[0:2, 0])

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
        ax_c = fig.add_subplot(gs[2:4, 0])

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


        ### Current Volume - active ###################################################################
        ax_c = fig.add_subplot(gs[0:3, 1:])

        # surface plot of active current volume
        im_c = ax_c.imshow(
            np.flipud(color_dimension_V_active),
            cmap=colmap_V_active,
            norm=norm_V_active,
            extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        )
        # contours of ceiling topography
        cont_H0 = ax_c.contour(X_grid, Y_grid, H0, colors="gray", alpha=0.25)
        # edge contour of active current
        cont_hb = ax_c.contour(
            X_grid,
            Y_grid,
            h * (h >= h_threshold),
            colors=clr_h,
            alpha=1.0,
            levels=[h_threshold],
            linewidths=2,
        )
        # contours of active current volume
        cont_V_active = ax_c.contour(
            X_grid,
            Y_grid,
            V_active * (h >= h_threshold),
            colors="black",
            alpha=0.5,
            levels=V_active_levels,
        )
        # contour labels
        cont_V_active_fmt = {}
        for l, s in zip(cont_V_active.levels, cont_V_active_strs):
            cont_V_active_fmt[l] = s
        ax_c.clabel(cont_H0, inline=True, fmt="H0 = %.2g", fontsize=10)
        ax_c.clabel(cont_V_active, inline=True, fmt=cont_V_active_fmt, fontsize=10)

        # Plot injection locations
        for k in range(0, n_inj_locs):
            ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")
        
        ax_c.set_title("Current Volume (active)")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

        # Add a colourbar to the right of the plot, and make it the same height as the plot
        divider = make_axes_locatable(ax_c)
        cax_c = divider.append_axes("right", size="5%", pad=0.05)

        cbar_Va = plt.colorbar(
            im_c,
            ax=ax_c,
            cmap=colmap_V_active,
            norm=norm_V_active,
            cax=cax_c,
            format=FuncFormatter(cbar_fmt),
        )

        ### Current Volume - trapped ##############################################################
        ax_c = fig.add_subplot(gs[3:, 1:])

        # surface plot of trapped current volume
        im_c = ax_c.imshow(
            np.flipud(color_dimension_V_trapped),
            cmap=colmap_V_trapped,
            norm=norm_V_trapped,
            extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        )
        # contours of ceiling topography
        cont_H0 = ax_c.contour(X_grid, Y_grid, H0, colors="gray", alpha=0.25)
        # edge contour of trapped current
        cont_hresb = ax_c.contour(
            X_grid,
            Y_grid,
            h_res * (h_res >= h_res_threshold),
            colors=clr_h_res,
            alpha=1.0,
            levels=[h_res_threshold],
            linewidths=2,
        )
        # contours of trapped current volume
        cont_V_trapped = ax_c.contour(
            X_grid,
            Y_grid,
            V_trapped * (h_res >= h_res_threshold),
            colors="black",
            alpha=0.5,
            levels=V_trapped_levels,
        )
        # contour labels
        cont_V_trapped_fmt = {}
        for l, s in zip(cont_V_trapped.levels, cont_V_trapped_strs):
            cont_V_trapped_fmt[l] = s
        ax_c.clabel(cont_H0, inline=True, fmt="H0 = %.2g", fontsize=10)
        ax_c.clabel(cont_V_trapped, inline=True, fmt=cont_V_trapped_fmt, fontsize=10)

        # Plot injection locations
        for k in range(0, n_inj_locs):
            ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

        ax_c.set_title("Current Volume (trapped)")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

        # Add a colourbar to the right of the plot, and make it the same height as the plot
        divider = make_axes_locatable(ax_c)
        cax_c = divider.append_axes("right", size="5%", pad=0.05)

        cbar_Vt = plt.colorbar(
            im_c,
            ax=ax_c,
            cmap=colmap_V_trapped,
            norm=norm_V_trapped,
            cax=cax_c,
            format=FuncFormatter(cbar_fmt),
        )


        ### Total Volume ##############################################################################
 
        ax_c = fig.add_subplot(gs[4:, 0])

        # Total active volume - output times
        ax_c.plot(plot_times, Volume_profile_active, ".", markersize=8, color=clr_h)
        # Total trapped volume - output times
        ax_c.plot(plot_times, Volume_profile_trapped, ".", markersize=8, color=clr_h_res)
        # Total active volume - every step
        ax_c.plot(Times, Volumes_active, "-.", color=clr_h, label="Mobile CO$_{2}$")
        # Total trapped volume - every step
        ax_c.plot(Times, Volumes_trapped, ":", color=clr_h_res, label="Trapped CO$_{2}$")
        # Total CO2 phase - every step
        ax_c.plot(Times, Volumes_active + Volumes_trapped, "k--", label="Total CO$_{2}$ phase")
        # Total CO2 injected - every step
        ax_c.plot(Times, Volumes_injected, "b-", label="Total CO$_{2}$ injected")


        # Current output time for this plot
        ax_c.axvline(x=t, color="gray")

        ax_c.set_title("Volume")
        ax_c.set_xlabel(r"$t$", fontsize=20)
        ax_c.set_ylabel(r"$V$", fontsize=20, labelpad=20, rotation=0)

        ax_c.legend(loc="upper left", fontsize=9)

        #Add a grid to the volume plot
        ax_c.grid(True, linestyle="-", alpha=0.5)

        #Tidy up the layout of subplots
        plt.tight_layout(rect=[0, 0, 1, 0.95])

        # Save this frame
        plt.savefig("./plots/Current_profile_plot/temp/t" + "{0:02d}".format(i) + ".png")
        plt.close()
