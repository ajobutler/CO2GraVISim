# Current_profile_plot.py
# This script produces plots from the output of CO2GraVISim,
# showing
# - the ceiling topography,
# - the local volumes of the mobile and trapped CO2
# - the total volumes of free and trapped CO2 over time

"""
TODO:- Fix my current bodge of adding in a blank column to get the spacing to work.
        Tight_layout doesn't get on well with this.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec
import warnings

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import argparse
import os
from pathlib import Path

# from multiprocessing import Pool
# from multiprocessing.pool import ThreadPool

from timeit import default_timer as timer

###################################################################################################
## Functions ######################################################################################
###################################################################################################


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


def colourmap_arrays(Array, colmap, minn, maxx):
    # This function produces various colourmap arrays used for the 2D and
    # 3D plots generated here.

    norm = colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap=colmap)
    m.set_array([])
    fcolors = m.to_rgba(Array)

    return fcolors, norm, m


def evenly_spaced_layers(N, n_layers):
    # This function produces a list of evenly spaced integers betwen 0 and N-1
    # that includes these end points. This is used for plotting layers within
    # the reservoir geometry

    if n_layers <= 1:
        return [0, N - 1]

    indices = np.round(np.linspace(0, N - 1, n_layers)).astype(int)
    # In case N-1 isn't included in the list, make sure it is
    indices[-1] = N - 1

    return indices


def generate_3D_animated_plot(
    Array,
    colmap,
    extrema_vals,
    X,
    Y,
    H0,
    B0,
    z,
    n_inj_locs,
    inj_loc_idxs,
    inj_grid_vals,
    str_title,
    elev_val,
    azim_val,
    Rotate_flag,
    n_plot,
    ax,
    cbax,
    fig,
    i,
):
    # This function produces a 3D plot of the reservoir topography: the ceiling layer, H0,
    # the basement layer, B0, and a series of interpolated layers in between. These layers
    # are coloured based on the values in a 3D array (Array) - either Porosity or
    # Permeability.

    min_val, max_val = extrema_vals

    # Formatting string for the contour labels
    cbar_fmt = lambda x, pos: "{:.3f}".format(x)

    n_layer_indices = len(layer_indices)
    line_handles = []
    surface_handles = []

    for i, k in enumerate(layer_indices):
        fcolors, _, m = colourmap_arrays(Array[:, :, k], colmap, min_val, max_val)

        z_order_val = 2 * (n_layer_indices - i)

        # Plot the topography of this layer, and colour it appropriately
        surface = ax.plot_surface(
            X,
            Y,
            H0 + z[k] * (B0 - H0),
            # edgecolor="none",
            facecolors=fcolors,  # np.transpose(fcolors, (1, 0, 2)),
            vmin=min_val,
            vmax=max_val,
            shade=False,
            # alpha=0.25,
            zorder=z_order_val,
        )
        surface_handles.append(surface)

    for i, _ in enumerate(layer_indices[:-1]):

        # z_order_val = 20
        z_order_val = 2 * (n_layer_indices - i) - 1

        # Plot injection locations
        if n_inj_locs == 1:
            idx_x, idx_y = map(int, inj_loc_idxs)
            k_u = layer_indices[i]
            k_l = layer_indices[i + 1]
            z_val_u = H0[idx_y, idx_x] + z[k_u] * (B0[idx_y, idx_x] - H0[idx_y, idx_x])
            z_val_l = H0[idx_y, idx_x] + z[k_l] * (B0[idx_y, idx_x] - H0[idx_y, idx_x])
            line = ax.plot(
                [inj_grid_vals[0, 0], inj_grid_vals[0, 0]],
                [inj_grid_vals[0, 1], inj_grid_vals[0, 1]],
                [z_val_u, z_val_l],
                "m",
                zorder=z_order_val,
            )
        else:
            for l in range(0, n_inj_locs):

                idx_x, idx_y = map(int, inj_loc_idxs[l, :])
                k_u = layer_indices[i]
                k_l = layer_indices[i + 1]
                z_val_u = H0[idx_y, idx_x] + z[k_u] * (
                    B0[idx_y, idx_x] - H0[idx_y, idx_x]
                )
                z_val_l = H0[idx_y, idx_x] + z[k_l] * (
                    B0[idx_y, idx_x] - H0[idx_y, idx_x]
                )

                line = ax.plot(
                    [inj_grid_vals[l, 0], inj_grid_vals[l, 0]],
                    [inj_grid_vals[l, 1], inj_grid_vals[l, 1]],
                    [z_val_u, z_val_l],
                    "m",
                    zorder=z_order_val,
                )
        # if n_inj_locs == 1:
        #     idx_x, idx_y = map(int, inj_loc_idxs)
        #     k_u = layer_indices[i]
        #     k_l = layer_indices[i + 1]
        #     z_val_u = H0[idx_x, idx_y] + z[k_u] * (B0[idx_x, idx_y] - H0[idx_x, idx_y])
        #     z_val_l = H0[idx_x, idx_y] + z[k_l] * (B0[idx_x, idx_y] - H0[idx_x, idx_y])
        #     line = ax.plot(
        #         [inj_grid_vals[0, 0], inj_grid_vals[0, 0]],
        #         [inj_grid_vals[0, 1], inj_grid_vals[0, 1]],
        #         [z_val_u, z_val_l],
        #         "m",
        #         zorder=z_order_val,
        #     )
        # else:
        #     for l in range(0, n_inj_locs):

        #         idx_x, idx_y = map(int, inj_loc_idxs[l, :])
        #         k_u = layer_indices[i]
        #         k_l = layer_indices[i + 1]
        #         z_val_u = H0[idx_x, idx_y] + z[k_u] * (
        #             B0[idx_x, idx_y] - H0[idx_x, idx_y]
        #         )
        #         z_val_l = H0[idx_x, idx_y] + z[k_l] * (
        #             B0[idx_x, idx_y] - H0[idx_x, idx_y]
        #         )

        #         line = ax.plot(
        #             [inj_grid_vals[l, 0] , inj_grid_vals[l, 0]],
        #             [inj_grid_vals[l, 1] , inj_grid_vals[l, 1]],
        #             [z_val_u, z_val_l],
        #             "m",
        #             zorder=z_order_val,
        #         )
        line_handles.append(line)

    # Since z points downwards
    ax.invert_zaxis()

    ax.set_title(str_title)
    ax.set_xlabel(r"$x$", fontsize=16)
    ax.set_ylabel(r"$y$", fontsize=16, rotation=0)
    ax.set_zlabel(r"$z$", fontsize=16, rotation=0)

    # Set the camera position, which moves to make the plot rotate if Rotate_flag = True
    print(f"{Rotate_flag = }")
    ax.view_init(
        elev=elev_val,
        azim=azim_val + Rotate_flag * 360 * (i / (n_plot - 1)),
        roll=None,
        vertical_axis="z",
        share=False,
    )

    # Add a colourbar
    fig.colorbar(
        m, ax=ax, cax=cbax, orientation="vertical", format=FuncFormatter(cbar_fmt)
    )


def generate_2D_plot(
    Array,
    threshold,
    extrema_vals,
    colmap,
    contour_clrs,
    contour_lvls,
    contour_strs,
    X_grid,
    Y_grid,
    H0,
    n_inj_locs,
    inj_grid_vals,
    str_title,
    ax,
    cbax,
):
    # This function produces a 2D plot of an array (e.g. mobile CO2 volume), along with relevant
    # contours for that array and for the ceiling topography

    min_val, max_val = extrema_vals

    # Formatting string for the contour labels
    cbar_fmt = lambda x, pos: "{:.3e}".format(x)

    # Masked array
    masked_array = np.ma.masked_where(Array < threshold, Array)
    _, norm, _ = colourmap_arrays(masked_array, colmap, min_val, max_val)

    # surface plot
    im_c = ax.imshow(
        np.flipud(masked_array),
        cmap=colmap,
        norm=norm,
        extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
    )
    # contours of ceiling topography
    cont_H0 = ax.contour(X_grid, Y_grid, H0, colors="grey", alpha=0.25)
    H0_labels = ax.clabel(cont_H0, inline=True, fmt="%.2f", fontsize=8)
    # H0_labels = ax.clabel(cont_H0, inline=True, fmt="H0 = %.2f", fontsize=6)

    # # Customize the label background
    # for label in H0_labels:
    #     label.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # edge contour of active current
    cont_hb = ax.contour(
        X_grid,
        Y_grid,
        masked_array,  # Array * (Array >= threshold),
        colors=contour_clrs,
        alpha=1.0,
        levels=[threshold],
        linewidths=2,
    )

    if max_val > threshold:
        # contours of active current volume
        cont_plot = ax.contour(
            X_grid,
            Y_grid,
            masked_array,  # Array * (Array >= threshold),
            colors="black",
            alpha=0.5,
            levels=contour_lvls,
        )
        # contour labels
        cont_fmt = {}
        for l, s in zip(cont_plot.levels, contour_strs):
            cont_fmt[l] = s

        ax.clabel(cont_plot, inline=True, fmt=cont_fmt, fontsize=10)

    # Plot injection locations
    for k in range(0, n_inj_locs):
        ax.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

    ax.set_title(str_title)
    ax.set_xlabel(r"$x$", fontsize=20)
    ax.set_ylabel(r"$y$", fontsize=20, rotation=0)

    cbar = plt.colorbar(
        im_c,
        ax=ax,
        cmap=colmap,
        norm=norm,
        cax=cbax,
        format=FuncFormatter(cbar_fmt),
    )

    # Add horizontal marker lines to the colourbar corresponding to the current contours
    cbar.ax.axhline(max_val, c="r")
    for _, L in enumerate(contour_lvls):
        cbar.ax.axhline(L, c="k", ls=":")


def load_common_files(filepath_A, filepath_B, filename, skiprows_val_A, skiprows_val_B):

    str_A = f"{filepath_A}/{filename}.txt"

    str_B = f"{filepath_B}/{filename}.txt"

    try:
        data = np.loadtxt(f"{filepath_A}/{filename}.txt", skiprows=skiprows_val_A)
    except FileNotFoundError:
        # File isn't in the filepath_A folder - look in filepath_B instead
        try:
            data = np.loadtxt(f"{filepath_B}/{filename}.txt", skiprows=skiprows_val_B)
        except FileNotFoundError:
            print(f"{filename} was not found in either {filepath_A} or {filepath_B}.")

    return data


def single_plot(i, fig, axs, cbaxs):
    # This function produces a single plot instance for a given value of i

    # Clear plot axes
    axs["poro_3D"].clear()
    axs["perm_3D"].clear()
    axs["active_volume"].clear()
    axs["trapped_volume"].clear()
    axs["volumes_lineplot"].clear()

    # Clear colourbar axes
    cbaxs["poro_3D"].clear()
    cbaxs["perm_3D"].clear()
    cbaxs["active_volume"].clear()
    cbaxs["trapped_volume"].clear()

    t = plot_times[i]

    # If there are no contours (because the function is zero everywhere), supress the warning message
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="No contour levels were found within the data range."
        )

        print("Plot " + str(i))

        print(" - Building Figure")

        # load height data for this frame
        V_mobile = V_mobile_array[:, :, i]
        V_trapped = V_trapped_array[:, :, i]

        ##### Figure setup #######################################################################

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

        ##### Porosity ############################################################################

        str_title = "Porosity"

        axs["poro_3D"].computed_zorder = False
        # This is required to make the zorder that I impose work correctly! (at least when viewed from above)
        # https://stackoverflow.com/questions/23188561/matplotlib-3d-plot-zorder-issue

        generate_3D_animated_plot(
            Porosity,
            colmap_poro,
            [0.0, 1.0],
            X,
            Y,
            H0,
            B0,
            z,
            n_inj_locs,
            inj_loc_idxs,
            inj_grid_vals,
            str_title,
            elev_val,
            azim_val,
            Rotate_flag,
            n_plot,
            axs["poro_3D"],
            cbaxs["poro_3D"],
            fig,
            i,
        )

        # # print("\n\nArtist features:")
        # # # Get all the artists in the plot and print their zorder
        # # for artist in axs["poro_3D"].get_children():
        # #     # Only print zorder for plot-related objects like lines, surfaces, etc.
        # #     if hasattr(artist, 'get_zorder'):
        # #         print(f"{artist.__class__.__name__} zorder: {artist.get_zorder()}")

        ##### Permeability ########################################################################

        str_title = "$\\log_{10}(Permeability)$"

        axs["perm_3D"].computed_zorder = False
        # This is required to make the zorder that I impose work correctly! (at least when viewed from above)
        # https://stackoverflow.com/questions/23188561/matplotlib-3d-plot-zorder-issue

        generate_3D_animated_plot(
            np.log10(Permeability),
            colmap_perm,
            [Log10_perm_min, Log10_perm_max],
            X,
            Y,
            H0,
            B0,
            z,
            n_inj_locs,
            inj_loc_idxs,
            inj_grid_vals,
            str_title,
            elev_val,
            azim_val,
            Rotate_flag,
            n_plot,
            axs["perm_3D"],
            cbaxs["perm_3D"],
            fig,
            i,
        )

        ##### Current Volume - active #############################################################

        V_act_max = np.max(V_mobile)
        V_mobile_levels = V_mobile_percentages * V_act_max

        str_title = "Current Volume (active)"
        generate_2D_plot(
            V_mobile,
            V_mobile_threshold,
            [0.0, V_act_max],
            colmap_V_mobile,
            clr_V,
            V_mobile_levels,
            cont_V_mobile_strs,
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            axs["active_volume"],
            cbaxs["active_volume"],
        )

        ##### Current Volume - trapped ############################################################

        V_trp_max = np.max(V_trapped)
        V_trapped_levels = V_trapped_percentages * V_trp_max

        str_title = "Current Volume (trapped)"
        generate_2D_plot(
            V_trapped,
            V_trapped_threshold,
            [0.0, V_trp_max],
            colmap_V_trapped,
            clr_V_res,
            V_trapped_levels,
            cont_V_trapped_strs,
            X_grid,
            Y_grid,
            H0,
            n_inj_locs,
            inj_grid_vals,
            str_title,
            axs["trapped_volume"],
            cbaxs["trapped_volume"],
        )

        ##### Total Volume ########################################################################

        ax_c = axs["volumes_lineplot"]

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

        ##### Save plot ###########################################################################

        # Tidy up the layout of subplots
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=2.0)
        # plt.tight_layout(rect=[0, 0, 1, 0.95])
        # gs.tight_layout(fig)
        # plt.tight_layout

        # plt.show()

        # Save this frame
        print(" - Saving Figure")
        plt.savefig(
            # "./plots/Current_profile_plot/temp/t" + "{0:02d}".format(i) + ".png",
            save_folder
            / f"t{i:02d}.png"
        )

        # plt.close()

    return True


def build_figure():
    fig, axs = plt.subplot_mosaic(
        [
            ["poro_3D", "active_volume", "trapped_volume"],
            ["perm_3D", "active_volume", "trapped_volume"],
            ["volumes_lineplot", "active_volume", "trapped_volume"],
        ],
        width_ratios=[0.75, 1.0, 1.0],
        # layout="constrained",
        # gridspec_kw={'hspace': 0.5, 'wspace': 0.75},
        # figsize=(15, 8),
        figsize=(12, 6),
    )
    # fig.set_constrained_layout_pads(w_pad=0.3, h_pad=0.1)#,
    #             # hspace=0.5, wspace=0.75)

    # Make porosity and permeability axes 3D
    axs["poro_3D"].remove()
    axs["poro_3D"] = plt.subplot(3, 3, 1, projection="3d")
    axs["perm_3D"].remove()
    axs["perm_3D"] = plt.subplot(3, 3, 4, projection="3d")

    # Associated colourbars
    cbaxs = {
        "poro_3D": None,
        "perm_3D": None,
        "active_volume": None,
        "trapped_volume": None,
    }

    # Add a colourbar to the right of the plot, and make it the same height as the plot
    cbaxs["poro_3D"] = inset_axes(
        axs["poro_3D"], width="5%", height="100%", loc="right", borderpad=-6
    )

    cbaxs["perm_3D"] = inset_axes(
        axs["perm_3D"], width="5%", height="100%", loc="right", borderpad=-6
    )

    divider = make_axes_locatable(axs["active_volume"])
    cbaxs["active_volume"] = divider.append_axes("right", size="5%", pad=0.05)

    divider = make_axes_locatable(axs["trapped_volume"])
    cbaxs["trapped_volume"] = divider.append_axes("right", size="5%", pad=0.05)

    return fig, axs, cbaxs


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
parser.add_argument(
    "--save_loc",
    type=str,
    default="./plots/Current_profile_plot/temp/",
    help="Folder path where the plots will be saved",
)
parser.add_argument(
    "--final", action="store_true", help="Activate to only produce the final plot"
)
parser.add_argument(
    "--view",
    type=float,
    nargs="+",
    default=[10.0, -35.0],
    help="Elevation and Azimuth for initial 3D camera position",
)
parser.add_argument(
    "--fixed", action="store_true", help="Activate to prevent 3D plots rotating"
)

# Parse the arguments
args = parser.parse_args()

# If the --batch option has been specifed, that should override
# the --input and --output options
if args.batch != None:
    # Extract the batch folder path, and remove any trailing quotation marks and slashes
    batch_run_folder = Path(args.batch)
    # Set the input and output paths relative to this
    input_folder = batch_run_folder / "Input"
    output_folder = batch_run_folder / "Output"
    common_files_folder = batch_run_folder / "../Common_files"
else:
    # Extract the input and output folder paths, and remove any trailing quotation marks and slashes
    input_folder = Path(args.input)
    output_folder = Path(args.output)
    common_files_folder = input_folder

save_folder = Path(args.save_loc)

final_flag = args.final

# Extract the camera position parameters
elev_val, azim_val = args.view
Rotate_flag = (
    not args.fixed
)  # easier to specify 'fixed' in input, but use 'rotate' as a flag higher up


# Ensure that the folders exist
# if not input_folder.is_dir():
#     print(f"Error: The input folder {input_folder} does not exist.")
#     exit(1)
# if not output_folder.is_dir():
#     print(f"Warning: The output folder {output_folder} does not exist.")
#     exit(1)
# if not common_files_folder.is_dir():
#     print(f"Warning: The common files folder {common_files_folder} does not exist.")
#     exit(1)
# if not save_folder.is_dir():
#     print(f"Warning: The save folder {save_folder} does not exist.")
#     exit(1)
for _, f in enumerate([input_folder, output_folder, common_files_folder, save_folder]):
    if not f.is_dir():
        print(f"Error: {f} does not exist.")
        exit(1)


print(f" Input folder:        {input_folder}")
print(f" Common Files folder: {common_files_folder}")
print(f" Output folder:       {output_folder}")
print(f" Save folder:         {save_folder}")

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


# Build the spatial grids
X_grid = np.arange(-(nx - 1) / 2.0, (nx - 1) / 2.0 + 1.0) * dx
Y_grid = np.arange(-(ny - 1) / 2.0, (ny - 1) / 2.0 + 1.0) * dy
X, Y = np.meshgrid(X_grid, Y_grid)


# load injection locations
inj_loc_idxs = np.loadtxt(f"{output_folder}/Other/injection_locations.txt")

if np.ndim(inj_loc_idxs) == 1:
    # If there's only one injection point, this needs to be done slightly differently
    n_inj_locs = 1
    inj_grid_vals = np.zeros((1, 2))
    inj_grid_vals[0, 0] = X_grid[int(inj_loc_idxs[0])]
    inj_grid_vals[0, 1] = Y_grid[int(inj_loc_idxs[1])]
else:
    shape_inj_locs = np.shape(inj_loc_idxs)
    n_inj_locs = shape_inj_locs[0]
    inj_grid_vals = np.zeros((n_inj_locs, 2))
    for k in range(0, n_inj_locs):
        inj_grid_vals[k, 0] = X_grid[int(inj_loc_idxs[k, 0])]
        inj_grid_vals[k, 1] = Y_grid[int(inj_loc_idxs[k, 1])]


# Reservoir Ceiling and Basement topography
H0 = load_common_files(common_files_folder, input_folder, "ceil_topo", 0, 0)
B0 = load_common_files(common_files_folder, input_folder, "base_topo", 0, 0)

# Porosity and Permeability fields
Porosity = read_in_3D_array(common_files_folder, input_folder, "porosity", nx, ny, nz)
Permeability = read_in_3D_array(
    common_files_folder, input_folder, "permeability", nx, ny, nz
)


## Load plot data into arrays #####################################################################

# Free and trapped CO2 thicknesses
h_array = np.zeros([ny, nx, n_plot])
h_res_array = np.zeros([ny, nx, n_plot])

# Free and trapped CO2 local volumes
V_mobile_array = np.zeros([ny, nx, n_plot])
V_trapped_array = np.zeros([ny, nx, n_plot])


# Load data from each of the output times and fill in the corresponding arrays
for i, t in enumerate(plot_times):

    V_mobile_array[:, :, i] = np.loadtxt(
        f"{output_folder}/Current_Volume/V" + "{0:02d}".format(i) + ".txt"
    )

    V_trapped_array[:, :, i] = np.loadtxt(
        f"{output_folder}/Current_Volume/V_res" + "{0:02d}".format(i) + ".txt"
    )


# Load the volume data recorded at each step, rather than just at the output times
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


# Calculate the mobile and trapped volumes for data saved at each plot time,
# rather than at each step
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


###################################################################################################
## Plot preparation ###############################################################################
###################################################################################################

min_H0 = np.min(H0)
max_H0 = np.max(H0)
range_H0 = max_H0 - min_H0

min_B0 = np.min(B0)
max_B0 = np.max(B0)
range_B0 = max_B0 - min_B0

D0 = B0 - H0
min_D0 = np.min(D0)
max_D0 = np.max(D0)
range_D0 = max_D0 - min_D0

print(f"{range_H0=}, {min_H0=}, {max_H0=}")
print(f"{range_B0=}, {min_B0=}, {max_B0=}")
print(f"{range_D0=}, {min_D0=}, {max_D0=}")

## Extreme values #################################################################################

# Print appropriate maximum and minimum values
max_c_mobile = np.amax(Volumes_mobile)
print("Maximum active volume is             " + str(max_c_mobile))

max_c_trapped = np.amax(Volumes_trapped)
print("Maximum trapped volume is            " + str(max_c_trapped))

# max_col = np.max(Volumes_mobile)
max_col = np.max(Volumes_injected)


# Local maxima from arrays
max_c_mobile_array = np.max(V_mobile_array)
max_c_trapped_array = np.max(V_trapped_array)

# Local maxima from iteration outputs
max_c_mobile_local = np.max(Volumes_mobile_local_max)
max_c_trapped_local = np.max(Volumes_trapped_local_max)

# Log10 min and max of permeability
Log10_perm_min = np.log10(np.amin(Permeability))
Log10_perm_max = np.log10(np.amax(Permeability))


##  Colour maps ###################################################################################

colmap_V_mobile = cm.BuPu  # cm.BuGn
colmap_V_trapped = cm.Reds  # cm.OrRd
colmap_P = cm.seismic

colmap_perm = cm.gray
colmap_poro = cm.gray

## Threshold values ###############################################################################

# Threshold thickness value for masking and bounding contour - used in determining the
# edge of the current
V_mobile_threshold = 1e-8
V_trapped_threshold = 1e-8
P_threshold = 1e-4

## Contours #######################################################################################

# Remaining percentage contours to plot
V_mobile_percentages = np.array(
    [0.1, 0.5, 0.9]
)  # Scale by current max value within loop below
cont_V_mobile_strs = ["10%", "50%", "90%"]

V_trapped_percentages = np.array([0.1, 0.5, 0.9])
cont_V_trapped_strs = ["10%", "50%", "90%"]


# Colours for the bounding contours of the free and trapped CO2 regions
clr_V = "green"
clr_V_res = "tab:brown"

## 3D layers ######################################################################################

# Layers to display in the 3D plots of reservoir geometry
z = np.linspace(0.0, 1.0, nz)
n_layers = 4
layer_indices = evenly_spaced_layers(nz, n_layers)


###################################################################################################
## Plot results ###################################################################################
###################################################################################################

start = timer()

# parallel_plot()

# if __name__ == '__main__':
#     with Pool(4) as p:
#         p.map(single_plot, range(n_plot))

# if __name__ == '__main__':
#     with ThreadPool(4) as p:
#         p.map(single_plot, range(n_plot))


fig, axs, cbaxs = build_figure()

# At the moment the layout for the first plot (final plot time) is slightly different to the rest,
# and this stands out when they're compiled as a gif. My quick hack for this is to generate this
# first plot twice, so that the second iteration of it should appear as the rest of the plots do.
single_plot(n_plot - 1, fig, axs, cbaxs)

# ## Main loop over each of the output plots
# # Start with the final plot time, so it can be viewed first
if final_flag:
    # Only produce the final plot
    n_end = n_plot - 2
else:
    # Produce every plot
    n_end = -1

# Plotting loop
for i in range(n_plot - 1, n_end, -1):
    single_plot(i, fig, axs, cbaxs)

end = timer()
print(f"\n Time take to run: {end - start} [s]")
