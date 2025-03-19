# Reservoir_preview.py
# This script plots the Topography, Permeability, and Porosity arrays, along with
# the injection locations, as currently specified in the Input folder

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import xml.etree.ElementTree as ET
from pathlib import Path
import argparse

###################################################################################################


def read_in_3D_array(filepath, nx, ny, nz):

    Array = np.zeros((nx, ny, nz))

    with open(filepath) as f:
        Lines = f.readlines()

        l_count = 0
        for k in range(0, nz):
            for j in range(0, ny):
                L = Lines[l_count].strip()
                L = L.split()
                Array[:, j, k] = [float(x) for x in L]
                l_count += 1

    return Array


def generate_colourmap(Array, colmap, minn, maxx):

    # create colormap according to h value, masked with a threshold value so that it doesn't colour the h=0 regions
    norm = colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap=colmap)
    m.set_array([])
    fcolors = m.to_rgba(Array)

    return fcolors


def evenly_spaced_layers(N, n_layers):
    # Built with help from ChatGPT 3.5

    if n_layers <= 1:
        return [0, N - 1]

    indices = np.round(np.linspace(0, N - 1, n_layers)).astype(int)

    indices[-1] = N - 1

    return indices


def get_grid_parameters(root, ns):

    # Find the <grid_parameters> block
    grid_params = root.find("cml:grid_parameters", ns)

    # Check the block exists, and extract the parameter values
    # We need to strip off the namespace tags that appears at the start
    # of each child when read in ( e.g. {http://www.xml-cml.org/schema}grid_parameters )
    if grid_params is not None:
        parameters = {child.tag.split("}")[-1]: child.text for child in grid_params}
    else:
        print("Error: <grid_parameters> block not found")

    # Check that necessary parameters are present
    required_keys = ["nx", "ny", "nz", "dx", "dy"]
    for key in required_keys:
        if key not in parameters:
            raise ValueError(f"Error: Can't read '{key}' from <grid_parameters>")

    # Convert values to appropriate types
    nx = int(parameters["nx"])
    ny = int(parameters["ny"])
    nz = int(parameters["nz"])
    dx = float(parameters["dx"])
    dy = float(parameters["dy"])

    return nx, ny, nz, dx, dy


def get_injection_locations(root, ns):

    # Find the injection location block
    # Find the <injection_locations> block
    injection_locations = root.find("cml:injection_locations", ns)

    # Extract (i, j) pairs
    locations = []
    if injection_locations is not None:
        for loc in injection_locations.findall("cml:location", ns):
            i = int(loc.find("cml:i", ns).text)
            j = int(loc.find("cml:j", ns).text)
            locations.append((i, j))
    else:
        raise ValueError("Error: <injection_locations> block not found")

    # Convert list to NumPy array
    locations_array = np.array(locations)
    n_inj_locs = locations_array.shape[0]

    # Print results
    print(f"Number of locations: {n_inj_locs}")
    print("Locations array:\n", locations_array)

    return n_inj_locs, locations_array


###################################################################################################


def main():

    ### Deal with inputs for input and output data folders ######################################

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Input data folder and XML file")

    # Add argument for input folder with default value
    parser.add_argument(
        "--input", type=str, default="./Input/", help="Input data folder path"
    )

    # Add arguments for input XML file
    parser.add_argument("--xml", type=str, help="Input XML file")

    # Parse the arguments
    args = parser.parse_args()

    # Extract the input and output folder paths, and remove any trailing quotation marks and slashes
    input_folder = Path(args.input)

    # Ensure that the folders exist
    if not input_folder.is_dir():
        print(f"Error: The input folder {input_folder} does not exist.")
        exit(1)

    print(f" Input folder : {input_folder}")

    if args.xml is None:
        print(f"Error: No XML path has been provided using --xml.")
        exit(1)

    XML_file = Path(args.xml)

    if not XML_file.is_file():
        print(f"Error: The input XML file {XML_file} does not exist.")
        exit(1)

    print(f" Input XML file: {XML_file}")

    #############################################################################################

    # Load in relevant XML data
    tree = ET.parse(XML_file)
    root = tree.getroot()

    ns = {"cml": "http://www.xml-cml.org/schema"}

    nx, ny, nz, dx, dy = get_grid_parameters(root, ns)

    # Print values
    print(f"nx: {nx}, ny: {ny}, nz: {nz}, dx: {dx}, dy: {dy}")

    n_inj_locs, inj_loc_idxs = get_injection_locations(root, ns)

    print(f"{inj_loc_idxs = }")
    print(f"{inj_loc_idxs[0] = }")

    #############################################################################################

    # # # load grid parameter values
    # # with open(f"{input_folder}/grid_parameters.txt") as f:
    # #     # Skip the initial comment lines
    # #     lines = (line for line in f if not line.startswith("--"))
    # #     parameters = np.loadtxt(lines)
    # #     # nx ny dx dy
    # #     nx, ny, nz, dx, dy = parameters
    # #     nx, ny, nz = [ int(i) for i in [nx,ny,nz] ] #Convert these to integers

    # Spatial grids
    X_grid = np.arange(-(nx - 1) / 2.0, (nx - 1) / 2.0 + 1.0) * dx
    Y_grid = np.arange(-(ny - 1) / 2.0, (ny - 1) / 2.0 + 1.0) * dy

    X, Y = np.meshgrid(X_grid, Y_grid)

    inj_grid_vals = np.zeros((n_inj_locs, 2))
    for k in range(0, n_inj_locs):
        inj_grid_vals[k, 0] = X_grid[int(inj_loc_idxs[k, 0])]
        inj_grid_vals[k, 1] = Y_grid[int(inj_loc_idxs[k, 1])]

    # # # load injection locations
    # # with open(f"{input_folder}/injection_locations.txt") as f:
    # #     # Skip the initial comment lines
    # #     lines = (line for line in f if not line.startswith("--"))
    # #     # lines now consists of the number of injection locations,
    # #     # followed by the index pairs of their locations
    # #     lines_list = list(lines)
    # #     n_inj_locs = int(lines_list[0].strip())
    # #     inj_locs_data = np.loadtxt(lines_list[1:])

    # Convert these index pairs into spatial locations for plotting
    # # inj_grid_vals = np.zeros((n_inj_locs, 2))

    # # if n_inj_locs == 1:
    # #     inj_grid_vals[0, 0] = X_grid[int(inj_locs_data[0])]
    # #     inj_grid_vals[0, 1] = Y_grid[int(inj_locs_data[1])]
    # # else:
    # #     for k in range(0, n_inj_locs):
    # #         inj_grid_vals[k, 0] = X_grid[int(inj_locs_data[k, 0])]
    # #         inj_grid_vals[k, 1] = Y_grid[int(inj_locs_data[k, 1])]

    # Read in Topography, Porosity, and Permeability
    H0 = np.loadtxt(f"{input_folder}/ceil_topo.txt")
    B0 = np.loadtxt(f"{input_folder}/base_topo.txt")
    D0 = B0 - H0

    Porosity = read_in_3D_array(f"{input_folder}/porosity.txt", nx, ny, nz)
    Permeability = read_in_3D_array(f"{input_folder}/permeability.txt", nx, ny, nz)

    # Print relevant stats about static reservoir properties
    print("")
    print("------------------------------")
    print(
        f"""Reservoir Vertical Extent [m]: \n
    - Min  = {np.amin(D0):f}, \n
    - Max  = {np.amax(D0):f}, \n
    - Mean = {np.mean(D0):f}
    """
    )
    print("------------------------------")
    print(
        f"""Porosity [-]: \n
    - Min  = {np.amin(Porosity):f}, \n
    - Max  = {np.amax(Porosity):f}, \n
    - Mean = {np.mean(Porosity):f}
    """
    )
    print("------------------------------")
    print(
        f"""Permeability [mD]: \n
    - Min  = {np.amin(Permeability):f}, \n
    - Max  = {np.amax(Permeability):f}, \n
    - Mean = {np.mean(Permeability):f}
    """
    )
    print("------------------------------")
    print("")

    # Specify the colour maps for the topography, permeability and porosity, and
    # the numerical ranges to use for them
    colmap_ceil = cm.seismic
    norm_ceil = colors.Normalize(np.min(H0), np.max(H0))

    colmap_base = cm.seismic
    norm_base = colors.Normalize(np.min(B0), np.max(B0))

    colmap_D0 = cm.Blues
    norm_D0 = colors.Normalize(np.min(D0), np.max(D0))

    colmap_perm = cm.gray
    norm_perm = colors.Normalize(0.0, np.max(Permeability))

    colmap_poro = cm.gray
    norm_poro = colors.Normalize(0.0, np.max(Porosity))

    # Format string to be used for tick labels in colourbars
    cbar_fmt = lambda x, pos: "{:.1f}".format(x)

    # Layers for plotting the 3D structure of porosity and permeability
    z = np.linspace(0.0, 1.0, nz)
    n_layers = 1  # 4
    layer_indices = evenly_spaced_layers(nz, n_layers)

    fig = plt.figure(figsize=(12, 8))
    # gs = GridSpec(2, 2, figure=fig)

    gs = GridSpec(2, 1, figure=fig)
    gs_top = gs[0].subgridspec(1,3)
    gs_bot = gs[1].subgridspec(1,2)

    # -- Ceiling Topography -----------------------------------------------------------------
    ax_c = fig.add_subplot(gs_top[0]) #gs[0, 0])

    # Contours of ceiling topography
    cont_H0 = ax_c.contour(X_grid, Y_grid, H0, colors="gray", alpha=0.25, levels=10)
    ax_c.imshow(
        np.flipud(H0),
        cmap=colmap_ceil,
        norm=norm_ceil,
        extent=[X.min(), X.max(), Y.min(), Y.max()],
    )

    # Plot the injection locations
    for k in range(0, n_inj_locs):
        ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "ko")

    # Labels for the contour maps
    # ax_c.clabel(cont_H0, inline=True, fmt="H0 = %.2g", fontsize=10)

    ax_c.set_title("Ceiling Topography")
    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

    # Add a colourbar
    m = cm.ScalarMappable(cmap=colmap_base, norm=norm_base)
    m.set_array(H0)
    axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad= -2)
    fig.colorbar(m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt))

    # -- Basement Topography ----------------------------------------------------------------
    ax_c = fig.add_subplot(gs_top[1]) #gs[0, 1])

    # Contours of basement topography
    cont_B0 = ax_c.contour(X_grid, Y_grid, B0, colors="gray", alpha=0.25, levels=10)
    ax_c.imshow(
        np.flipud(B0),
        cmap=colmap_base,
        norm=norm_base,
        extent=[X.min(), X.max(), Y.min(), Y.max()],
    )

    # Plot the injection locations
    for k in range(0, n_inj_locs):
        ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "ko")

    # Labels for the contour maps
    ax_c.clabel(cont_B0, inline=True, fmt="B0 = %.2g", fontsize=10)

    ax_c.set_title("Basement Topography")
    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

    # Add a colourbar
    m = cm.ScalarMappable(cmap=colmap_base, norm=norm_base)
    m.set_array(B0)
    axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad= -2)
    fig.colorbar(m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt))


    # -- Reservoir Vertical Extent ----------------------------------------------------------------
    ax_c = fig.add_subplot(gs_top[2]) #gs[0, 1])

    # Contours of the vertical extent of the storage interval
    cont_D0 = ax_c.contour(X_grid, Y_grid, D0, colors="gray", alpha=0.25, levels=10)
    ax_c.imshow(
        np.flipud(D0),
        cmap=colmap_D0,
        norm=norm_D0,
        extent=[X.min(), X.max(), Y.min(), Y.max()],
    )

    # Plot the injection locations
    for k in range(0, n_inj_locs):
        ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "ko")

    # Labels for the contour maps
    ax_c.clabel(cont_D0, inline=True, fmt="D0 = %.2g", fontsize=10)

    ax_c.set_title("Reservoir Vertical Extent")
    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

    # Add a colourbar
    m = cm.ScalarMappable(cmap=colmap_D0, norm=norm_D0)
    m.set_array(D0)
    axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad= -2)
    fig.colorbar(m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt))


    # -- Porosity --------------------------------------------------------------------------

    ax_c = fig.add_subplot(gs_bot[0], projection="3d") # gs[1, 0]

    for k in layer_indices:
        fcolors_poro = generate_colourmap(Porosity[:, :, k], colmap_poro, np.min(Porosity), np.max(Porosity))

        ax_c.plot_surface(
            X,
            Y,
            H0 + z[k] * (B0 - H0),
            edgecolor="none",
            facecolors=np.transpose(fcolors_poro, (1, 0, 2)),
            vmin=0.0,
            vmax=1.0,
            shade=True,
            alpha=0.25,
        )

        # # Plot injection locations
        # for l in range(0, n_inj_locs):
        #     z_val = H0[int(inj_loc_idxs[l, 1]), inj_loc_idxs[l, 0]] + z[k] * (
        #         B0[int(inj_loc_idxs[l, 1]), inj_loc_idxs[l, 0]]
        #         - H0[int(inj_loc_idxs[l, 1]), inj_loc_idxs[l, 0]]
        #     )
        #     print(f"{z_val = }")
        #     ax_c.plot(
        #         [inj_grid_vals[l, 0]],
        #         [inj_grid_vals[l, 1]],
        #         [0.95 * z_val],
        #         "mx",
        #     )

    # Plot injection locations
    for l in range(0, n_inj_locs):
        idxs = inj_loc_idxs[l, :]
        grid_vals = inj_grid_vals[l, :]
        ax_c.plot(
            [grid_vals[0], grid_vals[0]],
            [grid_vals[1], grid_vals[1]],
            [H0[idxs[1], idxs[0]], B0[idxs[1], idxs[0]]],
            "m",
            lw=3,
        )

    ax_c.set_title("Porosity")
    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)
    ax_c.set_zlabel(r"$z$", fontsize=20, rotation=0)

    ax_c.invert_zaxis()

    # ax_c.view_init(
    #     elev=elev_val,
    #     azim=azim_val + Rotate_flag * 360 * (i / (len(plot_times) - 1)),
    #     roll=None,
    #     vertical_axis="z",
    #     share=False,
    # )

    m = cm.ScalarMappable(cmap=colmap_poro, norm=norm_poro)
    m.set_array(H0)

    axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad=-6)

    fig.colorbar(m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt))

    # -- Permeability ----------------------------------------------------------------

    ax_c = fig.add_subplot(gs_bot[1], projection="3d") #gs[1, 1]

    for k in layer_indices:
        fcolors_perm = generate_colourmap(
            Permeability[:, :, k],
            colmap_perm,
            np.min(Permeability),
            np.max(Permeability),
        )

        ax_c.plot_surface(
            X,
            Y,
            H0 + z[k] * (B0 - H0),
            edgecolor="none",
            facecolors=np.transpose(fcolors_perm, (1, 0, 2)),
            vmin=0.0,
            vmax=1.0,
            shade=False,
            alpha=0.25,
        )

    # Plot injection locations
    for l in range(0, n_inj_locs):
        idxs = inj_loc_idxs[l, :]
        grid_vals = inj_grid_vals[l, :]
        ax_c.plot(
            [grid_vals[0], grid_vals[0]],
            [grid_vals[1], grid_vals[1]],
            [H0[idxs[1], idxs[0]], B0[idxs[1], idxs[0]]],
            "m",
            lw=3,
        )

    ax_c.set_title("Permeability")
    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)
    ax_c.set_zlabel(r"$z$", fontsize=20, rotation=0)

    ax_c.invert_zaxis()

    # ax_c.view_init(
    #     elev=elev_val,
    #     azim=azim_val + Rotate_flag * 360 * (i / (len(plot_times) - 1)),
    #     roll=None,
    #     vertical_axis="z",
    #     share=False,
    # )

    m = cm.ScalarMappable(cmap=colmap_perm, norm=norm_perm)
    m.set_array(H0)

    axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad=-6)

    fig.colorbar(m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt))

    # # plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    gs.tight_layout(fig)

    # Save this frame
    plt.savefig("./plots/Reservoir_preview/Reservoir_preview.png")

    plt.show()

    # plt.close()


if __name__ == "__main__":
    main()
