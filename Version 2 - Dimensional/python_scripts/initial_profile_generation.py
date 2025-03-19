# initial_profile_generation.py
# This script generates initial profiles for mobile CO2 thickness (h)
# and residually trapped CO2 thickness (h_res) for use in CO2GraVISim


import numpy as np
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# -------------------------------------------------------------------------------------------------
# -- Functions ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


# build_a_semi-elipsoid
def build_semielipsoid(x_centre, y_centre, r_x, r_y, height, x, y, nx, ny):

    # F = np.zeros((nx,ny))

    z_squared = (height**2) * (
        1.0 - ((x - x_centre) / r_x) ** 2 - ((y - y_centre) / r_y) ** 2
    )

    F = np.where(z_squared >= 0, np.sqrt(z_squared * (z_squared >= 0)), 0.0)
    # F = z_squared

    return F


# -------------------------------------------------------------------------------------------------
def write_array_to_file(Array, target_folder, filename):
    # Built with help from ChatGPT 3.5

    s = np.shape(Array)
    # nx = s[0]
    ny = s[1]

    # file_path = f"{target_folder}/{filename}.txt"
    file_path = Path(target_folder) / f"{filename}.txt"

    if Array.ndim > 2:
        nz = s[2]
        with open(file_path, "w") as f:
            for k in range(0, nz):
                for j in range(0, ny):
                    f.write(
                        " ".join(map(lambda x: f"{np.float32(x):.6e}", Array[:, j, k]))
                        + "\n"
                    )

    else:
        # nz = 1
        with open(file_path, "w") as f:
            for j in range(0, ny):
                f.write(
                    " ".join(map(lambda x: f"{np.float32(x):.6e}", Array[:, j])) + "\n"
                )          


# -------------------------------------------------------------------------------------------------
def generate_initial_profiles(nx,ny,dx,dy,h_init_params,target_folder,Plot_flag):

    #############################################################################################

    # Spatial grids

    print(f"{nx = }, {ny = }, {dx = }, {dy = }")

    X_grid = np.arange(-(nx - 1) / 2.0, (nx - 1) / 2.0 + 1.0) * dx
    Y_grid = np.arange(-(ny - 1) / 2.0, (ny - 1) / 2.0 + 1.0) * dy

    Y, X = np.meshgrid(Y_grid, X_grid)

    # -- Build initial profiles -------------------------------------------------------------------

    print("Building initial profiles for h, h_init, and P.")

    x_c, y_c, r_x, r_y, h_max = h_init_params

    print(
    f"Parameter values are \n\
    {x_c   = }, \n\
    {y_c   = }, \n\
    {r_x   = }, \n\
    {r_y   = }, \n\
    {h_max = }."
    )

    # Mobile current thicknes
    h = build_semielipsoid(x_c, y_c, r_x, r_y, h_max, X, Y, nx, ny)
    # Make sure there are no -0. entries, just in case that's an issue down the line
    h[h <= 0.0] = 0.0

    # Residually trapped current thickness
    h_res = np.zeros_like(h)

    # Non-hydrostatic ambient pressure
    P = np.zeros_like(h)

    # -- Save data --------------------------------------------------------------------------------

    print(f"Saving initial profiles to {target_folder}.")

    write_array_to_file(h, target_folder, "h_init")
    write_array_to_file(h_res, target_folder, "h_res_init")
    write_array_to_file(P, target_folder, "P_init")

    # -- Plot -------------------------------------------------------------------------------------

    if Plot_flag:

        print("Plotting initial profiles for h, h_init, and P.")

        # # Convert injection location index pairs into spatial locations for plotting
        # inj_grid_vals = np.zeros((n_inj_locs, 2))

        # if n_inj_locs == 1:
        #     inj_grid_vals[0, 0] = X_grid[int(inj_loc_idxs[0])]
        #     inj_grid_vals[0, 1] = Y_grid[int(inj_loc_idxs[1])]
        # else:
        #     for k in range(0, n_inj_locs):
        #         inj_grid_vals[k, 0] = X_grid[int(inj_loc_idxs[k, 0])]
        #         inj_grid_vals[k, 1] = Y_grid[int(inj_loc_idxs[k, 1])]

        fig = plt.figure(figsize=(12, 6), constrained_layout=True)
        gs = GridSpec(1, 3, figure=fig)

        ## h plot
        ax_c = fig.add_subplot(gs[0, 0], projection="3d")

        ax_c.plot_surface(
            X,
            Y,
            h,
            edgecolor="none",
            shade=False,
        )

        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20)
        ax_c.set_title(r"$h$", fontsize=20)

        ## h_res plot
        ax_c = fig.add_subplot(gs[0, 1], projection="3d")

        ax_c.plot_surface(
            X,
            Y,
            h_res,
            edgecolor="none",
            shade=False,
        )

        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20)
        ax_c.set_title(r"$h_{res}$", fontsize=20)

        ## P plot
        ax_c = fig.add_subplot(gs[0, 2], projection="3d")

        ax_c.plot_surface(
            X,
            Y,
            P,
            edgecolor="none",
            shade=False,
        )

        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20)
        ax_c.set_title(r"$P$", fontsize=20)

        plt.show()

# -------------------------------------------------------------------------------------------------
# -- Main Code ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


def main():

    ### Deal with inputs for input and output data folders ######################################

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Input folder and Elliptical Parameters for initial CO2 plume profile.")

    # Add arguments for input folder with default value
    parser.add_argument(
        "--input", type=str, default="./Input/", help="Input data folder path"
    )

    # Grid parameters
    parser.add_argument(
        "--grid",
        type=float,
        nargs="+",
        default=[200, 350, 20, 2.25, 2.25],
        help="Grid parameters nx, ny, nz, dx, and dy. \
                    The default values are [200, 350, 20, 2.25, 2.25].",
    )

    # Add parameters for the elliptical CO2 plume
    parser.add_argument(
        "--params",
        type=float,
        nargs="+",
        default=[-20.0, -40.0, 30.0, 25.0, 0.0],
        help="Parameters for the initial elliptical CO2 plume.\
            Default values are (x_c, y_c, r_x, r_y, h_max) = [-20.0, -40.0, 30.0, 25.0, 0.0].",
    )

    # Flag for whether to plot
    parser.add_argument("--plot", action="store_true")

    # Parse the arguments
    args = parser.parse_args()

    # Extract the input and output folder paths, and remove any trailing quotation marks and slashes
    input_folder = Path(args.input)

    # Ensure that the folders exist
    if not input_folder.is_dir():
        print(f"Error: The input folder {input_folder} does not exist.")
        exit(1)

    print(f" Input folder : {input_folder}")

    # Grid parameters 
    grid_params = args.grid

    # Parameter values for the elliptical CO2 profile
    h_init_params = args.params

    print(f"{h_init_params = }")

    #############################################################################################

    # # load grid parameter values
    # with open(f"{input_folder}/grid_parameters.txt") as f:
    #     # Skip the initial comment lines
    #     lines = (line for line in f if not line.startswith("--"))
    #     parameters = np.loadtxt(lines)
    #     # nx ny dx dy
    #     nx, ny, nz, dx, dy = parameters
    #     nx, ny, nz = [int(i) for i in [nx, ny, nz]]  # Convert these to integers

    nx, ny, nz, dx, dy = grid_params
    nx, ny, nz = [int(i) for i in [nx, ny, nz]]  # Convert these to integers

    # # # load injection locations
    # # with open(f"{input_folder}/injection_locations.txt") as f:
    # #     # Skip the initial comment lines
    # #     lines = (line for line in f if not line.startswith("--"))
    # #     # lines now consists of the number of injection locations,
    # #     # followed by the index pairs of their locations
    # #     lines_list = list(lines)
    # #     n_inj_locs = int(lines_list[0].strip())
    # #     inj_loc_idxs = np.loadtxt(lines_list[1:])

    Plot_flag = args.plot


    generate_initial_profiles(nx,ny,dx,dy,h_init_params,input_folder,Plot_flag)

    


if __name__ == "__main__":
    main()
