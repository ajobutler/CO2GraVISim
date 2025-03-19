# input_generation.py
# This script generates inputs to be used by CO2GraVISim


import numpy as np
import warnings
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
from matplotlib import colors

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from Write_XML_input import generate_XML

# -------------------------------------------------------------------------------------------------
# -- Functions ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------
def build_block(nx, ny, phi_a, phi_b, X_grid, Y_grid):

    poro_h = np.zeros((nx, ny))

    # Build a low-porosity rectangle a prescribed nondimensional distance away from the origin
    for j in range(0, ny):
        for i in range(0, nx):

            if X_grid[i] >= -3 and X_grid[i] <= 3 and Y_grid[j] >= 2 and Y_grid[j] <= 5:
                poro_h[i, j] = phi_a
            else:
                poro_h[i, j] = phi_b

    return poro_h


# -------------------------------------------------------------------------------------------------
def build_channel(nx, ny, phi_a, phi_b, X_grid):

    poro_h = np.zeros((nx, ny))

    # Build a low-porosity rectangle a prescribed nondimensional distance away from the origin
    for j in range(0, ny):
        for i in range(0, nx):

            if X_grid[i] >= 2 and X_grid[i] <= 5:
                poro_h[i, j] = phi_a
            else:
                poro_h[i, j] = phi_b

    return poro_h


# -------------------------------------------------------------------------------------------------
def build_checkboard(nx, ny, phi_a, phi_b, tile_nx, tile_ny):

    poro_h = np.zeros((nx, ny))

    phi_sum = phi_a + phi_b
    phi_diff = phi_a - phi_b

    num_tiles_x = int(np.ceil((nx - 1) / tile_nx))
    num_tiles_y = int(np.ceil((ny - 1) / tile_ny))

    print(
        "Number of tiles in x and y: "
        + str(num_tiles_x)
        + ", "
        + str(num_tiles_y)
        + "\n"
    )

    # Start by constructing a checkerboard pattern of +/- 1 values
    # Build up the first column
    for i in range(0, num_tiles_x - 1):
        # All but the last row, first column
        poro_h[i * tile_nx : (i + 1) * tile_nx, 0:tile_ny] = (-1.0) ** (i)
    # Last row, first column
    poro_h[(num_tiles_x - 1) * tile_nx : nx, 0:tile_ny] = (-1.0) ** (num_tiles_x - 1)

    # Build the remaining columns by alternating the sign from the first column.
    # First do this for the 2nd to penultimate columns
    for j in range(1, num_tiles_y - 1):
        poro_h[:, j * tile_ny : (j + 1) * tile_ny] = (-1.0) ** (j) * poro_h[
            :, 0:tile_ny
        ]

    # Now do this for the final column
    # We may not have a full column here, so fill in the rest of this part of the array
    poro_h[:, (num_tiles_y - 1) * tile_ny : ny] = (
        (-1.0) ** (num_tiles_y - 1)
    ) * poro_h[:, 0 : (ny - (num_tiles_y - 1) * tile_ny)]

    # Switch from +/- 1 values to phia/phib values
    poro_h = 0.5 * (phi_sum + phi_diff * poro_h)

    return poro_h


# -------------------------------------------------------------------------------------------------
def build_v_linear(z, phi_ratio):

    nz = len(z)

    # z ranges fom 0 to 1 across the reservoir
    poro_v = 1.0 - (1.0 - phi_ratio) * z

    return poro_v


# -------------------------------------------------------------------------------------------------
def build_v_compaction(z, phi_ratio):

    nz = len(z)

    if phi_ratio >= 1.0:
        poro_v = np.ones((nz,))
    else:
        L = -z[-1] / np.log(phi_ratio)
        poro_v = np.exp(-z / L)

    return poro_v


# -------------------------------------------------------------------------------------------------
def build_v_layers(z, phi_ratio, width_1, width_2):

    nz = len(z)

    poro_v = np.zeros((nz,))

    n_1 = int(np.floor(nz * width_1))
    n_2 = int(np.floor(nz * width_2))

    # 'b' layers
    poro_v[:] = phi_ratio
    # 'a' layers
    for k in range(0, nz, (n_1 + n_2)):
        poro_v[k : k + n_1] = 1.0

    print(n_1, n_2, nz)
    print(np.max(poro_v), np.min(poro_v))

    return poro_v


# -------------------------------------------------------------------------------------------------
def build_full_poro_array(nx, ny, nz, poro_h, poro_v):

    poro = np.zeros((nx, ny, nz))

    for k in range(0, nz):
        poro[:, :, k] = poro_h[:, :] * poro_v[k]

    return poro


# -------------------------------------------------------------------------------------------------
def build_topography(nx, ny, X_grid, Y_grid, slope_params, bump_params, ceil_mean, base_mean):

    ceil_topo = topography_fn(
        nx, ny, X_grid, Y_grid, slope_params[0:2], bump_params[0:5], ceil_mean
    )
    base_topo = topography_fn(nx, ny, X_grid, Y_grid, slope_params[2:], bump_params[5:], base_mean)

    # Make sure that the base and ceiling and sufficiently far apart
    # base_topo = base_topo - np.min(np.abs(base_topo - ceil_topo)) + 1.0

    # base_topo = base_topo + 1.0  # Add a consistent offset - bumps could intersect now!
    # base_topo = base_topo + 10.0  # Add a consistent offset - bumps could intersect now!

    base_topo = np.maximum(
        ceil_topo + 1e-3, base_topo
    )  # Make sure there aren't any intersections

    return ceil_topo, base_topo


# -------------------------------------------------------------------------------------------------
def topography_fn(nx, ny, X_grid, Y_grid, C_slope, C_bump,mean_depth):

    angle_x, angle_y = [C_slope[i] for i in (0, 1)]
    amp, L_x, L_y, W_x, W_y = [C_bump[i] for i in (0, 1, 2, 3, 4)]

    topo = np.zeros((nx, ny))

    for i in range(0, nx):
        for j in range(0, ny):
            topo[i, j] = (
                np.tan(angle_x * np.pi / 180.0) * X_grid[i]  # dx * (i - (nx - 1) / 2)
                + np.tan(angle_y * np.pi / 180.0) * Y_grid[j]  # dy * (j - (ny - 1) / 2)
                + amp
                * (
                    np.cos(2 * np.pi * X_grid[i] / L_x)  # / X_grid[-1]
                    * np.exp(-0.5 * ((X_grid[i] / W_x) ** 2))
                )
                * (
                    np.cos(2 * np.pi * Y_grid[j] / L_y)  # / Y_grid[-1]
                    * np.exp(-0.5 * ((Y_grid[j] / W_y) ** 2))
                )
            )

    topo = topo + mean_depth

    
    return topo


# -------------------------------------------------------------------------------------------------
def generate_colourmap(Array, colmap, minn, maxx):

    # create colormap according to h value, masked with a threshold value so that it doesn't colour
    # the h=0 regions
    norm = colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap=colmap)
    m.set_array([])
    fcolors = m.to_rgba(Array)

    return fcolors


# -------------------------------------------------------------------------------------------------
def evenly_spaced_layers(N, n_layers):
    # Built with help from ChatGPT 3.5

    if n_layers <= 1:
        return [0, N - 1]

    indices = np.round(np.linspace(0, N - 1, n_layers)).astype(int)

    indices[-1] = N - 1

    return indices


# -------------------------------------------------------------------------------------------------
def write_array_to_file(Array, target_folder, filename):
    # Built with help from ChatGPT 3.5

    s = np.shape(Array)
    # nx = s[0]
    ny = s[1]

    file_path = f"{target_folder}/{filename}.txt"

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
def write_default_XML(nx,ny,nz,dx,dy,input_folder):

    ### flow parameters
    flow_parameters_dict = {
        "rho_c": 650.0,
        "rho_a_unsat": 1030.0,
        "rho_a_sat": 1042.0,
        "mu_c": 7e-5,
        "mu_a": 9e-4,
        "s_c_r": 0.3,
        "s_a_i": 0.18,
        "krn_mobile": 0.5,
        "krw_residual": 1.0,
        "C_sat": 0.04,
        "g": 9.81,
        "D_mol": 2e-9,
        "perm_ratio": 0.1,
        "omega_conv": 1.0,
    }

    ### Grid parameters
    grid_parameters_dict = {"nx": nx, "ny": ny, "nz": nz, "dx": dx, "dy": dy}

    ### Boundary Conditions
    BC_dict = {
        "current_thickness": {"north": 1, "east": 1, "south": 1, "west": 1},
        "ambient_pressure": {"north": 1, "east": 1, "south": 1, "west": 1},
    }

    ### Plot times
    Plot_times_dict = {
        "times": [
            0.00000,
            365.00000,
            730.00000,
            1096.00000,
            1461.00000,
            1826.00000,
            2191.00000,
            2557.00000,
            2922.00000,
            3287.00000,
            3652.00000,
            4018.00000,
            4383.00000
        ]
    }


    ### Injection information
    Injection_dict = {
        # 'locations' : [[100, 175]], #This has to be an array rather than a list!
        "locations": [
            [int(nx/2), int(ny/2)],
        ],  # This has to be an array rather than a list!
        "Flux_times": [0.0, 365.0, 730.0],
        "Flux_vals": [
            [0.002, 0.001, 0.00],
        ],  # This has to be an array rather than a list!
    }

    ## Generate the XML document
    # Name it after the target Input folder, unless this is the default Input folder
    if (input_folder.name == 'Input'):
        XML_filepath = input_folder / "default_XML.xml"
    else:
        XML_filepath = input_folder / f"{input_folder.name}.xml"

    generate_XML(
        grid_parameters_dict,
        BC_dict,
        flow_parameters_dict,
        Injection_dict,
        Plot_times_dict,
        XML_filepath,
    )

# -------------------------------------------------------------------------------------------------
# -- Main Code ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


def main():

    ### Deal with inputs for input and output data folders ######################################

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Input and Output data folders")

    # Add arguments for input folder with default value
    parser.add_argument(
        "--input", type=str, default="./Input/", help="Input data folder path"
    )

    # Flag for whether to plot
    parser.add_argument("--plot", action="store_true", help="Display the reservoir properties once generated")

    # Horizontal porosity structure
    parser.add_argument(
        "--ptype_h",
        type=str,
        default="Uniform",
        help="Horizontal spatial structure for the porosity: can be Uniform, Block, Channel, \
            or Checkerboard. The default is Uniform.",
    )

    # Parameters for horizontal porosity structure
    parser.add_argument(
        "--pparams_h",
        type=float,
        nargs="+",
        default=[0.75, 0.3, 5, 5],
        help="Parameters for the horizontal porosity structure. \
                These are phi_a, phi_b, tile_nx, tile_ny (only Checkerboard requires the last two). \
                Default values are [0.75, 0.3, 5, 5].",
    )

    # Vertical porosity structure
    parser.add_argument(
        "--ptype_v",
        type=str,
        default="Uniform",
        help="Vertical spatial structure for the porosity: can be Uniform, Linear, Compaction, \
            or Layers. The default is Uniform.",
    )

    # Parameters for vertical porosity structure
    parser.add_argument(
        "--pparams_v",
        type=float,
        nargs="+",
        default=[0.6, 5.0, 1.2],
        help="Parameters for the vertical porosity structure. \
                These are phi_ratio, length_scale_1, length_scale_2 \
                    (only Compaction and Layers require length_scale_1, and only Layers requires length_scale_2). \
                Default values are [0.6, 5., 1.2].",
    )

    # Slope parameters for topography
    parser.add_argument(
        "--slope",
        type=float,
        nargs="+",
        default=[0.1, -0.2, 0.1, -0.2],
        help="Slope parameters for the topography. \
                These are ceil_angle_x, ceil_angle_y, base_angle_x, base_angle_y. \
                    The default values are [0.1, -0.2, 0.1, -0.2].",
    )

    # Bump parameters for topography
    parser.add_argument(
        "--bump",
        type=float,
        nargs="+",
        default=[0.0, 9.0, 15.0, 30.0, 40.0, 0.0, 5.0, 3.0, 10.0, 9.0],
        help="Bump parameters for the topography. \
                These are ceil_bump_amp, ceil_bump_lambda_x, ceil_bump_lambda_y, ceil_bump_decay_x, ceil_bump_decay_y, \
                    base_bump_amp, base_bump_lambda_x, base_bump_lambda_y, base_bump_decay_x, base_bump_decay_y. \
                    The default values are [0.0, 9.0, 15.0, 30.0, 40.0, 0.0, 5.0, 3.0, 10.0, 9.0].",
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

    # Generate default xml file
    parser.add_argument("--default_XML", action="store_true", help="Generate a default XML input file in the target Input folder")


    # Parse the arguments
    args = parser.parse_args()

    # Extract the input and output folder paths, and remove any trailing quotation marks and slashes
    input_folder = Path(args.input)

    # Ensure that the folders exist
    if not input_folder.is_dir():
        print(f"Error: The input folder {input_folder} does not exist.")
        exit(1)

    print(f" Input folder : {input_folder}")

    Plot_flag = args.plot

    Poro_h_type = args.ptype_h
    Poro_h_params = args.pparams_h

    Poro_v_type = args.ptype_v
    Poro_v_params = args.pparams_v

    slope_params = args.slope

    bump_params = args.bump

    grid_params = args.grid

    XML_flag = args.default_XML

    #############################################################################################

    # load grid parameter values
    # # with open(f"{input_folder}/grid_parameters.txt") as f:
    # #     # Skip the initial comment lines
    # #     lines = (line for line in f if not line.startswith("--"))
    # #     parameters = np.loadtxt(lines)
    # #     # nx ny dx dy
    # #     nx, ny, nz, dx, dy = parameters
    # #     nx, ny, nz = [int(i) for i in [nx, ny, nz]]  # Convert these to integers

    nx, ny, nz, dx, dy = grid_params
    nx, ny, nz = [int(i) for i in [nx, ny, nz]]  # Convert these to integers

    # # load injection locations
    # with open(f"{input_folder}/injection_locations.txt") as f:
    #     # Skip the initial comment lines
    #     lines = (line for line in f if not line.startswith("--"))
    #     # lines now consists of the number of injection locations,
    #     # followed by the index pairs of their locations
    #     lines_list = list(lines)
    #     n_inj_locs = int(lines_list[0].strip())
    #     inj_loc_idxs = np.loadtxt(lines_list[1:])

    n_inj_locs = 1
    inj_loc_idxs = [nx / 2, ny / 2]

    input_gen(
        input_folder,
        nx,
        ny,
        nz,
        dx,
        dy,
        n_inj_locs,
        inj_loc_idxs,
        Poro_h_type,
        Poro_h_params,
        Poro_v_type,
        Poro_v_params,
        slope_params,
        bump_params,
        Plot_flag,
        XML_flag,
    )


def input_gen(
    input_folder,
    nx,
    ny,
    nz,
    dx,
    dy,
    n_inj_locs,
    inj_loc_idxs,
    Poro_h_type,
    Poro_h_params,
    Poro_v_type,
    Poro_v_params,
    slope_params,
    bump_params,
    Plot_flag,
    XML_flag,
):

    # Make these lowercase to avoid issues when comparing for the cases below
    Poro_h_type = Poro_h_type.lower()

    Poro_v_type = Poro_v_type.lower()

    # -- Dimensional Scales -----------------------------------------------------------------------

    H_scale = 1050.0 # [m]
    D_scale = 30.2   # [m]
    k_scale = 121.4  # [mD]


    # -- Reservoir Topography ---------------------------------------------------------------------

    print("nx, ny, nz, dx, dy: " + str([nx, ny, nz, dx, dy]))

    # Spatial grids
    X_grid = np.arange(-(nx - 1) / 2.0, (nx - 1) / 2.0 + 1.0) * dx
    Y_grid = np.arange(-(ny - 1) / 2.0, (ny - 1) / 2.0 + 1.0) * dy

    Y, X = np.meshgrid(Y_grid, X_grid)

    print('\n\n')
    print('Spatial grid:')
    print(f' - {nx}*{ny} with {nz} layers')
    print(f' - Physical extent: x_max = {X_grid[-1]}[m], y_max = {Y_grid[-1]}[m]')
    print('\n\n')

    print("Slope Parameters: " + str(slope_params))
    print("Bump Parameters:  " + str(bump_params))

    #Mean depths for ceiling and basement
    ceil_mean = H_scale
    base_mean = H_scale + D_scale

    ceil_topo, base_topo = build_topography(
        nx, ny, X_grid, Y_grid, slope_params, bump_params, ceil_mean, base_mean
    )

    print(f"{np.min(ceil_topo) = }, {np.max(ceil_topo) = }, {np.mean(ceil_topo) = }")

    # -- Porosity ---------------------------------------------------------------------------------

    poro_h = np.zeros((nx, ny))

    print("\n - Horizontal Structure:")

    match Poro_h_type:

        case "uniform":
            print("You chose Uniform \n")

            if len(Poro_h_params) < 1:
                print("The values of phi_a must be specified for this case. \n")
                print("parameters provided: " + str(Poro_h_params) + "\n")
                exit

            phi_a = Poro_h_params[0]

            poro_h[:, :] = phi_a

        case "block":
            print("You chose Low-Porosity Rectangle \n")

            # parameters: phi_a, phi_b

            if len(Poro_h_params) < 2:
                print(
                    "The values of phi_a and phi_b must be specified for this case. \n"
                )
                print("parameters provided: " + str(Poro_h_params) + "\n")
                exit

            phi_a, phi_b = Poro_h_params[0:2]
            print(f"{phi_a = }, {phi_b = }")

            poro_h = build_block(nx, ny, phi_a, phi_b, X_grid, Y_grid)

        case "channel":
            print("You chose High-Porosity Channel \n")

            # parameters: phi_a, phi_b

            if len(Poro_h_params) < 2:
                print(
                    "The values of phi_a and phi_b must be specified for this case. \n"
                )
                print("parameters provided: " + str(Poro_h_params) + "\n")
                exit

            phi_a, phi_b = Poro_h_params[0:2]
            print(f"{phi_a = }, {phi_b = }")

            poro_h = build_channel(nx, ny, phi_a, phi_b, X_grid)

        case "checkerboard":
            print("You chose Checkerboard \n")

            # parameters: phi_a, phi_b, tile_nx, tile_ny

            if len(Poro_h_params) < 4:
                print(
                    "The values of phi_a, phi_b, tile_nx, and tile_ny must be specified for this case. \n"
                )
                print("parameters provided: " + str(Poro_h_params) + "\n")
                exit

            phi_a, phi_b = Poro_h_params[0:2]
            tile_nx = int(np.round(Poro_h_params[2]))  # Make sure this is an integer
            tile_ny = int(np.round(Poro_h_params[3]))

            print(
                "[phi_ratio, tile_nx, tile_ny] = "
                + str([phi_a, phi_b, tile_nx, tile_ny])
                + "\n"
            )

            poro_h = build_checkboard(nx, ny, phi_a, phi_b, tile_nx, tile_ny)

        case _:
            print(
                "You chose "
                + str(Poro_h_type)
                + ". Please choose from Uniform, Block, Channel, or Checkerboard."
            )

    # -- Vertical structure -----------------------------------------------------------------------

    poro_v = np.zeros((nz,))

    z = np.linspace(0.0, 1.0, nz)

    print("\n - Vertical Structure:")

    match Poro_v_type:

        case "uniform":
            print("You chose Uniform \n")

            poro_v[:] = 1.0

        case "linear":
            print("You chose Linear \n")

            # parameters: phi_ratio

            if len(Poro_v_params) < 1:
                print("The value of phi_ratio must be specified for this case. \n")
                print("parameters provided: " + str(Poro_v_params) + "\n")
                exit

            phi_ratio = Poro_v_params[0]

            print("phi_ratio = " + str(phi_ratio))

            poro_v = build_v_linear(z, phi_ratio)

        case "compaction":
            print("You chose Compaction \n")

            # parameters: phi_ratio

            if len(Poro_v_params) < 1:
                print("The value of phi_ratio must be specified for this case. \n")
                print("parameters provided: " + str(Poro_v_params) + "\n")
                exit

            phi_ratio = Poro_v_params[0]

            print("phi_ratio = " + str(phi_ratio))

            poro_v = build_v_compaction(z, phi_ratio)

        case "layers":
            print("You chose Layers \n")

            # parameters: phi_ratio, width_a, width_b

            if len(Poro_v_params) < 3:
                print(
                    "The values of phi_ratio, width_a, and width_b must be specified for this case. \n"
                )
                print("parameters provided: " + str(Poro_v_params) + "\n")
                exit

            phi_ratio, width_a, width_b = Poro_v_params[0:3]

            print("phi_ratio = " + str(phi_ratio))
            print("width_a   = " + str(width_a))
            print("width_b   = " + str(width_b))

            poro_v = build_v_layers(z, phi_ratio, width_a, width_b)


    # Set Z_layers for output later
    Z_layers = np.reshape(z, [1,1,nz])


    # Combine horizontal and vertical structures into the full porosity array
    poro = build_full_poro_array(nx, ny, nz, poro_h, poro_v)

    print(f"{np.min(poro) = }, {np.max(poro) = }")

    # -- Permeability -----------------------------------------------------------------------------

    # Build the permeability from the porosity using the Kozeny-Carman relation

    perm = np.zeros((nx, ny, nz))
    perm[:, :, :] = (poro**3) * (1.0 - poro) ** (-2)

    # Rescale so that the mean value is k_scale
    perm_nd_mean = np.mean(perm)
    perm = (k_scale / perm_nd_mean) * perm

    # -- Value scales for plotting ----------------------------------------------------------------

    poro_min = np.amin(poro)
    poro_max = np.amax(poro)

    Log10_perm_min = np.log10(np.amin(perm))
    Log10_perm_max = np.log10(np.amax(perm))

    z_min = np.min(ceil_topo)
    z_max = np.max(base_topo)

    # -- Save data --------------------------------------------------------------------------------

    write_array_to_file(poro, input_folder, "porosity")

    write_array_to_file(perm, input_folder, "permeability")

    write_array_to_file(ceil_topo, input_folder, "ceil_topo")

    write_array_to_file(base_topo, input_folder, "base_topo")

    write_array_to_file(Z_layers, input_folder, "z_layers")


    # -- Generate default XML input file ----------------------------------------------------------

    if (XML_flag):
        write_default_XML(nx,ny,nz,dx,dy,input_folder)

    # -- Plot -------------------------------------------------------------------------------------

    if Plot_flag:

        # Convert injection location index pairs into spatial locations for plotting
        inj_grid_vals = np.zeros((n_inj_locs, 2))

        if n_inj_locs == 1:
            inj_grid_vals[0, 0] = X_grid[int(inj_loc_idxs[0])]
            inj_grid_vals[0, 1] = Y_grid[int(inj_loc_idxs[1])]
        else:
            for k in range(0, n_inj_locs):
                inj_grid_vals[k, 0] = X_grid[int(inj_loc_idxs[k, 0])]
                inj_grid_vals[k, 1] = Y_grid[int(inj_loc_idxs[k, 1])]

        # Specify the colour maps for the permeability and porosity, and
        # the numerical ranges to use for them
        colmap_perm = cm.viridis
        colmap_poro = cm.plasma  # cm.gray

        norm_perm = colors.Normalize(np.log10(np.amin(perm)), np.log10(np.amax(perm)))
        norm_poro = colors.Normalize(0.0, 1.0)

        # Format string to be used for tick labels in colourbars
        cbar_fmt = lambda x, pos: "{:.1f}".format(x)

        n_layers = 7
        layer_indices = evenly_spaced_layers(nz, n_layers)

        fig = plt.figure(figsize=(12, 6), constrained_layout=True)
        gs = GridSpec(2, 2, figure=fig)

        ###########################################################################################

        ax_c = fig.add_subplot(gs[0, 0])

        # surface plot of Porosity
        im_c = ax_c.imshow(
            np.flipud(poro_h[:, :]),
            cmap=colmap_poro,
            norm=norm_poro,
            extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        )

        ax_c.set_title("Horizontal Porosity Structure")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

        # contours of ceiling topography
        cont_H0 = ax_c.contour(
            X_grid, Y_grid, np.transpose(ceil_topo), colors="gray", alpha=0.25
        )

        # Plot injection locations
        for k in range(0, n_inj_locs):
            ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

        m = cm.ScalarMappable(cmap=colmap_poro, norm=norm_poro)
        m.set_array(ceil_topo)

        axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad=-6)

        fig.colorbar(
            m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt)
        )

        ax_c.clabel(cont_H0, inline=True, fmt="H0 = %.2g", fontsize=10)

        ###########################################################################################

        ax_c = fig.add_subplot(gs[0, 1])

        # Vertical structure of porosity
        ax_c.plot(poro_v, z)
        # im_c = ax_c.imshow(
        #     np.flipud( poro_h[:,:] ),
        #     cmap=colmap_poro,
        #     norm=norm_poro,
        #     extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
        # )
        ax_c.set_xlim([0, 1.1])

        ax_c.invert_yaxis()

        ax_c.set_title("Vertical Porosity Structure")
        ax_c.set_xlabel(r"$\phi(z)$", fontsize=20)
        ax_c.set_ylabel(r"$z$", fontsize=20, rotation=0)

        ###########################################################################################

        ax_c = fig.add_subplot(gs[1, 0], projection="3d")

        # # Plot injection locations
        # for k in range(0, n_inj_locs):
        #     ax_c.plot(
        #         [inj_grid_vals[k, 0], inj_grid_vals[k, 0]],
        #         [inj_grid_vals[k, 1], inj_grid_vals[k, 1]],
        #         [z_min, z_max],
        #         "m",
        #         linewidth=2
        #     )

        for k in layer_indices:
            fcolors_poro = generate_colourmap(poro[:, :, k], colmap_poro, 0.0, 1.0)

            ax_c.plot_surface(
                X,
                Y,
                ceil_topo + z[k] * (base_topo - ceil_topo),
                edgecolor="none",
                facecolors=fcolors_poro,
                vmin=0.0,
                vmax=1.0,
                shade=False,
            )

            # Plot injection locations
            if n_inj_locs == 1:
                z_val = ceil_topo[int(inj_loc_idxs[0]), int(inj_loc_idxs[1])] + z[k] * (
                    base_topo[int(inj_loc_idxs[0]), int(inj_loc_idxs[1])]
                    - ceil_topo[int(inj_loc_idxs[0]), int(inj_loc_idxs[1])]
                )
                ax_c.plot([inj_grid_vals[0, 0]], [inj_grid_vals[0, 1]], [z_val], "kx")
            else:
                for l in range(0, n_inj_locs):
                    z_val = ceil_topo[
                        int(inj_loc_idxs[l, 0]), int(inj_loc_idxs[l, 1])
                    ] + z[k] * (
                        base_topo[int(inj_loc_idxs[l, 0]), int(inj_loc_idxs[l, 1])]
                        - ceil_topo[int(inj_loc_idxs[l, 0]), int(inj_loc_idxs[l, 1])]
                    )
                    ax_c.plot(
                        [inj_grid_vals[l, 0]],
                        [inj_grid_vals[l, 1]],
                        [0.95 * z_val],
                        "kx",
                    )

        ax_c.set_title("Porosity")
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)
        ax_c.set_zlabel(r"$z$", fontsize=20, rotation=0)

        ax_c.invert_zaxis()

        m = cm.ScalarMappable(cmap=colmap_poro, norm=norm_poro)
        m.set_array(ceil_topo)

        axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad=-6)

        fig.colorbar(
            m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt)
        )

        ###########################################################################################

        ax_c = fig.add_subplot(gs[1, 1], projection="3d")

        for k in layer_indices:
            fcolors_perm = generate_colourmap(
                np.log10(perm[:, :, k]), colmap_perm, Log10_perm_min, Log10_perm_max
            )

            ax_c.plot_surface(
                X,
                Y,
                ceil_topo + z[k] * (base_topo - ceil_topo),
                edgecolor="none",
                facecolors=fcolors_perm,
                vmin=Log10_perm_min,
                vmax=Log10_perm_max,
                shade=False,
            )

        ax_c.set_title(
            r"$\log_{10}(Permeability)$",
        )
        ax_c.set_xlabel(r"$x$", fontsize=20)
        ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)
        ax_c.set_zlabel(r"$z$", fontsize=20, rotation=0)

        ax_c.invert_zaxis()

        m = cm.ScalarMappable(cmap=colmap_perm, norm=norm_perm)
        m.set_array([np.log10(np.squeeze(perm[:, :, 0]))])

        axins = inset_axes(ax_c, width="5%", height="100%", loc="right", borderpad=-6)

        fig.colorbar(
            m, cax=axins, orientation="vertical", format=FuncFormatter(cbar_fmt)
        )

        ###########################################################################################

        # plt.tight_layout(pad=0.4, w_pad=1.5, h_pad=3.0)
        # plt.tight_layout()

        plt.show()


if __name__ == "__main__":
    main()
