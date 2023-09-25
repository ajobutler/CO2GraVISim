# Slices_plot.py
# This script produces plots of the active current thickness h,
# trapped current thickness h_res, and ambient pressure P. 
# The plots produced are the outline as if viewed either along the
# x or y direction (or an equivalent outline in the case of P).
# Along with plots of the ceiling topography, the porosity, and 
# the permeability, for context.


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors


# Read in plot times
plot_times = np.loadtxt("./Output/Other/plot_times.txt")


# load parameter values
parameters = np.loadtxt(
    "./Output/Other/parameter_values.txt"
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
inj_locs_data = np.loadtxt("./Output/Other/injection_locations.txt")

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


# Read in Topography, Porosity, and Permeability
H0              = np.loadtxt("./Input/ceil_topo.txt")
B0              = np.loadtxt("./Input/base_topo.txt")
Porosity        = np.loadtxt("./Input/porosity.txt")
Permeability    = np.loadtxt("./Input/permeability.txt")

# Specify the colour maps for the permeability and porosity, and 
# the numerical ranges to use for them
colmap_perm = cm.gray
colmap_poro = cm.gray

norm_perm = colors.Normalize(0., np.max(Permeability))
norm_poro = colors.Normalize(0., np.max(Porosity))


# Load plot data into arrays
h_array         = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])
h_res_array     = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])
V_active_array  = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])
P_array         = np.zeros([np.shape(H0)[0], np.shape(H0)[1], len(plot_times)])


# Load data from each of the output times and fill in the corresponding arrays
for i, t in enumerate(plot_times):
    # load height data for this frame
    if i == 0:
        h_array[:, :, i] = np.zeros([np.shape(H0)[0], np.shape(H0)[1]])
    else:
        h_array[:, :, i] = np.loadtxt("./Output/Current_Thickness/h" + "{0:02d}".format(i) + ".txt")

    V_active_array[:, :, i] = Porosity * (1.0 - s_a_i) * h_array[:, :, i]
    h_res_array[:, :, i] = np.loadtxt(
        "./Output/Current_Thickness/h_res" + "{0:02d}".format(i) + ".txt"
    )
    P_array[:, :, i] = np.loadtxt(
        "./Output/Current_Pressure/P" + "{0:02d}".format(i) + ".txt"
    )


# Maximum Current Thickness and Current Volume
max_current_thickness = np.amax(h_array)
print("Maximum current thickness is " + str(max_current_thickness))

max_h_res = np.max(h_res_array)
print("Maximum trapped current thickness is " + str(max_h_res))

max_V_active = np.amax(V_active_array)
print("Maximum current volume is " + str(max_V_active))


max_amb_pressure = np.max(P_array)
print("Maximum ambient pressure is " + str(max_amb_pressure))
min_amb_pressure = np.min(P_array)
print("Minimum ambient pressure is " + str(min_amb_pressure))

max_abs_amb_pressure = np.max(np.abs(P_array))

# Formatting string for the contour labels
cbar_fmt = lambda x, pos: "{:.2f}".format(x)


## Main loop over each of the output plots
for i, t in enumerate(plot_times):
    print("Plot " + str(i))

    # load height data for this frame
    h           = h_array[:, :, i]
    h_res       = h_res_array[:, :, i]
    V_active    = V_active_array[:, :, i]
    P           = P_array[:, :, i]

    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(3, 3, figure=fig)

    str_title = "[{:02d} of {:02d}]: t = {:3.2e}".format(i, len(plot_times) - 1, t)

    fig.suptitle(
        str_title,
        fontsize=16,
        x=0.01,
        y=0.01,
        horizontalalignment="left",
        verticalalignment="bottom",
        bbox=dict(facecolor="none", edgecolor="black"),
    )

    # -- Topography ----------------------------------------------------------------
    ax_c = fig.add_subplot(gs[0, 0])

    #Contours of ceiling topography
    ax_c.contour(X_grid, Y_grid, H0, colors="gray", alpha=0.25, levels=4)

    # Plot injection locations
    for k in range(0, n_inj_locs):
        ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

    ax_c.set_title("Ceiling Topography")
    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

    # -- Porosity ------------------------------------------------------------------
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

    # -- Permeability ------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[0, 2])

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

    # -- Current Thickness (outline along x) -------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])

    #Plot outline of mobile current (h) and trapped current (h+h_res)
    ax_c.plot(X_grid, np.max(h, 0), color="blue")
    ax_c.plot(X_grid, np.max(h + h_res, 0), color="red", linestyle="dashed")
    ax_c.set_ylim(0.0, 1.0)
    ax_c.invert_yaxis()

    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$\max_{y}\{h(x,y,t)\}$", fontsize=20)

    # -- Current Volume (outline along x) ---------------------------------------------
    ax_c = fig.add_subplot(gs[1, 1])

    #Plot equivalent outline of mobile and trapped current volumes
    ax_c.plot(X_grid, np.max(V_active, 0), color="blue")
    ax_c.set_ylim(0.0, 1.0)
    ax_c.invert_yaxis()

    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$\max_{y}\{V(x,y,t)\}$", fontsize=20)

    # -- Pressure (outline along x) ---------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 2])

    #Plot equivalent outline of ambient pressure
    ax_c.plot(X_grid, np.max(P, 0), color="blue")
    ax_c.set_ylim(min_amb_pressure, max_amb_pressure)

    ax_c.set_xlabel(r"$x$", fontsize=20)
    ax_c.set_ylabel(r"$\max_{y}\{P(x,y,t)\}$", fontsize=20)

    # -- Current Thickness (outline along y) -----------------------------------------
    ax_c = fig.add_subplot(gs[2, 0])

    #Plot outline of mobile current (h) and trapped current (h+h_res)
    ax_c.plot(Y_grid, np.max(h, 1), color="blue")
    ax_c.plot(Y_grid, np.max(h + h_res, 1), color="red", linestyle="dashed")
    ax_c.set_ylim(0.0, 1.0)
    ax_c.invert_yaxis()

    ax_c.set_xlabel(r"$y$", fontsize=20)
    ax_c.set_ylabel(r"$\max_{x}\{h(x,y,t)\}$", fontsize=20)

    # -- Current Volume (outline along y) --------------------------------------------
    ax_c = fig.add_subplot(gs[2, 1])

    #Plot equivalent outline of mobile and trapped current volumes
    ax_c.plot(Y_grid, np.max(V_active, 1), color="blue")
    ax_c.set_ylim(0.0, 1.0)
    ax_c.invert_yaxis()

    ax_c.set_xlabel(r"$y$", fontsize=20)
    ax_c.set_ylabel(r"$\max_{x}\{V(x,y,t)\}$", fontsize=20)

    # -- Pressure (outline along y) -------------------------------------------------
    ax_c = fig.add_subplot(gs[2, 2])

    #Plot equivalent outline of ambient pressure
    ax_c.plot(Y_grid, np.max(P, 1), color="blue")
    ax_c.set_ylim(min_amb_pressure, max_amb_pressure)

    ax_c.set_xlabel(r"$y$", fontsize=20)
    ax_c.set_ylabel(r"$\max_{x}\{P(x,y,t)\}$", fontsize=20)


    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    # Save this frame
    plt.savefig("./plots/Slices_plot/temp/t" + "{0:02d}".format(i) + ".png")
    plt.close()
