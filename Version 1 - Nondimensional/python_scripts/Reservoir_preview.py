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


# load grid parameter values
with open("./Input/grid_parameters.txt") as f:
    # Skip the initial comment lines
    lines = (line for line in f if not line.startswith("--"))
    parameters = np.loadtxt(lines)
    # nx ny dx dy
nx = int(parameters[0])
ny = int(parameters[1])
dx = parameters[2]
dy = parameters[3]


# Spatial grids
X_grid = np.arange(-(nx - 1) / 2.0, (nx - 1) / 2.0 + 1.0) * dx
Y_grid = np.arange(-(ny - 1) / 2.0, (ny - 1) / 2.0 + 1.0) * dy

X, Y = np.meshgrid(X_grid, Y_grid)

# load injection locations
with open("./Input/injection_locations.txt") as f:
    # Skip the initial comment lines
    lines = (line for line in f if not line.startswith("--"))
    # lines now consists of the number of injection locations,
    # followed by the index pairs of their locations
    lines_list = list(lines)
    n_inj_locs = int(lines_list[0].strip())
    inj_locs_data = np.loadtxt(lines_list[1:])

# Convert these index pairs into spatial locations for plotting
inj_grid_vals = np.zeros((n_inj_locs, 2))

if n_inj_locs == 1:
    inj_grid_vals[0, 0] = X_grid[int(inj_locs_data[0])]
    inj_grid_vals[0, 1] = Y_grid[int(inj_locs_data[1])]
else:
    for k in range(0, n_inj_locs):
        inj_grid_vals[k, 0] = X_grid[int(inj_locs_data[k, 0])]
        inj_grid_vals[k, 1] = Y_grid[int(inj_locs_data[k, 1])]


# Read in Topography, Porosity, and Permeability
H0 = np.loadtxt("./Input/ceil_topo.txt")
B0 = np.loadtxt("./Input/base_topo.txt")
Porosity = np.loadtxt("./Input/porosity.txt")
Permeability = np.loadtxt("./Input/permeability.txt")


# Specify the colour maps for the permeability and porosity, and
# the numerical ranges to use for them
colmap_perm = cm.gray
colmap_poro = cm.gray

norm_perm = colors.Normalize(0.0, np.max(Permeability))
norm_poro = colors.Normalize(0.0, np.max(Porosity))

# Format string to be used for tick labels in colourbars
cbar_fmt = lambda x, pos: "{:.1f}".format(x)


fig = plt.figure(figsize=(10, 6))
gs = GridSpec(2, 2, figure=fig)


# -- Ceiling Topography -----------------------------------------------------------------
ax_c = fig.add_subplot(gs[0, 0])

# Contours of ceiling topography
cont_H0 = ax_c.contour(X_grid, Y_grid, H0, colors="gray", alpha=0.25, levels=10)
ax_c.imshow(
    np.flipud(H0),
    cmap=cm.seismic,
    norm=colors.Normalize(np.min(H0), np.max(H0)),
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


# -- Basement Topography ----------------------------------------------------------------
ax_c = fig.add_subplot(gs[1, 0])

# Contours of basement topography
cont_B0 = ax_c.contour(X_grid, Y_grid, B0, colors="gray", alpha=0.25, levels=10)
ax_c.imshow(
    np.flipud(B0),
    cmap=cm.seismic,
    norm=colors.Normalize(np.min(B0), np.max(B0)),
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


# -- Porosity --------------------------------------------------------------------------
ax_c = fig.add_subplot(gs[0, 1])

# Surface plot of the Porosity array
im_c = ax_c.imshow(
    np.flipud(Porosity),
    cmap=colmap_poro,
    norm=norm_poro,
    extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
)

# Plot the injection locations
for k in range(0, n_inj_locs):
    ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")


ax_c.set_title("Porosity")
ax_c.set_xlabel(r"$x$", fontsize=20)
ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

# Add a colourbar, and make it the same height as the surface plot
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

# -- Permeability ----------------------------------------------------------------
ax_c = fig.add_subplot(gs[1, 1])

# Surface plot of the Permeability array
im_c = ax_c.imshow(
    np.flipud(Permeability),
    cmap=colmap_perm,
    norm=norm_perm,
    extent=[X_grid[0], X_grid[-1], Y_grid[0], Y_grid[-1]],
)

# Plot the injection locations
for k in range(0, n_inj_locs):
    ax_c.plot(inj_grid_vals[k, 0], inj_grid_vals[k, 1], "mx")

ax_c.set_title("Permeability")
ax_c.set_xlabel(r"$x$", fontsize=20)
ax_c.set_ylabel(r"$y$", fontsize=20, rotation=0)

# Add a colourbar, and make it the same height as the surface plot
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


# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

# Save this frame
plt.savefig("./plots/Reservoir_preview/Reservoir_preview.png")

plt.show()


# plt.close()
