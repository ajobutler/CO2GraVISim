# Volume_plot.py (04/09/24)
# This script plots the volumes of the various CO2 components from the
# results produced by a run of CO2GraVISim

import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
from matplotlib.gridspec import GridSpec


### Deal with inputs for input and output data folders ######################################

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
    help="If set to e.g. ./batch_runs/run_1, this sets Input to ./batch_runs/run_1/Input and Output to ./batch_runs/run_1/Output",
)


# Parse the arguments
args = parser.parse_args()

if args.batch:
    input_folder = Path(args.batch) / "Input"  # f"{args.batch.rstrip(r'/\\\'"')}/Input"
    output_folder = Path(args.batch) / "Output"
else:
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

#############################################################################################


# load parameter values
parameters = np.loadtxt(
    output_folder / "Other/parameter_values.txt"
)  # nx ny dx dy M Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve
nx, ny, nz, dx, dy, M, Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve = parameters
nx, ny, nz = [int(i) for i in [nx, ny, nz]]  # Convert these to integers


# Load the volume data recorded at each step, rather than just at the output times
Volume_data = np.loadtxt(output_folder / "Other/Volumes.txt", skiprows=1)
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


# Time scale and units
if Times[-1] > 730.0:
    # Express in terms of standard years
    Time_scale = 365.25
    Time_unit = "years"
    Times = Times / Time_scale
else:
    Time_scale = 1.0
    Time_unit = "days"

# fig, ax_c = plt.subplots(nrows=3, ncols=1, figsize=(12, 6), constrained_layout=True)

# create a figure
fig = plt.figure(figsize=(12, 6))

# create grid for different subplots
spec = GridSpec(ncols=1, nrows=2, height_ratios=[2, 1])

fig.suptitle("Volumes [$m^{3}$]", fontsize=20)

### Total Volume ##########################################################################

ax_V = fig.add_subplot(spec[0])

# Colours for the bounding contours of the free and trapped CO2 regions
clr_h = "green"
clr_h_res = "tab:brown"

# Total CO2 injected - every step
ax_V.plot(Times, Volumes_injected, "b-", label="Total CO$_{2}$ injected")

# Total mobile volume - every step
ax_V.plot(Times, Volumes_mobile, "-.", color=clr_h, label="Mobile CO$_{2}$")
# Total trapped volume - every step
ax_V.plot(Times, Volumes_trapped, ":", color=clr_h_res, label="Trapped CO$_{2}$")

# Total dissolved volume
ax_V.plot(Times, Volumes_dissolution, "m-", label="Dissolved CO$_{2}$")


# # ax_V.plot(
# #     Times,
# #     Volumes_total,
# #     "c--",
# #     label="Incl. volume lost through boundaries",
# # )

# Volume lost through the boundaries
ax_V.plot(
    Times,
    Volumes_boundary,
    "c--",
    label="Volume lost through boundaries",
)


# Total free-phase CO2 volume - every step
ax_V.plot(
    Times,
    Volumes_mobile + Volumes_trapped,
    "k--",
    label="Free-phase CO$_{2}$ (Mo. + Tr.)",
)
# Total CO2 phase (including dissolved) - every step
ax_V.plot(
    Times,
    # Volumes_mobile + Volumes_trapped + Volumes_dissolution,
    Volumes_total,
    "r:",
    # label="Total CO$_{2}$ (Mo. + Tr. + Ds.)",
    label="Total CO$_{2}$",
)

# # Current output time for this plot
# ax_V.axvline(x=t, color="gray")

# ax_V.set_xlabel(f"$t$ [{Time_unit}]", fontsize=20)
ax_V.set_ylabel(r"$V$", fontsize=20, labelpad=20, rotation=0)

# ax_V.legend(loc="upper left", fontsize=9)
ax_V.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)

# Add a grid to the volume plot
ax_V.grid(True, linestyle="-", alpha=0.5)


### Deviation from target volume #############################################################

ax_dV = fig.add_subplot(spec[1], sharex=ax_V)

# # # V_total = (
# # #     Volumes_total + Volumes_boundary
# # # )  # Volumes_mobile + Volumes_trapped + Volumes_dissolution
V_inj_max = np.max(Volumes_injected)
delta_V = Volumes_total - Volumes_injected

delta_V_max = np.max(np.abs(delta_V))
delta_V_percentage = 100.0 * (delta_V_max / V_inj_max)

ax_dV.plot(
    Times,
    delta_V,
    "r",
    label=r"$\Delta V = V_{Mo.+Tr.+Ds.+Bd.} - V_{inj}$",
)

ax_dV.hlines(
    [-delta_V_max, delta_V_max],
    Times[0],
    Times[-1],
    label=f"$\\|\\Delta V\\|_{{max}} = {delta_V_max:4.3f} = {delta_V_percentage:.2f}\\% \\times V_{{inj}}^{{max}}$",
    colors="k",
    linestyles="--",
)

# ax_dV.set_title(r"$V_{inj} - V_{M+T+D}$")
ax_dV.set_xlabel(f"$t$ [{Time_unit}]", fontsize=20)
ax_dV.set_ylabel(r"$\Delta V$", fontsize=20, labelpad=20, rotation=0)

# ax_dV.legend(loc="best", fontsize=9, framealpha=0.5)
ax_dV.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)

# Add a grid to the volume plot
ax_dV.grid(True, linestyle="-", alpha=0.5)

### Display ##################################################################################

spec.tight_layout(fig, rect=[0, 0, 1, 0.97])

plt.show()
# # # plt.savefig("./plots/Volume.png")
