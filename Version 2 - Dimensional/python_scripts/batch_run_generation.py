# batch_run_generation.py (17/03/24)
# This python script can be used to generate the various folders and
# input files used for a batch run of CO2GraVISim

# Built with help from ChatGPT 3.5

import argparse
from pathlib import Path
import shutil
from datetime import datetime
from itertools import product

from brg.brg_BC_parameters import BC_variations, N_BC
from brg.brg_flow_parameters import flow_parameter_variations, N_flow
from brg.brg_grid_parameters import grid_parameter_variations, N_grid
from brg.brg_injection_parameters import injection_variations, N_injection
from brg.brg_plot_times import plot_times_variations, N_plot
from brg.brg_poro_and_perm_parameters import poro_and_perm_variations, N_poro_and_perm
from brg.brg_topography_parameters import topo_variations, N_topo
from brg.brg_initial_profiles import initial_profile_variations, N_initial

from input_generation import input_gen
from Write_XML_input import generate_XML
from initial_profile_generation import generate_initial_profiles

# from Generate_Plot_times import pt_generate


def archive_existing_runs(target_directory):
    # Generate a unique name for the archive directory based on the current datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_directory = target_directory / f"Archived/old_runs_{timestamp}"

    # Create the archive directory if there are any folders to move
    created_archive = False

    if target_directory.is_dir():
        for item in target_directory.iterdir():
            # If the item is a directory and its name starts with 'run_'
            if item.is_dir() and (
                item.name.startswith("run_") or item.name.startswith("Common_files")
            ):
                if not created_archive:
                    # Make the archive directory, and update the flag so this is only done once
                    archive_directory.mkdir(parents=True, exist_ok=True)
                    created_archive = True
                # Move the run_N directory to the archive directory
                # item contains target_directory as well here, so we need to use item.name
                shutil.move(target_directory / item.name, archive_directory / item.name)
                print(f"Moved {item} to {archive_directory}")


def create_common_files_folder(common_files_dir):
    # Create Common Files folder
    if not common_files_dir.exists():
        common_files_dir.mkdir(parents=True)


def create_run_folders(target_directory, num_runs):

    # Make sure the target directory exists, if not, create it
    if not target_directory.exists():
        target_directory.mkdir(parents=True)

    for i in range(1, num_runs + 1):
        # Create the run_N directory
        run_directory = target_directory / f"run_{i}"
        run_directory.mkdir(parents=True, exist_ok=True)

        # Inside each run_N directory, create Input and Output directories
        input_directory = run_directory / "Input"
        output_directory = run_directory / "Output"

        input_directory.mkdir(parents=True, exist_ok=True)
        output_directory.mkdir(parents=True, exist_ok=True)

        # Create folders within Output directories
        (output_directory / "Current_Pressure").mkdir(parents=True, exist_ok=True)
        (output_directory / "Current_Thickness").mkdir(parents=True, exist_ok=True)
        (output_directory / "Current_Volume").mkdir(parents=True, exist_ok=True)
        (output_directory / "Other").mkdir(parents=True, exist_ok=True)

        print(f"Created: {run_directory}\n\tWith subdirectories: Input, Output")


def main():

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Input and Output data folders")

    # Flag for whether to just preview the number of variations for each input type
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Activate to only display the number of variations of each input type",
    )

    # Flag for whether to archive any existing files in the target location
    parser.add_argument(
        "--archive_off",
        action="store_false",
        help="Activate to skip archiving any batch-run folders already present in the target location",
    )

    # Add arguments for input and output folder with default values
    parser.add_argument(
        "--folder",
        type=str,
        default="./batch_runs/",
        help="Folder path for generating the batch runs. Default is ./batch_runs/",
    )

    # Parse the arguments
    args = parser.parse_args()
    Preview_flag = args.preview
    Archive_flag = args.archive_off

    target_dir = Path(args.folder.rstrip("/\\'\""))
    common_files_dir = target_dir / "Common_files"
    #############################################################################################
    #############################################################################################

    # Total number of runs
    # N_runs = N_flow * N_topo * N_poro_perm * N_flux
    N_runs = (
        N_BC
        * N_flow
        * N_grid
        * N_injection
        * N_plot
        * N_poro_and_perm
        * N_topo
        * N_initial
    )

    print("Variations:")
    print(
        f"""
        {N_BC            = }
        {N_flow          = }
        {N_grid          = }
        {N_injection     = }
        {N_plot          = }
        {N_poro_and_perm = }
        {N_topo          = }
        {N_initial       = }
            """
    )

    print(f"Total number of runs: N = {N_runs}\n\n")

    if Preview_flag:

        print("\n-- [PREVIEW] --\n")

    else:
        # Run the actual calculations and file generation

        # Combinations of indices to access Topography parameters and Porosity and Permeability parameters
        index_array = list(
            product(
                range(0, N_BC),
                range(0, N_flow),
                range(0, N_grid),
                range(0, N_injection),
                range(0, N_plot),
                range(0, N_poro_and_perm),
                range(0, N_topo),
                range(0, N_initial),
            )
        )

        ###########################################################################################
        ### Archive any existing folders ##########################################################
        ###########################################################################################
        if Archive_flag:
            # Archive any existing run_N directories
            print("\n\nArchiving existing Files:")
            archive_existing_runs(target_dir)

        ###########################################################################################
        ### Common Files ##########################################################################
        ###########################################################################################
        print("\n\nGenerating Common Files:")

        # Create the common files folder, if it doesn't already exist
        create_common_files_folder(common_files_dir)

        ### Flags indicating which files vary and which are common between all of the runs
        ## - XML
        # We'll create an XML document for each run, since these are relative small and are named after the run

        ## - Initial profiles
        if N_grid * N_initial > 1:
            # Initial profiles vary between runs
            initial_vary_flag = True
        else:
            # Initial profiles are the same for every run: generate once here and send the output to the common_files folder
            initial_vary_flag = False
            print("-- sending initial profiles to common_files folder")

            IP = initial_profile_variations[0]
            h_init_params = [IP["x_c"], IP["y_c"], IP["r_x"], IP["r_y"], IP["h_max"]]

            generate_initial_profiles(
                grid_parameter_variations[0]["nx"],
                grid_parameter_variations[0]["ny"],
                grid_parameter_variations[0]["dx"],
                grid_parameter_variations[0]["dy"],
                h_init_params,
                common_files_dir,
                Plot_flag=False,
            )

        ## - Porosity, Permeability, and Topography
        # If either N_poro_and_perm > 1 or N_topo > 1, then we will have to run input_generation.py as part of the loop below
        if N_grid * N_poro_and_perm * N_topo > 1:
            # Porosity, permeability, and/or topography differ between runs: generate as part of the loop below
            poro_topo_vary_flag = True
        else:
            # Porosity, permeability, and topography are the same for every run: generate once here and send the output
            # to the common_files folder
            poro_topo_vary_flag = False

            print(
                "-- sending topography, porosity, and permeability to common_files folder"
            )
            input_gen(
                common_files_dir,
                grid_parameter_variations[0]["nx"],
                grid_parameter_variations[0]["ny"],
                grid_parameter_variations[0]["nz"],
                grid_parameter_variations[0]["dx"],
                grid_parameter_variations[0]["dy"],
                len(injection_variations[0]["locations"]),
                injection_variations[0]["locations"],
                poro_and_perm_variations[0]["ptype_h"],
                poro_and_perm_variations[0]["pparams_h"],
                poro_and_perm_variations[0]["ptype_v"],
                poro_and_perm_variations[0]["pparams_v"],
                topo_variations[0]["slope"],
                topo_variations[0]["bump"],
                Plot_flag=False,  # Plot_Flag
                XML_flag=False,
            )

        ###########################################################################################
        ## Create new run folders #################################################################
        ###########################################################################################

        ## Create the run directories, each containing the relevant Input and Output folders
        print("\n\n-- Creating folders --")
        create_run_folders(target_dir, N_runs)

        ###########################################################################################
        ## Main Loop ##############################################################################
        ###########################################################################################

        for i in range(0, N_runs):

            print(f"\n [ -- n = {i+1} -- ] \n")

            ## Relevant filepaths for this loop
            run_directory = target_dir / f"run_{i+1}"
            run_Input_directory = run_directory / f"Input"
            XML_filepath = run_Input_directory / f"run_{i+1}.xml"

            ## Relevant indices for each input set for this loop
            (
                BC_idx,
                flow_idx,
                grid_idx,
                injection_idx,
                plot_idx,
                poro_and_perm_idx,
                topo_idx,
                initial_idx,
            ) = index_array[i][:]

            # Shorthands for relevant parameter dictonaries
            BC = BC_variations[BC_idx]
            Fl = flow_parameter_variations[flow_idx]
            Gr = grid_parameter_variations[grid_idx]
            In = injection_variations[injection_idx]
            Pl = plot_times_variations[plot_idx]

            ## Write XML input file
            generate_XML(
                Gr,
                BC,
                Fl,
                In,
                Pl,
                XML_filepath,
            )

            ## Run generate_initial_profiles only if the initial profiles vary between runs
            if initial_vary_flag:
                IP = initial_profile_variations[initial_idx]
                h_init_params = [
                    IP["x_c"],
                    IP["y_c"],
                    IP["r_x"],
                    IP["r_y"],
                    IP["h_max"],
                ]

                print("-- initial profiles")
                generate_initial_profiles(
                    Gr["nx"],
                    Gr["ny"],
                    Gr["dx"],
                    Gr["dy"],
                    h_init_params,
                    run_Input_directory,
                    Plot_flag=False,
                )

            ## Run input_gen only if porosity, permeability, or topography vary between the runs
            if poro_topo_vary_flag:

                Po = poro_and_perm_variations[poro_and_perm_idx]
                To = topo_variations[topo_idx]

                print("-- topography, porosity, and permeability")
                input_gen(
                    run_Input_directory,
                    Gr["nx"],
                    Gr["ny"],
                    Gr["nz"],
                    Gr["dx"],
                    Gr["dy"],
                    len(In["locations"]),
                    In["locations"],
                    Po["ptype_h"],
                    Po["pparams_h"],
                    Po["ptype_v"],
                    Po["pparams_v"],
                    To["slope"],
                    To["bump"],
                    Plot_flag=False,
                    XML_flag=False,
                )


if __name__ == "__main__":
    main()
