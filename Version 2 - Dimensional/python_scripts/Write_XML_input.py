# Write_XML_input.py (14/01/25)
# This python script generates the XML input document for CO2GraVISim for the specified input parameter values

import numpy as np
import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from pathlib import Path


# -------------------------------------------------------------------------------------------------
# -- Functions ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


def generate_XML(
    grid_parameters_dict,
    BC_dict,
    flow_parameters_dict,
    Injection_dict,
    Plot_times_dict,
    XML_filepath,
):
    # This function generates the full XML document from the respective XML blocks and saves it to the specified path

    ## XML document setup
    root = ET.Element("cml", xmlns="http://www.xml-cml.org/schema")

    ## Build relevant XML blocks
    # # build_dataset(root, dataset_dict)
    build_grid_parameters(root, grid_parameters_dict)
    build_BCs(root, BC_dict)
    build_flow_parameters(root, flow_parameters_dict)
    # # build_injection_interval_parameters(root, Injection_dict)
    build_injection_locations(root, Injection_dict)
    build_injection_schedule(root, Injection_dict)
    build_plot_times(root, Plot_times_dict)

    ## File generation comment
    XML_Comment(root, 0, "")
    XML_Comment(
        root,
        0,
        f"Input file for CO2GraVISim generated at {datetime.today().strftime('%Y-%m-%d %H:%M:%S')}",
    )
    XML_Comment(root, 0, "")

    ## Write to XML document with the appropriate declerations at the start of the document
    xml_declaration = '<?xml version="1.0" encoding="UTF-8"?>\n'
    stylesheet_declaration = '<?xml-stylesheet href="http://www.eminerals.org/XSLT/display.xsl" type="text/xsl"?>\n'
    xml_body = prettify_xml(root)

    final_xml = (
        xml_declaration + stylesheet_declaration + xml_body[xml_body.find("<cml") :]
    )

    with open(XML_filepath, "w") as file:
        file.write(final_xml)


def XML_Comment(element, level, text):
    # This function adds the set text as a comment at the specified level
    # and in the given element
    comment = ET.Comment(text)
    element.insert(level, comment)


def XML_index(element, parent):
    # This function finds the index of an XML element within its parent
    # parent = element.getparent()  # Available in lxml; not in xml.etree.ElementTree
    # if parent is None:
    #     parent = root  # Assuming root if getparent is unavailable

    # Find the index of the target element within its parent
    index = list(parent).index(element)

    return index


# # def build_dataset(root, dataset_dict):
# #     """This function modifies root to add in the dataset child element,
# #     generating a block of the form
# #         <dataset>
# #             <filepath>...</filepath>
# #         </dataset>

# #     """
# #     dataset = ET.SubElement(root, "dataset")
# #     filepath = ET.SubElement(dataset, "filepath")
# #     filepath.text = f"{dataset_dict['filepath']}"


def build_grid_parameters(root, grid_parameters_dict):
    """This function modifies root to add in the grid_parameters child element,
    generating a block of the form
        <grid_parameters>
            <nx>...</nx>
        ...
        </grid_parameters>
    """
    grid_parameters = ET.SubElement(root, "grid_parameters")
    for _, (key, val) in enumerate(grid_parameters_dict.items()):
        entry = ET.SubElement(grid_parameters, key)
        entry.text = f"{val}"

    XML_Comment(root, XML_index(grid_parameters, root), " Grid parameters ")
    XML_Comment(root, XML_index(grid_parameters, root), " Units: dx, dy [m] ")


def build_BCs(root, BC_dict):
    """This function modifies root to add in the boundary_conditions child element,
    generating a block of the form
        <boundary_conditions>
            <current_thickness>
                    <north>...</north>
                    ...
            </current_thickness>

            <ambient_pressure>
                    ...
            </ambient_pressure>
        </boundary_conditions>
    """
    boundary_conditions = ET.SubElement(root, "boundary_conditions")

    BC_current = ET.SubElement(boundary_conditions, "current_thickness")
    for _, (key, val) in enumerate(BC_dict["current_thickness"].items()):
        entry = ET.SubElement(BC_current, key)
        entry.text = f"{val}"

    BC_pressure = ET.SubElement(boundary_conditions, "ambient_pressure")
    for _, (key, val) in enumerate(BC_dict["ambient_pressure"].items()):
        entry = ET.SubElement(BC_pressure, key)
        entry.text = f"{val}"

    XML_Comment(root, XML_index(boundary_conditions, root), " Boundary Conditions ")
    XML_Comment(
        root,
        XML_index(boundary_conditions, root),
        " Flags: 1 for Dirichlet BC (f=0), 2 for Neumann BC (df/dn = 0) ",
    )


def build_flow_parameters(root, flow_parameters_dict):
    """This function modifies root to add in the flow_parameters child element
    generating a block of the form
        <flow_parameters>
            <rho_c>...</rho_c>
            ...
        </flow_parameters>
    """

    Comments_dict = {
        "rho_c": " CO2 density [kg m^-3] ",
        "rho_a_unsat": " Unsaturated ambient density [kg m^-3] ",
        "rho_a_sat": " Saturated ambient density [kg m^-3] ",
        "mu_c": " CO2 viscosity [Pa s] ",
        "mu_a": " Ambient viscosity [Pa s] ",
        "s_c_r": " CO2 residual trapping saturation [-] ",
        "s_a_i": " Ambient irreducible saturation [-] ",
        "krn_mobile": " Relative permeability of CO2 in mobile region [-] ",
        "krw_residual": " Relative permeability of ambient in trapping region [-] ",
        "C_sat": " Volume fraction of CO2 in saturated ambient [-] ",
        "g": " Gravitational acceleration [m s^-2] ",
        "D_mol": " Molecular diffusivity of CO2 in ambient [m^2 s^-1] ",
        "perm_ratio": " Ratio of vertical to horizontal absolute permeability [-] ",
        "omega_conv": " Control prefactor for convective dissolution [-] ",
    }

    flow_parameters = ET.SubElement(root, "flow_parameters")
    for _, (key, val) in enumerate(flow_parameters_dict.items()):
        entry = ET.SubElement(flow_parameters, key)
        entry.text = f"{val}"
        XML_Comment(
            flow_parameters, XML_index(entry, flow_parameters), Comments_dict[key]
        )

    XML_Comment(root, XML_index(flow_parameters, root), " Flow Parameters ")


# def build_injection_interval_parameters(root, Injection_dict):
#     """This function modifies root to add in the injection_interval_parameters child element
#     generating a block of the form
#         <injection_interval_parameters>
#             <z_idx_ceil>...</z_idx_ceil>
#             <z_idx_base>...</z_idx_base>
#         </injection_interval_parameters>
#     """

#     injection_interval_parameters = ET.SubElement(root, "injection_interval_parameters")
#     z_idx_ceil = ET.SubElement(injection_interval_parameters, "z_idx_ceil")
#     z_idx_ceil.text = f"{Injection_dict['z_idx_ceil']}"
#     z_idx_base = ET.SubElement(injection_interval_parameters, "z_idx_base")
#     z_idx_base.text = f"{Injection_dict['z_idx_base']}"

#     XML_Comment(
#         root,
#         XML_index(injection_interval_parameters, root),
#         " Injection Interval Parameters ",
#     )


def build_injection_locations(root, Injection_dict):
    """This function modifies root to add in the injection_locations child element
    generating a block of the form
        <injection_locations>
            <location id="loc1">
                <i>...</i>
                <j>...</j>
            </location>
        ...
        </injection_locations>
    """

    injection_locations = ET.SubElement(root, "injection_locations")
    n_loc = len(Injection_dict["locations"])
    for k in range(0, n_loc):
        location = ET.SubElement(injection_locations, "location", id=f"loc{k+1}")
        i = ET.SubElement(location, "i")
        i.text = f"{Injection_dict['locations'][k][0]}"
        j = ET.SubElement(location, "j")
        j.text = f"{Injection_dict['locations'][k][1]}"

    XML_Comment(root, XML_index(injection_locations, root), " Injection Locations ")
    XML_Comment(
        root,
        XML_index(injection_locations, root),
        " i,j indices of each injection location, with id tags used by <injection_schedule> ",
    )


def build_injection_schedule(root, Injection_dict):
    """This function modifies root to add in the injection_schedule child element
    generating a block of the form
        <injection_schedule>
            <schedule>
                <time>0.0</time>
                <location id="loc1">
                    <rate>...</rate>
                </location>

                <location id="loc2">
                    <rate>...</rate>
                </location>
            ...
            </schedule>

            <schedule>
            ...
            </schedule>
            ...
        </injection_schedule>
    """

    injection_schedule = ET.SubElement(root, "injection_schedule")
    n_times = len(Injection_dict["Flux_times"])
    n_loc = len(Injection_dict["locations"])
    for nt in range(0, n_times):
        schedule = ET.SubElement(injection_schedule, "schedule")
        time = ET.SubElement(schedule, "time")
        time.text = f"{Injection_dict['Flux_times'][nt]}"
        for nl in range(0, n_loc):
            location = ET.SubElement(schedule, "location", id=f"loc{nl+1}")
            rate = ET.SubElement(location, "rate")
            rate.text = f"{Injection_dict['Flux_vals'][nl][nt]}"

    XML_Comment(root, XML_index(injection_schedule, root), " Injection Schedule ")
    XML_Comment(
        root,
        XML_index(injection_schedule, root),
        " Times when the injection flux changes, and the corresponding reservoir-level flux values for each injection location as labelled in <injection_locations> ",
    )
    XML_Comment(
        root,
        XML_index(injection_schedule, root),
        " Units: times [days], fluxes [m^3 s^-1] ",
    )


def build_plot_times(root, Plot_times_dict):
    """This function modifies root to add in the plot_times child element
    generating a block of the form
        <plot_times>
            <time>...</time>
            ...
        </plot_times>
    """

    plot_times = ET.SubElement(root, "plot_times")
    plot_times_vals = Plot_times_dict["times"]
    for k, t in enumerate(plot_times_vals):
        time = ET.SubElement(plot_times, "time")
        time.text = f"{t}"

    XML_Comment(root, XML_index(plot_times, root), " Plot Times ")
    XML_Comment(root, XML_index(plot_times, root), " Units: days ")


def prettify_xml(elem):
    # This function improves the aesthetics of the generated XML string: it will now appear on
    # multiple lines rather than one, and there will be a blank line left after each code block
    # https://stackoverflow.com/questions/15418509/python-and-elementtree-writing-one-long-line-of-output

    ## Apply minidom's improvements
    rough_string = ET.tostring(elem, "utf-8")
    # rough_string = ET.tostring(elem, encoding="unicode", method="xml")
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="\t")

    ## Add blank lines after each code block
    # Split into lines and process
    lines = pretty_xml.split("\n")
    result = []
    for i, line in enumerate(lines):
        result.append(line)
        # Add a blank line if this line is a closing tag and the next line is an opening tag
        if line.strip().startswith("</") and i + 1 < len(lines):
            next_line = lines[i + 1].strip()
            if next_line.startswith("<") and not next_line.startswith("</"):
                result.append("")
    return "\n".join(result)


# -------------------------------------------------------------------------------------------------
# -- Main Code ------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------


def main():

    ### Deal with inputs for input and output data folders ######################################

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Input data folder and XML filename")

    # Add argument for input folder with default value
    parser.add_argument(
        "--input", type=str, default="./Input/", help="Input data folder path"
    )

    # Add argument for the XML file name with default value
    parser.add_argument(
        "--name", type=str, default="XML_input_file", help="Name for the XML file"
    )

    # Parse the arguments
    args = parser.parse_args()

    # Extract the input and output folder paths, and remove any trailing quotation marks and slashes
    input_folder = Path(args.input)

    # Ensure that the input folder exists
    if not input_folder.is_dir():
        print(f"Error: The input folder {input_folder} does not exist.")
        exit(1)

    XML_file_name = args.name
    XML_filename_path = Path(XML_file_name)

    #Check if a file extension has been included in the filename
    #This approach is a little crude - using path.with_suffix() might be better.
    #Would then have to worry about double extensions like .tar.gz
    if XML_filename_path.suffix == '.xml':
        #Remove this file extension (then add it back on later)
        XML_file_name = XML_filename_path.stem
    elif XML_filename_path.suffix:
        #There is an extension and it isn't .xml
        raise ValueError(f"Unexpected file extension included: {XML_filename_path.suffix}")
    
    
    ## ----- Filepath for the generated XML document -----
    XML_filepath = input_folder / f"{XML_file_name}.xml"

    ## Parameter values ###########################################################################

    # # ## ----- dataset -----
    # # dataset_dict = {
    # #     "filepath": r"C:\Users\blabl\Documents\Work\Cambridge_Postdoc\Aquifer_Simulator\CO2GraVISim\standalone_dimensional\Input"
    # # }

    ## ----- flow parameters -----
    flow_parameters_dict = {
        "rho_c": 810.0,  # CO2 density [kg m^-3]
        "rho_a_unsat": 1030.0,  # Ambient density [kg m^-3]
        "rho_a_sat": 1042.0,  # Saturated ambient density [kg m^-3]
        "mu_c": 7e-5,  # CO2 density [Pa s]
        "mu_a": 9e-4,  # Ambient density [Pa s]
        "s_c_r": 0.36,  # Residual CO2 saturation in trapping region [-]
        "s_a_i": 0.2,  # Irreducible ambient saturation in mobile CO2 region [-]
        "krn_mobile": 1.0,  # Relative permeability of the non-wetting phase (CO2) in the mobile region [-]
        "krw_residual": 1.0,  # Relative permeability of the wetting phase (ambient) in the trapping region [-]
        "C_sat": 0.04,  # Volume fraction of CO2 dissolved in ambient [-]
        "g": 9.81,  # gravitational acceleration [m s^-2]
        "D_mol": 2e-9,  # Molecular diffusivity of CO2 in ambient [m^2 s^-1]
        "perm_ratio": 0.1,  # Ratio of vertical to horizontal absolute permeability [-]
        "omega_conv": 0.0,  # Control prefactor for convective dissolution term [-]
    }

    ## ----- Grid parameters -----
    grid_parameters_dict = {
        "nx": 200,  # Number of grid points in x
        "ny": 350,  # Number of grid points in y
        "nz": 20,  # Number of grid points in z
        "dx": 2.25,  # Grid spacing in x
        "dy": 2.25,  # Grid spacing in y
    }

    ## ----- Boundary Conditions -----
    # Boundary condition flags for the current thickness (h) and ambient pressure (P)
    # on the 'North', 'East', 'South', and 'West' boundaries of the domain.
    # Flag values:
    #  - 1 : Dirichlet BC f = 0
    #  - 2 : Neumann BC df/dn = 0
    BC_dict = {
        "current_thickness": {"north": 1, "east": 1, "south": 1, "west": 1},
        "ambient_pressure": {"north": 1, "east": 1, "south": 1, "west": 1},
    }

    ## ----- Plot times -----
    # Times at which to save the full solution, measured in days
    # since the start of the simulation
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
            4383.00000,
            4748.00000,
            5113.00000,
        ],
    }

    ## ----- Injection information -----
    # Information related to the injection.
    Injection_dict = {
        # # "z_idx_ceil": 1,  # z index of the upper sealing layer (caprock) of the storage interval
        # # "z_idx_base": 20,  # z index of the lower sealing layer (basement) of the storage interval
        "locations": [
            [120, 175],
            [80, 175],
        ],  # Array of (i,j) pairs for each of the injection wells.
        # NOTE: This has to be an array rather than a list!
        "Flux_times": [
            0.0,
            365.0,
            730.0,
        ],  # List of times [days] at which the injection flux changes
        "Flux_vals": [
            [0.002, 0.00, 0.00],
            [0.00, 0.001, 0.0],
        ],  # Array of the correspondig flux values [m^3 s^-1] at each of these flux times for each of the injection locations.
        # Each row corresponds to an injection location.
        # NOTE: This has to be an array rather than a list!
    }

    ## Generate the XML document ##################################################################
    generate_XML(
        grid_parameters_dict,
        BC_dict,
        flow_parameters_dict,
        Injection_dict,
        Plot_times_dict,
        XML_filepath,
    )

    print(f"Writing {XML_filepath} with default parameter values.")


if __name__ == "__main__":
    main()
