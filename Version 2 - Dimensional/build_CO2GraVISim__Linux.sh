#!/bin/bash

# Exit on any error
set -e


# Program output Directory
program_dir="."


# Compiler
Compiler="gfortran"

# Selection of Compiler flags to choose from
# 1) Flags for debugging
# CompilerFlags="-g -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -fopenmp"
# 2) Flags for an optimised run
CompilerFlags="-O3 -fopenmp"
# 3) Flags for a vanilla run
# CompilerFlags=""

# For the ifx compiler, '-fopenmp' needs to be replaced by '-Qopenmp'

echo
echo
echo "========================================================="
echo "Compiler"
echo "---------------------------------------------------------"
echo
echo "- Compiler"
echo "$Compiler"
echo
echo "- Compiler flags"
echo "$CompilerFlags"
echo
echo "========================================================="



# Source Files Directory
src_dir="./src_files"

# Build Directory
build_dir="$src_dir/build"
mkdir -p "$build_dir"

# .Mod Directory
dot_mod_dir="$build_dir"

# FoX Directory
FoX_dir="../../FoX/fox-master/objs"
#FoX_dir=""

#Check that a value has been set for FoX_dir
if [[ -z "$FoX_dir" ]]; then
    #FoX_dir is empty
    echo "The FoX directory variable FoX_dir has not been set." >&2
    echo "FoX must first be compiled and FoX_dir set to the objs folder (e.g. /FoX/fox-master/objs) so that CO2GraVISim can compile correctly." >&2
    exit 1
fi
	

# Libraries
# lib_dir="$src_dir/libraries"
lib_dir="$lib_dir -L$FoX_dir/lib"
lib_list=(
    FoX_wkml
    FoX_dom
    FoX_sax
    FoX_wcml
    FoX_wxml
    FoX_common
    FoX_utils
    FoX_fsys
)

# Build the library list for the linker
lib_link_flag="$lib_dir"
for lib in "${lib_list[@]}"; do
    lib_link_flag="$lib_link_flag -l$lib"
done


# Incorporate FoX Files
CompilerFlags="$CompilerFlags -I$FoX_dir/finclude"

# List of module files to compile - these need to be in order of increasing dependency
modules_list=(
    CO2GraVISim_global
    CO2GraVISim_XML_input
    CO2GraVISim_timer
    CO2GraVISim_err_calc
    CO2GraVISim_input_parameters
    CO2GraVISim_logging
    CO2GraVISim_vertical_structure
    CO2GraVISim_iter_calculations
    CO2GraVISim_h_routines
    CO2GraVISim_h_res_routines
    CO2GraVISim_pressure_routines
    CO2GraVISim_solver
)

# List of programs to compile
programs_list=(
    CO2GraVISim_single_run
)

echo
echo
echo "========================================================="
echo "Directories and Files"
echo "---------------------------------------------------------"
echo
echo "- Source Directory"
echo "$src_dir"
echo
echo "- Build Directory"
echo "$build_dir"
echo
echo "- .MOD File Directory"
echo "$dot_mod_dir"
echo
echo "- Library Directory"
echo "$lib_dir"
echo
echo "- Program Output Directory"
echo "$program_dir"
echo
echo "- Modules"
for module in "${modules_list[@]}"; do
    echo "$module"
done
echo
echo "- Programs"
for program in "${programs_list[@]}"; do
    echo "$program"
done
echo
echo "- Libraries"
for library in "${lib_list[@]}"; do
    echo "$library"
done
echo
echo "========================================================="

# Compile Fortran modules
echo
echo
echo "========================================================="
echo "Compiling modules"
echo "---------------------------------------------------------"
for module in "${modules_list[@]}"; do
    echo
    echo "-- Compiling $module.f90"

    $Compiler $CompilerFlags -J"$dot_mod_dir" -c "$src_dir/$module.f90" -o "$build_dir/$module.o"
done
echo
echo "========================================================="

# Compile Fortran programs
echo
echo
echo "========================================================="
echo "Compiling programs"
echo "---------------------------------------------------------"
for program in "${programs_list[@]}"; do
    echo
    echo "-- Compiling $program.f90"

    $Compiler $CompilerFlags -J"$dot_mod_dir" -c "$src_dir/$program.f90" -o "$build_dir/$program.o" $lib_link_flag
done
echo
echo "========================================================="

# Link object files
echo
echo
echo "========================================================="
echo "Linking Object Files"
echo "---------------------------------------------------------"
for program in "${programs_list[@]}"; do
    echo
    echo "-- Creating $program"

    $Compiler $CompilerFlags -o "$program_dir/$program" "$build_dir"/*.o $lib_link_flag
done
echo
echo "========================================================="

echo
echo
echo "--- Build process complete ---"
echo
echo
