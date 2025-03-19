#!/bin/sh
# build_all - this batch file compiles all of the executables involved in my 
# Confined Current reservoir simulator

### Compiler
Compiler="gfortran"

### Flags for debugging
# CompilerFlags="-Og -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace"
### Flags for an optimised run
CompilerFlags="-O3"
### Flags for a normal run
# CompilerFlags= ""

echo "Compiling"

cd src_files


echo "-- Compiling CO2GraVISim_global.mod --"
${Compiler} ${CompilerFlags} -c CO2GraVISim_global.f90

echo "-- Compiling CO2GraVISim_injection_profiles.mod --"
${Compiler} ${CompilerFlags} -c CO2GraVISim_injection_profiles.f90

echo "-- Compiling CO2GraVISim_input_parameters.mod --"
${Compiler} ${CompilerFlags} -c CO2GraVISim_input_parameters.f90 CO2GraVISim_global.f90

echo "-- Compiling CO2GraVISim_solver.mod --"
${Compiler} ${CompilerFlags} -c CO2GraVISim_solver.f90 CO2GraVISim_input_parameters.f90 CO2GraVISim_injection_profiles.f90 CO2GraVISim_global.f90

echo "-- Compiling CO2GraVISim_InputGen --"
${Compiler} ${CompilerFlags} -c CO2GraVISim_InputGen.f90 CO2GraVISim_global.f90
${Compiler} ${CompilerFlags} CO2GraVISim_InputGen.f90 CO2GraVISim_global.f90 -o ../CO2GraVISim_InputGen

echo "-- Compiling CO2GraVISim_single_run --"
${Compiler} ${CompilerFlags} -c CO2GraVISim_single_run.f90 CO2GraVISim_global.f90 CO2GraVISim_solver.f90 CO2GraVISim_input_parameters.f90 CO2GraVISim_injection_profiles.f90
${Compiler} ${CompilerFlags} CO2GraVISim_single_run.f90 CO2GraVISim_global.f90 CO2GraVISim_solver.f90 CO2GraVISim_input_parameters.f90 CO2GraVISim_injection_profiles.f90 -o ../CO2GraVISim_single_run
cd ..

echo "Compiling Complete."
