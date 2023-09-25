:: build_CO2GraVISim - this batch file compiles all of the executables 
:: involved in the CO2GraVISim program

@Echo OFF

:: Compiler
set Compiler=gfortran

:: Flags for debugging
@REM set CompilerFlags=-Og -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
:: Flags for an optimised run
set CompilerFlags=-O3
:: Flags for a normal run
@REM set CompilerFlags=


ECHO Compiling

cd src_files

ECHO -- Compiling CO2GraVISim_global.mod --
%Compiler% %CompilerFlags% -c CO2GraVISim_global.f90

ECHO -- Compiling CO2GraVISim_injection_profiles.mod --
%Compiler% %CompilerFlags% -c CO2GraVISim_injection_profiles.f90

ECHO -- Compiling CO2GraVISim_input_parameters.mod --
%Compiler% %CompilerFlags% -c CO2GraVISim_input_parameters.f90 CO2GraVISim_global.f90

ECHO -- Compiling CO2GraVISim_solver.mod --
%Compiler% %CompilerFlags% -c CO2GraVISim_solver.f90 CO2GraVISim_input_parameters.f90 CO2GraVISim_injection_profiles.f90 CO2GraVISim_global.f90

ECHO -- Compiling CO2GraVISim_InputGen.exe --
%Compiler% %CompilerFlags% -c CO2GraVISim_InputGen.f90 CO2GraVISim_global.f90
%Compiler% %CompilerFlags% CO2GraVISim_InputGen.f90 CO2GraVISim_global.f90 -o ../CO2GraVISim_InputGen.exe

ECHO -- Compiling CO2GraVISim_single_run.exe --
%Compiler% %CompilerFlags% -c CO2GraVISim_single_run.f90 CO2GraVISim_global.f90 CO2GraVISim_solver.f90 CO2GraVISim_input_parameters.f90 CO2GraVISim_injection_profiles.f90
%Compiler% %CompilerFlags% CO2GraVISim_single_run.f90 CO2GraVISim_global.f90 CO2GraVISim_solver.f90 CO2GraVISim_input_parameters.f90 CO2GraVISim_injection_profiles.f90 -o ../CO2GraVISim_single_run.exe

cd ..

ECHO Compiling Complete.