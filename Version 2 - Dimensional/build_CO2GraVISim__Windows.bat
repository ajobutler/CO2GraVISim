:: build_CO2GraVISim - this batch file compiles all of the executables 
:: involved in the CO2GraVISim program

@Echo OFF
setlocal enabledelayedexpansion

:: Program output Directory
set "program_dir=."


:: Compiler
set Compiler=gfortran

:: Selection of Compiler flags to choose from
:: 1) Flags for debugging
rem set CompilerFlags=-g -O0 -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace -fopenmp
:: 2) Flags for an optimised run
set CompilerFlags=-O3 -fopenmp
:: 3) Flags for a vanilla run
rem set CompilerFlags=

:: For the ifx compiler, '-fopenmp' needs to be replaced by '-Qopenmp'

echo.
echo.
echo =========================================================
echo Compiler
echo ---------------------------------------------------------
echo.
echo - Compiler
echo %Compiler%
echo.
echo - Compiler flags
echo %CompilerFlags%
echo.
echo =========================================================



:: Source Files Directory
set "src_dir=.\src_files"

:: Build Directory
set "build_dir=%src_dir%\build"
if not exist %build_dir% mkdir %build_dir%

:: .Mod Directory
set "dot_mod_dir=%build_dir%"


if "%Compiler%"=="ifx" (
    set dot_mod_flag=-module:"%dot_mod_dir%"
    :: For Intel's ifx compiler, we have to include some compiler-specific libraries. Modify the below as appropriate.
    rem set LIB="C:\Program Files (x86)\Intel\oneAPI\compiler\2023.1.0\lib; C:\Program Files (x86)\Intel\oneAPI\compiler\2023.1.0\windows\compiler\lib; C:\Program Files (x86)\Intel\oneAPI\compiler\2023.1.0\windows\lib"
) else (
    set dot_mod_flag=-J"%dot_mod_dir%"
)
echo %LIB%

:: FoX Directory
set "FoX_dir=.\src_files\FoX_files"

:: Libraries
set "lib_dirs=-L%src_dir%\libraries"
set "lib_dirs=%lib_dirs% -L%FoX_dir%\lib"


set "lib_list="
set "lib_list=%lib_list% -lFoX_wkml"
set "lib_list=%lib_list% -lFoX_dom"
set "lib_list=%lib_list% -lFoX_sax"
set "lib_list=%lib_list% -lFoX_wcml"
set "lib_list=%lib_list% -lFoX_wxml"
set "lib_list=%lib_list% -lFoX_common"
set "lib_list=%lib_list% -lFoX_utils"
set "lib_list=%lib_list% -lFoX_fsys"

set "lib_link_flag=%lib_dirs% %lib_list%"




::Incorporate FoX Files
set CompilerFlags=%CompilerFlags% -I%FoX_dir%\finclude


:: List of module files to compile - these need to be in order
:: of increasing dependency
set modules_list=CO2GraVISim_global 
set modules_list=%modules_list% CO2GraVISim_XML_input
set modules_list=%modules_list% CO2GraVISim_timer
set modules_list=%modules_list% CO2GraVISim_err_calc
set modules_list=%modules_list% CO2GraVISim_input_parameters 
set modules_list=%modules_list% CO2GraVISim_logging
set modules_list=%modules_list% CO2GraVISim_vertical_structure 
set modules_list=%modules_list% CO2GraVISim_iter_calculations
set modules_list=%modules_list% CO2GraVISim_h_routines
set modules_list=%modules_list% CO2GraVISim_h_res_routines
set modules_list=%modules_list% CO2GraVISim_pressure_routines
set modules_list=%modules_list% CO2GraVISim_solver 

:: List of programs to compile
set programs_list=CO2GraVISim_single_run

echo.
echo.
echo =========================================================
echo Directories and Files
echo ---------------------------------------------------------
echo.
echo - Source Directory
echo %src_dir%
echo.
echo - Build Directory
echo %build_dir%
echo.
echo - .MOD File Directory
echo %dot_mod_dir%
echo.
echo - Library Directory
echo %lib_dirs%
echo.
echo - Program Output Directory
echo %program_dir%
echo.
echo - Modules
for %%m in (%modules_list%) do (
    set "module=%%m"
    echo !module!
)
echo.
echo - Programs
for %%p in (%programs_list%) do (
    set "program=%%p"
    echo !program!
)
echo.
echo - Libraries
for %%l in (%lib_list%) do (
    set "library=%%l"
    echo !library!
)
echo.
echo =========================================================


:: Compile Fortran modules
echo.
echo.
echo =========================================================
echo Compiling modules
echo ---------------------------------------------------------
for %%m in (%modules_list%) do (
    set "module=%%m"
    echo.
    echo -- Compiling !module!.f90

    %Compiler% %CompilerFlags% %dot_mod_flag% -c "%src_dir%\!module!.f90" -o "%build_dir%\!module!.o"
    if errorlevel 1 (
        echo Compilation failed for !module!.f90
        exit /b 1
    )
)
echo.
echo =========================================================


:: Compile Fortran programs
echo.
echo.
echo =========================================================
echo Compiling programs
echo ---------------------------------------------------------
for %%p in (%programs_list%) do (
    set "program=%%p"
    echo.
    echo -- Compiling !program!.f90

    %Compiler% %CompilerFlags% %dot_mod_flag% -c "%src_dir%\!program!.f90" -o "%build_dir%\!program!.o" %lib_link_flag%
    if errorlevel 1 (
        echo Compilation failed for !program!.f90
        exit /b 1
    )
)
echo.
echo =========================================================


:: Link object files
echo.
echo.
echo =========================================================
echo Linking Object Files
echo ---------------------------------------------------------
for %%p in (%programs_list%) do (
    set "program=%%p"
    echo.
    echo -- Creating !program!

    %Compiler% %CompilerFlags% -o "%program_dir%\!program!.exe" %build_dir%\*.o %lib_link_flag%

    if errorlevel 1 (
        echo Compilation failed for !program!
        exit /b 1
    )
)
echo.
echo =========================================================

echo.
echo.
echo --- Build process complete ---
echo.
echo.

