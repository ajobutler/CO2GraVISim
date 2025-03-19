:: CO2GraVISim_batch_manage (26/06/24)
:: This script organises and runs several instances of CO2GraVISim_batch_run
:: Built with the help of ChatGPT!


@echo off
setlocal enabledelayedexpansion

set "baseFolder=.\batch_runs\run_"


:: Set the number of sets of runs that we want to perform at the same time
:: Check for command-line inputs for the number of instances
if "%~1"=="" (
    set "instances=4"
) else (
    set "instances=%1"
)

:: We first need to find the total number of runs to perform, by iterating through
:: the run folders that are present

:: Initialize NumRuns to 0
set /a "NumRuns=0"

:: Scan directories to find the highest number N in /run_N/
for /d %%d in (%baseFolder%*) do (

    set "folderName=%%~nxd"
    set "number=!folderName:*_=!"
   
    :: Convert to number for comparison
    set /a "num=!number!"
    if !num! GTR !NumRuns! (
        set /a "NumRuns=!num!"
    )
)

:: Make sure that instances <= NumRuns
if %instances% gtr %NumRuns% (
    set instances=%NumRuns%
)

:: Calculate the typical size of each chunk
set /a ChunkSize=NumRuns/instances
set /a remainder=NumRuns%%instances

:: Print a summary of the batch distribution process
echo:
echo -----------------------------------------------------------
echo Splitting %NumRuns% runs into %instances% groups of %ChunkSize%, with remainder %remainder%
echo -----------------------------------------------------------
echo:

:: Initialise
set start_index=1
set stop_index=0

for /L %%i in (1,1,%instances%) do (

    set /a stop_index=!start_index!+%ChunkSize%-1

    @REM If there is remainder still to distribute, add one to the stop index
    @REM and remove it from the remainder
    if !remainder! gtr 0 (
        set /a stop_index+=1
        set /a remainder-=1
    )   

    @REM Make sure that we don't overshoot (should be covered by the above)
    if %%i==%instances% (
        set stop_index=%NumRuns%
    )

    start cmd /k "CO2GraVISim_batch_run.bat !start_index! !stop_index!"


    set /a start_index=!stop_index!+1
)

endlocal