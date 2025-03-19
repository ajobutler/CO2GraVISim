:: CO2GraVISim_batch_run (16/03/24)
:: This script performs a batch run of CO2GraVISim over the inputs stored in
:: batch_runs/run_*/Input, with the resulting outputs stored in 
:: batch_runs/run_*/Output

:: Generated with ChatGPT 3.5's help!


@echo off
setlocal enabledelayedexpansion



:: Function to get the current date and time
for /f %%a in ('wmic os get localdatetime ^| find "."') do set datetime=%%a
set "start_timestamp=%datetime:~0,4%-%datetime:~4,2%-%datetime:~6,2% %datetime:~8,2%:%datetime:~10,2%:%datetime:~12,2%"
set /a start_time=((1%datetime:~0,4%-10000)*31556952 + (1%datetime:~4,2%-100)*2629746 + (1%datetime:~6,2%-100)*86400 + (1%datetime:~8,2%-100)*3600 + (1%datetime:~10,2%-100)*60 + (1%datetime:~12,2%-100))

:: Print start timestamp
echo Start time: %start_timestamp%



set "baseFolder=.\batch_runs\run_"
set "processScript=.\CO2GraVISim_single_run.exe"

:: Other flags to add when calling processScript
:: For the unconfined version of the solver, add -unconfined
:: set "runFlags="
set "runFlags=-unconfined"


:: We first need to find the total number of runs to perform, by iterating through
:: the run folders that are present

:: Initialize maxNum to 0
set /a "maxNum=0"

:: Scan directories to find the highest number N in /run_N/
for /d %%d in (%baseFolder%*) do (

    set "folderName=%%~nxd"
    set "number=!folderName:*_=!"
   
    :: Convert to number for comparison
    set /a "num=!number!"
    if !num! GTR !maxNum! (
        set /a "maxNum=!num!"
    )
)

:: Check for command-line inputs for the loop indices
if "%~1"=="" (
    set "start_index=1"
) else (
    set "start_index=%1"
)

if "%~2"=="" (
    set "stop_index=%maxNum%"
) else (
    set "stop_index=%2"
)

:: Make sure that start_index <= stop_index
if %start_index% gtr %stop_index% (
    echo Start index is greater than Stop index: [%start_index%, %stop_index%]
    echo Halting.
    Exit /b
)

:: Make sure the endpoint is at most the number of run folders
if %stop_index% gtr %maxNum% (
    set "%stop_index% = %maxNum%"
)

:: Number of runs to perform, and a counter for the loop
set /a "run_total=%stop_index%-%start_index%+1"
set "run_index=1"

:: Write start and stop indices to standard out, for clarity
echo:
echo:
echo --------------------------------------
echo Start Index:              %start_index%
echo Stop Index:               %stop_index%
echo Number of runs to do:     %run_total%
echo Number of runs available: %maxNum%
echo --------------------------------------
echo:
echo:

:: Main Loop over batch-run folders
if %maxNum% gtr 0 (
    for /l %%i in (%start_index%, 1, %stop_index%) do (

        echo: 
        echo [ --- Run !run_index! of %run_total% --- ]
        echo: 

        :: Relevant Input and Output folders
        set "batchPath=%baseFolder%%%i"
        set "XMLPath=!batchPath!\Input\run_%%i.xml"

        :: Check if the Input and Output directories exist
        if exist "!batchPath!" (
                :: Call the process function/script with Input and Output paths as arguments
                echo %processScript% -batch !batchPath! -xml !XMLPath! !runFlags!
                call %processScript% -batch !batchPath! -xml !XMLPath! !runFlags!
        )

        set /a "run_index+=1"
    )
) else (
    echo No suitable \run_N\ folders found.
)


:: Get end time
for /f %%a in ('wmic os get localdatetime ^| find "."') do set datetime=%%a
set "end_timestamp=%datetime:~0,4%-%datetime:~4,2%-%datetime:~6,2% %datetime:~8,2%:%datetime:~10,2%:%datetime:~12,2%"
set /a end_time=((1%datetime:~0,4%-10000)*31556952 + (1%datetime:~4,2%-100)*2629746 + (1%datetime:~6,2%-100)*86400 + (1%datetime:~8,2%-100)*3600 + (1%datetime:~10,2%-100)*60 + (1%datetime:~12,2%-100))



:: Compute elapsed time in seconds
set /a elapsed=end_time-start_time

:: Convert elapsed time to HH:MM:SS format
set /a hours=elapsed / 3600
set /a minutes=(elapsed %% 3600) / 60
set /a seconds=elapsed %% 60

:: Print start and end timestamps and Elapsed time
echo ****************************************************************
echo Start time  : %start_timestamp%
echo End time    : %end_timestamp%
echo Elapsed time: !hours!h !minutes!m !seconds!s
echo ****************************************************************






endlocal
