:: CO2GraVISim_batch_run (16/03/24)
:: This is performs a batch run of CO2GraVISim over the inputs stored in
:: batch_runs/run_*/Input, with the resulting outputs stored in 
:: batch_runs/run_*/Output

:: Generated with ChatGPT 3.5's help!


@echo off
setlocal enabledelayedexpansion

set "baseFolder=.\batch_runs\run_"
set "processScript=.\CO2GraVISim_single_run.exe"


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

:: Now maxNum contains the highest number found
if %maxNum% gtr 0 (
    for /l %%i in (1, 1, %maxNum%) do (

        echo: 
        echo [ --- Run %%i of %maxNum% --- ]
        echo: 

        :: Relevant Input and Output folders
        set "inputPath=%baseFolder%%%i\Input"
        set "outputPath=%baseFolder%%%i\Output"

        @REM call %processScript% -input !inputPath! -output !outputPath!
        @REM copy ".\Output\Other\Volumes.txt" "!outputPath!\Other\Volumes.txt"

        :: Check if the Input and Output directories exist
        if exist "!inputPath!" (
            if exist "!outputPath!" (
                :: Call the process function/script with Input and Output paths as arguments
                echo %processScript% -input !inputPath! -output !outputPath!
                call %processScript% -input !inputPath! -output !outputPath!

                :: Currently the Volumes values recorded at every successful iteration are stored in
                :: the default output folder. I now need to copy that across to the output folder
                :: for this run
                copy ".\Output\Other\Volumes.txt" "!outputPath!\Other\Volumes.txt"
            )
        )
    )
) else (
    echo No suitable \run_N\ folders found.
)

endlocal
