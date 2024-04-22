#!/bin/bash

# CO2GraVISim_batch_run (22/04/24)
# This is performs a batch run of CO2GraVISim over the inputs stored in
# batch_runs/run_*/Input, with the resulting outputs stored in 
# batch_runs/run_*/Output

# Generated with ChatGPT 3.5's help!

baseFolder="./batch_runs/run_"
processScript="./CO2GraVISim_single_run"


# We first need to find the total number of runs to perform, by iterating through
# the run folders that are present

# Initialize maxNum to 0
maxNum=0

# Scan directories to find the highest number N in /run_N/
for dir in "${baseFolder}"*/; do
    folderName=$(basename "$dir")
    number="${folderName##*_}"

    # Convert to number for comparison
    num=$((10#$number))
    if [ "$num" -gt "$maxNum" ]; then
        maxNum="$num"
    fi
done

# Now maxNum contains the highest number found
if [ "$maxNum" -gt 0 ]; then
    for ((i=1; i<="$maxNum"; i++)); do
        echo ""
        echo "[ --- Run $i of $maxNum --- ]"
        echo ""

        # Relevant Input and Output folders
        inputPath="${baseFolder}${i}/Input"
        outputPath="${baseFolder}${i}/Output"

        # Check if the Input and Output directories exist
        if [ -d "$inputPath" ] && [ -d "$outputPath" ]; then
            # Call the process function/script with Input and Output paths as arguments
            echo "$processScript" -input "$inputPath" -output "$outputPath"
            "$processScript" -input "$inputPath" -output "$outputPath"

            # Copy Volumes.txt to the output folder for this run
            cp "./Output/Other/Volumes.txt" "$outputPath/Other/Volumes.txt"
        fi
    done
else
    echo "No suitable /run_N/ folders found."
fi
