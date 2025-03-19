#!/bin/bash

# CO2GraVISim_batch_run (22/04/24)
# This is performs a batch run of CO2GraVISim over the inputs stored in
# batch_runs/run_*/Input, with the resulting outputs stored in 
# batch_runs/run_*/Output

# Generated with ChatGPT 3.5's help!


# Function to get the current date and time
start_time=$(date +%s)
start_timestamp=$(date +"%Y-%m-%d %H:%M:%S")

# Print start timestamp
echo "Start time: $start_timestamp"



baseFolder="./batch_runs/run_"
processScript="./CO2GraVISim_single_run"


# Other flags to add when calling processScript
# For the unconfined version of the solver, add -unconfined
# runFlags=""
runFlags="-unconfined"



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

# Check for command-line inputs for the loop indices
if [ -z "$1" ]; then
    start_index=1
else
    start_index=$1
fi

if [ -z "$2" ]; then
    stop_index=$maxNum
else
    stop_index=$2
fi

# Make sure that start_index <= stop_index
if [ "$start_index" -gt "$stop_index" ]; then
    echo "Start index is greater than Stop index: [$start_index, $stop_index]"
    echo "Halting."
    exit 1
fi

# Make sure the endpoint is at most the number of run folders
if [ "$stop_index" -gt "$maxNum" ]; then
    stop_index=$maxNum
fi

# Number of runs to perform, and a counter for the loop
run_total=$((stop_index - start_index + 1))
run_index=1



# Write start and stop indices to standard out, for clarity
echo
echo
echo "--------------------------------------"
echo "Start Index:              $start_index"
echo "Stop Index:               $stop_index"
echo "Number of runs to do:     $run_total"
echo "Number of runs available: $maxNum"
echo "--------------------------------------"
echo
echo

# Main Loop over batch-run folders
if [ "$maxNum" -gt 0 ]; then
    for ((i = start_index; i <= stop_index; i++)); do
        echo ""
        echo "[ --- Run $run_index of $run_total --- ]"
        echo ""

        # Relevant Input and Output folders
        batchPath="${baseFolder}${i}"
        XMLPath="${batchPath}\Input\run_${i}.xml"

        # Check if the Input and Output directories exist
        if [ -d "$batchPath" ]; then
            # Call the process function/script with Input and Output paths as arguments
            echo "$processScript" -batch "$batchPath" -xml "$XMLPath" "$runFlags"
            "$processScript" -batch "$batchPath" -xml "$XMLPath" "$runFlags"
        fi

        run_index=$((run_index + 1))
    done
else
    echo "No suitable /run_N/ folders found."
fi



# Get end time
end_time=$(date +%s)
end_timestamp=$(date +"%Y-%m-%d %H:%M:%S")
echo "End time: $end_timestamp"

# Compute elapsed time in seconds
elapsed=$((end_time - start_time))

# Convert elapsed time to HH:MM:SS format
hours=$((elapsed / 3600))
minutes=$(((elapsed % 3600) / 60))
seconds=$((elapsed % 60))

# Print start and end timestamps and Elapsed time
echo ****************************************************************
echo Start time  : $start_timestamp%
echo End time    : $end_timestamp
echo Elapsed time: !${hours}h ${minutes}m ${seconds}s
echo ****************************************************************
