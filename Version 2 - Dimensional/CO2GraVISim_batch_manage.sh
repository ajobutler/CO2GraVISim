#!/bin/bash

# CO2GraVISim_batch_manage (26/06/24)
# This script organises and runs several instances of CO2GraVISim_batch_run
# Built with the help of ChatGPT!

# TODO: Adapt this to identify which terminal command it should use, based on the OS?

# Set the base folder
baseFolder="./batch_runs/run_"

# Set the number of sets of runs that we want to perform at the same time
instances=4

# Initialize NumRuns to 0
NumRuns=0

# Scan directories to find the highest number N in /run_N/
for d in "$baseFolder"*; do
    if [ -d "$d" ]; then
        folderName=$(basename "$d")
        number="${folderName##*_}"
        
        # Convert to number for comparison
        num=$((number))
        if [ "$num" -gt "$NumRuns" ]; then
            NumRuns=$num
        fi
    fi
done

# Calculate the typical size of each chunk
ChunkSize=$((NumRuns / instances))
remainder=$((NumRuns % instances))

# Print a summary of the batch distribution process
echo
echo "-----------------------------------------------"
echo "Splitting $NumRuns runs into $instances groups of $ChunkSize, with remainder $remainder"
echo "-----------------------------------------------"
echo

# Initialize
start_index=1
stop_index=0

for i in $(seq 1 $instances); do
    stop_index=$((start_index + ChunkSize - 1))

    # If there is remainder still to distribute, add one to the stop index
    # and remove it from the remainder
    if [ $remainder -gt 0 ]; then
        stop_index=$((stop_index + 1))
        remainder=$((remainder - 1))
    fi

    # Make sure that we don't overshoot (should be covered by the above)
    # Uncomment if you want to handle overshoot cases
    if [ $i -eq $instances ]; then
        stop_index=$NumRuns
    fi

    
    current_dir=$(pwd)

    ## A Linux version (do I need to cd in here as well?)
    # gnome-terminal -- bash -c "./CO2GraVISim_batch_run.sh $start_index $stop_index; exec bash"

    ## Apple version
    osascript -e 'tell application "Terminal" to do script "cd '$current_dir' && bash CO2GraVISim_batch_run.sh '$start_index' '$stop_index'"'

    start_index=$((stop_index + 1))
done
