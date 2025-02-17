#!/bin/bash

# Set the path to the SLiM executable
SLIM_EXEC="C:/msys64/mingw64/bin/slim.exe"

# Simulation paths
SIM_SCRIPT05="./Simulation/Neutral.slim"
SIM_SCRIPT01="./Simulation/Neutral_001MR.slim"
SIM_SCRIPT1="./Simulation/Neutral_01MR.slim"

# Output directories
OUTPUT_DIR05="./R_analysis/data/neutral_005"
OUTPUT_DIR01="./R_analysis/data/neutral_001"
OUTPUT_DIR1="./R_analysis/data/neutral_01"

# Number of replicates
NUM_REPLICATES=2

# Create output directories
mkdir -p "$OUTPUT_DIR05"
mkdir -p "$OUTPUT_DIR01"
mkdir -p "$OUTPUT_DIR1"

# Function to run a single replicate
run_replicate() {
    local sim_script=$1
    local output_dir=$2
    local sim_label=$3
    local replicate_num=$4
    local seed=$5

    echo "Running replicate $replicate_num with seed $seed for $sim_label..."
    local output_file="$output_dir/${sim_label}_rep${replicate_num}_log_${seed}.txt"
    "$SLIM_EXEC" -s "$seed" "$sim_script" > "$output_file"

    if [[ $? -ne 0 ]]; then
        echo "Replicate $replicate_num for $sim_label failed."
    else
        # Rename and move VCF files
        for vcf_file in *.vcf; do
            mv "$vcf_file" "$output_dir/${sim_label}_rep${replicate_num}_seed${seed}.vcf" 2>/dev/null
        done
    fi
}

# Function to manage tasks with a maximum of 3 in parallel
manage_tasks() {
    local task_count=0
    for sim_script in "$SIM_SCRIPT05" "$SIM_SCRIPT01" "$SIM_SCRIPT1"; do
        local sim_label
        local output_dir

        # Assign output directory and label based on simulation script
        case $sim_script in
            "$SIM_SCRIPT05") sim_label="neutral005"; output_dir="$OUTPUT_DIR05";;
            "$SIM_SCRIPT01") sim_label="neutral001"; output_dir="$OUTPUT_DIR01";;
            "$SIM_SCRIPT1") sim_label="neutral01"; output_dir="$OUTPUT_DIR1";;
        esac

        # Run all replicates for the current simulation script
        for (( i=1; i<=NUM_REPLICATES; i++ )); do
            local seed=$(( RANDOM + i ))
            run_replicate "$sim_script" "$output_dir" "$sim_label" "$i" "$seed" &
            task_count=$(( task_count + 1 ))

            # Wait if the maximum number of tasks is reached
            if (( task_count >= 3 )); then
                wait -n  # Wait for at least one task to finish
                task_count=$(( task_count - 1 ))  # Decrement task count
            fi
        done
    done

    # Wait for all remaining tasks to finish
    wait
}

# Run the simulations
manage_tasks

echo "All replicates completed. Results stored in respective directories."
