#!/bin/bash

# Set the path to the SLiM executable
SLIM_EXEC="C:/msys64/mingw64/bin/slim.exe"
SIM_SCRIPT="C:/Users/User/Desktop/Master project/Simulation/Neutral.slim"
NUM_REPLICATES=100
OUTPUT_DIR="./R_analysis/data/neutral_raw"
SAMPLE_DIR="./R_analysis/data/processed"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$SAMPLE_DIR"

for (( i=1; i<=NUM_REPLICATES; i++ ))
do
    REPLICATE_DIR="$OUTPUT_DIR/replicate_$i"
    mkdir -p "$REPLICATE_DIR"

    SEED=$(( RANDOM + i ))
    echo "Running replicate $i with seed $SEED..."
    "$SLIM_EXEC" -s "$SEED" "$SIM_SCRIPT" > "$REPLICATE_DIR/log.txt"

    if [[ $? -ne 0 ]]; then
        echo "Replicate $i failed."
    else
        mv *.vcf "$REPLICATE_DIR/" 2>/dev/null
		mv *.txt "$REPLICATE_DIR/" 2>/dev/null
        python preprocessing.py "$REPLICATE_DIR" "$SAMPLE_DIR"
        echo "Replicate $i preprocessing completed."
    fi
done

echo "All replicates completed. Results stored in $SAMPLE_DIR."
