#!/bin/bash

# Set the path to the SLiM executable
SLIM_EXEC="C:/msys64/mingw64/bin/slim.exe"
SIM_SCRIPT="C:/Users/SwissHardware/Desktop/Master-Project/Simulation/Neutral.slim"
NUM_REPLICATES=100
OUTPUT_DIR="./R_analysis/data/neutral"


mkdir -p "$OUTPUT_DIR"


for (( i=1; i<=NUM_REPLICATES; i++ ))
do
    SEED=$(( RANDOM + i ))
    echo "Running replicate $i with seed $SEED..."
    "$SLIM_EXEC" -s "$SEED" "$SIM_SCRIPT" > "$OUTPUT_DIR/$i_log_$SEED.txt"

    if [[ $? -ne 0 ]]; then
        echo "Replicate $i failed."
    else
        mv *.vcf "$OUTPUT_DIR/" 2>/dev/null
    fi
done

echo "All replicates completed. Results stored in $OUTPUT_DIR."
