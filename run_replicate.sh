#!/bin/bash

# Set the path to the SLiM executable
SLIM_EXEC="C:/msys64/mingw64/bin/slim.exe"

# Set the simulation script file
SIM_SCRIPT="C:\Users\User\Desktop\Master project\Simulation\Neutral.slim"

# Number of replicates to run
NUM_REPLICATES=100

# Output directory for result
OUTPUT_DIR=".\R_analysis\data\neutral"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Debug: Confirm output directory
echo "Output directory set to: $OUTPUT_DIR"

# Loop replicates
for (( i=1; i<=NUM_REPLICATES; i++ ))
do
	echo $i
	# Generate a unique seed for each replicate
	SEED=$(( RANDOM + i))
	
	REPLICATE_DIR="$OUTPUT_DIR\replicate_$i"
	mkdir -p "$REPLICATE_DIR"
	
	# Run the simulation with the seed and redirect Output
	echo "Running replicate $i with seed $SEED..."
	"$SLIM_EXEC" -s "$SEED" "$SIM_SCRIPT" > "$REPLICATE_DIR/log.txt"
	
	if [[ $? -ne 0 ]]; then
        echo "Replicate $i failed."
    else
        echo "Replicate $i completed successfully."
	
	# Move VCF files to the replicate directory
    mv *.vcf "$REPLICATE_DIR/" 2>/dev/null
    fi
done

echo "All replicates are done. Results are in the $OUTPUT_DIR directory."
	