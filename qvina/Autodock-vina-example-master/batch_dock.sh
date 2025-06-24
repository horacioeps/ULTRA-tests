#!/bin/bash

# Batch docking script for QuickVina
# Usage: ./batch_dock.sh

echo "Starting batch docking with QuickVina..."
echo "=============================================="

# Create results directory
mkdir -p batch_results

# Initialize results file
echo "Ligand_Name,Best_Affinity_kcal_mol,Second_Best_kcal_mol,Third_Best_kcal_mol,Runtime_seconds" > batch_results/docking_results.csv

# Function to run docking for a single ligand
dock_ligand() {
    local ligand_file=$1
    local ligand_name=$(basename "$ligand_file" .pdbqt)
    
    echo "Docking $ligand_name..."
    
    # Start timing
    start_time=$(date +%s.%N)
    
    # Run QuickVina
    ../bin/vina --receptor pocket.pdbqt \
                --ligand "$ligand_file" \
                --center_x 139 --center_y 145 --center_z 171 \
                --size_x 25 --size_y 25 --size_z 25 \
                --num_modes 3 \
                --exhaustiveness 4 \
                --out "batch_results/${ligand_name}_out.pdbqt" \
                --log "batch_results/${ligand_name}_log.txt" \
                > "batch_results/${ligand_name}_output.txt" 2>&1
    
    # End timing
    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc -l)
    
    # Extract best three affinities from log
    if [ -f "batch_results/${ligand_name}_log.txt" ]; then
        affinities=$(grep "^\s*[1-3]\s" "batch_results/${ligand_name}_log.txt" | awk '{print $2}' | tr '\n' ',' | sed 's/,$//')
        if [ -z "$affinities" ]; then
            affinities="N/A,N/A,N/A"
        fi
    else
        affinities="ERROR,ERROR,ERROR"
    fi
    
    # Add to results CSV
    echo "$ligand_name,$affinities,$runtime" >> batch_results/docking_results.csv
    
    echo "  -> Completed $ligand_name (${runtime}s)"
}

# Export function for parallel execution
export -f dock_ligand

# Find all ligand files and run docking
ligand_files=(ligand_database/*.pdbqt)
if [ ${#ligand_files[@]} -eq 0 ]; then
    echo "No ligand files found in ligand_database/"
    echo "Creating sample ligands..."
    ./create_sample_ligands.sh
fi

# Run docking for each ligand
for ligand_file in ligand_database/*.pdbqt; do
    if [ -f "$ligand_file" ]; then
        dock_ligand "$ligand_file"
    fi
done

echo ""
echo "Batch docking completed!"
echo "Results saved to: batch_results/docking_results.csv"
echo ""
echo "Summary table:"
echo "=============="
column -t -s ',' batch_results/docking_results.csv 