#!/bin/bash

# Create sample ligands script
echo "Creating sample ligands for batch docking..."

# Create ligand database directory
mkdir -p ligand_database

# Common drug molecules with SMILES strings
declare -A molecules=(
    ["ibuprofen"]="CC(C)Cc1ccc(C(C)C(=O)O)cc1"
    ["aspirin"]="CC(=O)OC1=CC=CC=C1C(=O)O"
    ["acetaminophen"]="CC(=O)Nc1ccc(O)cc1"
    ["caffeine"]="CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    ["lisinopril"]="CCCCN1CCCC1C(=O)N2CCCC2C(=O)N3CCC(CC3)C(=O)O"
    ["metformin"]="CN(C)C(=N)NC(=N)N"
    ["atorvastatin"]="CC(C)c1c(C(=O)Nc2ccccc2F)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O"
    ["warfarin"]="CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"
    ["diazepam"]="CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"
    ["morphine"]="CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5"
)

# Check if obabel is available
if ! command -v obabel &> /dev/null; then
    echo "Open Babel (obabel) not found. Installing via brew..."
    brew install open-babel
fi

# Function to create ligand from SMILES
create_ligand() {
    local name=$1
    local smiles=$2
    
    echo "Creating $name..."
    
    # Convert SMILES to SDF with 3D coordinates
    echo "$smiles" | obabel -ismi -osdf --gen3d -O "ligand_database/${name}.sdf" 2>/dev/null
    
    # Convert SDF to PDBQT
    if [ -f "ligand_database/${name}.sdf" ]; then
        obabel "ligand_database/${name}.sdf" -opdbqt -O "ligand_database/${name}.pdbqt" \
               --partialcharge gasteiger 2>/dev/null
        # Clean up intermediate file
        rm -f "ligand_database/${name}.sdf"
        
        if [ -f "ligand_database/${name}.pdbqt" ]; then
            echo "  ✓ Created $name.pdbqt"
        else
            echo "  ✗ Failed to create $name.pdbqt"
        fi
    else
        echo "  ✗ Failed to create $name.sdf"
    fi
}

# Create all ligands
for name in "${!molecules[@]}"; do
    smiles="${molecules[$name]}"
    create_ligand "$name" "$smiles"
done

# Also copy the original ligand for comparison
cp ligand-b.pdbqt ligand_database/original_ligand.pdbqt

echo ""
echo "Created ligands:"
ls -la ligand_database/*.pdbqt | wc -l | xargs echo "Total ligands:"
ls ligand_database/*.pdbqt 