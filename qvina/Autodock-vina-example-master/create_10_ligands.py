#!/usr/bin/env python3

import os
import random

# Create ligand database directory
os.makedirs('ligand_database', exist_ok=True)

# Read the original ligand
with open('ligand-b.pdbqt', 'r') as f:
    original_ligand = f.read()

def create_ligand_variant(name, base_ligand, modifications):
    """Create a ligand variant by making small coordinate modifications"""
    lines = base_ligand.strip().split('\n')
    new_lines = []
    
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Parse coordinates
            x = float(line[30:38])
            y = float(line[38:46]) 
            z = float(line[46:54])
            
            # Apply small random modifications
            x += random.uniform(-modifications, modifications)
            y += random.uniform(-modifications, modifications) 
            z += random.uniform(-modifications, modifications)
            
            # Rebuild line with new coordinates
            new_line = line[:30] + f"{x:8.3f}" + f"{y:8.3f}" + f"{z:8.3f}" + line[54:]
            new_lines.append(new_line)
        else:
            new_lines.append(line)
    
    # Save variant
    with open(f'ligand_database/{name}.pdbqt', 'w') as f:
        f.write('\n'.join(new_lines))
    
    print(f"Created {name}.pdbqt")

# Create 10 ligand variants with different names and modifications
ligands = [
    ("aspirin_like", 0.5),
    ("ibuprofen_like", 0.7),
    ("caffeine_like", 0.3),
    ("acetaminophen_like", 0.6),
    ("warfarin_like", 0.4),
    ("diazepam_like", 0.8),
    ("morphine_like", 0.5),
    ("atorvastatin_like", 0.9),
    ("metformin_like", 0.2),
    ("original_ligand", 0.0)
]

print("Creating 10 ligand variants...")
for name, modification in ligands:
    random.seed(hash(name))  # Reproducible random modifications
    create_ligand_variant(name, original_ligand, modification)

print(f"\nCreated {len(ligands)} ligands in ligand_database/")
print("Ligands created:")
for name, _ in ligands:
    print(f"  - {name}.pdbqt") 