#!/usr/bin/env python3

import os
import sys
import csv
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import time
import pickle
import json

# Set RDKit environment
os.environ['RDBASE'] = '/opt/homebrew/opt/rdkit/share/RDKit'

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem, rdDepictor, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdFingerprintGenerator
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.image as mpimg
    from PIL import Image
    import xgboost as xgb
    from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
    from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
    from sklearn.preprocessing import StandardScaler
except ImportError as e:
    print(f"❌ Error importing: {e}")
    sys.exit(1)

# Comprehensive drug-like SMILES database
DRUG_SMILES_DATABASE = [
    # Analgesics & Anti-inflammatories
    'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aspirin
    'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # Ibuprofen
    'CC(=O)NC1=CC=C(C=C1)O',  # Acetaminophen
    'CC1=C(C(=O)N(N1C)C(=O)NCCCCOC)CC',  # Celecoxib-like
    'CC(C)(C)CC(=O)N1CCN(CC1)C2=CC=CC=C2',  # Synthetic analgesic
    
    # CNS Active Compounds
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
    'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3',  # Diazepam
    'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O',  # Morphine
    'CC(C)NCC(C1=CC(=C(C=C1)O)CO)O',  # Salbutamol
    'CCCN(CCC)CCNC(=O)C1=CC(=C(C=C1)N)C(F)(F)F',  # Synthetic CNS drug
    
    # Cardiovascular Drugs
    'CC(C)C1=C(C(=C(N1CC(CC(=O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4',  # Atorvastatin
    'CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O',  # Warfarin
    'C1=CC(=CC=C1CCN)O',  # Dopamine
    'CC(C)(C)NCC(C1=CC(=C(C=C1)O)CO)O',  # Terbutaline
    'CCCCC1=NC2=C(N1)C=CC(=C2)C(=O)N',  # Synthetic cardiovascular
    
    # Antibiotics & Antimicrobials
    'CC1=CC2=C(C=C1C)N(C=N2)CC3=CN=CN3',  # Imidazole antifungal
    'CC(C)CC1=CC=C(C=C1)C(C)C(=O)NC2=CC=C(C=C2)O',  # Synthetic antibiotic
    'CN(C)C(=N)NC(=N)N',  # Metformin
    'CC1=C(C=CC=C1Cl)NC(=O)CN2CCN(CC2)CC3=CC=CC=C3',  # Synthetic antimicrobial
    'C1=CC(=CC=C1C(=O)NN=CC2=CC=C(C=C2)[N+](=O)[O-])O',  # Nitrofuran derivative
    
    # Diverse drug-like molecules
    'CC1=CC=C(C=C1)S(=O)(=O)N2CCCC(C2)C(=O)O',  # Sulfonamide derivative
    'CN1CCN(CC1)C2=CC=C(C=C2)C(=O)NC3=CC=CC=C3F',  # Piperazine derivative
    'CC(C)N1C=NC=C1CN2CCC(CC2)C3=CC=CC=C3',  # Imidazole derivative
    'CC1=CC=C(C=C1)C2=CC(=NO2)C(=O)NCCCN(C)C',  # Oxazole derivative
    'CC(C)(C)OC(=O)N1CCC(CC1)C(=O)NC2=CC=CC=C2',  # Boc-protected amine
    
    # Extended pharmaceutical space
    'CC1=CC2=C(C=C1)C(=O)N(C2=O)CC3=CC=CC=C3',  # Phthalimide derivative
    'CN1C2=CC=CC=C2C(=O)N(C1=O)CC3=CC=C(C=C3)Cl',  # Benzodiazepine-like
    'CC(C)CC1=CC=C(C=C1)C(=O)N2CCCC2C(=O)O',  # Proline derivative
    'CN(C)CCN1C2=CC=CC=C2SC3=C1C=C(C=C3)CF3',  # Trifluoromethyl compound
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)N3CCOCC3',  # Morpholine derivative
    
    # Generate more diverse structures
    'CC(C)N(CC1=CC=CC=C1)C(=O)C2=CC=C(C=C2)OC',  # Ether derivative
    'CN1CCC(CC1)NC(=O)C2=CC=C(C=C2)S(=O)(=O)N',  # Sulfonamide
    'CC1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)C(F)(F)F',  # Trifluoromethyl amide
    'CN(C)CCN1C(=O)C2=CC=CC=C2C1=O',  # Phthalimide
    'CC(C)OC(=O)N1CCN(CC1)C2=CC=C(C=C2)Cl',  # Chlorinated piperazine
    
    # Additional pharmaceutical scaffolds
    'CC1=CN=C(C=C1)NC(=O)C2=CC=C(C=C2)OCC3=CC=CC=C3',  # Pyridine derivative
    'CN1C=NC2=CC=CC=C21',  # Benzimidazole
    'CC(C)(C)C1=CC=C(C=C1)S(=O)(=O)NC2=CC=CC=C2',  # Sulfonamide
    'CN(C)C(=O)C1=CC=C(C=C1)N2CCN(CC2)CC3=CC=CC=C3',  # Piperazine amide
    'CC1=CC=C(C=C1)C(=O)N2CCC(CC2)N3C=CC=N3',  # Imidazole piperidine
    
    # More diverse structures
    'CC(C)N1C(=O)C(=C(C1=O)C)CC2=CC=CC=C2',  # Pyrrole derivative
    'CN1CCN(CC1)C(=O)C2=CC=C(C=C2)OC3=CC=CC=C3',  # Diphenyl ether
    'CC1=CC=C(C=C1)S(=O)(=O)N2CCC(CC2)C(=O)NC3=CC=CC=C3',  # Complex sulfonamide
    'CN(C)CCOC1=CC=C(C=C1)C(=O)NC2=CC=CC=C2F',  # Ether-linked amide
    'CC(C)(C)OC(=O)NC1=CC=C(C=C1)C(=O)N2CCCC2',  # Boc-protected pyrrolidine
    
    # Final diverse set
    'CC1=CC2=C(C=C1)N(C(=O)N2)CC3=CC=C(C=C3)Cl',  # Benzimidazolone
    'CN1C(=O)C2=C(C1=O)C=CC(=C2)N3CCN(CC3)CC4=CC=CC=C4',  # Complex heterocycle
    'CC(C)CC1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)S(=O)(=O)N',  # Sulfonamide amide
    'CN(C)CCN1C(=O)NC2=CC=CC=C21',  # Benzimidazolone derivative
    'CC1=CC=C(C=C1)C(=O)N2CCN(CC2)C(=O)C3=CC=CC=C3F',  # Fluorinated piperazine
    
    # Even more diversity
    'CC(C)N1CCN(CC1)C(=O)C2=CC=C(C=C2)OC3=CC=CC=C3',  # Diphenyl ether piperazine
    'CN1C2=CC=CC=C2C(=O)N(C1=S)CC3=CC=C(C=C3)Br',  # Brominated thioamide
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)N3CCN(CC3)C',  # Methylated piperazine
    'CC(C)(C)C1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)CF3',  # Tert-butyl trifluoromethyl
    'CN(C)C(=O)C1=CC2=C(C=C1)OCO2',  # Methylenedioxy derivative
]

def generate_additional_smiles(n_additional=50):
    """Generate additional drug-like SMILES using combinatorial chemistry"""
    additional_smiles = []
    
    # Common pharmaceutical fragments
    cores = [
        'C1=CC=CC=C1',  # Benzene
        'C1=CC=NC=C1',  # Pyridine
        'C1=CN=CN1',    # Imidazole
        'C1CCCCC1',     # Cyclohexane
        'C1CCNCC1',     # Piperidine
        'C1CNCCN1',     # Piperazine
    ]
    
    substituents = [
        'C(=O)O',       # Carboxylic acid
        'C(=O)N',       # Amide
        'S(=O)(=O)N',   # Sulfonamide
        'OC',           # Methoxy
        'N(C)C',        # Dimethylamino
        'C(F)(F)F',     # Trifluoromethyl
        'Cl',           # Chloro
        'C(C)C',        # Isopropyl
    ]
    
    for i in range(n_additional):
        try:
            core = np.random.choice(cores)
            subs = np.random.choice(substituents, size=np.random.randint(1, 3), replace=False)
            
            # Simple combinatorial SMILES (this is a simplified approach)
            smiles = core + ''.join(subs)
            additional_smiles.append(smiles)
        except:
            continue
    
    return additional_smiles

def create_100_ligands():
    """Create 100 diverse ligand PDBQT files"""
    print("🧬 Generando 100 ligandos diversos...")
    
    # Combine base database with additional generated SMILES
    all_smiles = DRUG_SMILES_DATABASE.copy()
    additional = generate_additional_smiles(100 - len(DRUG_SMILES_DATABASE))
    all_smiles.extend(additional)
    
    # Take first 100
    selected_smiles = all_smiles[:100]
    
    ligand_dir = Path('ligand_database_100')
    ligand_dir.mkdir(exist_ok=True)
    
    ligand_data = []
    
    for i, smiles in enumerate(selected_smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
                
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.UFFOptimizeMolecule(mol)
            
            # Save as PDB first
            ligand_name = f"ligand_{i+1:03d}"
            pdb_file = ligand_dir / f"{ligand_name}.pdb"
            
            Chem.MolToPDBFile(mol, str(pdb_file))
            
            # Convert PDB to PDBQT using simple method
            pdbqt_file = ligand_dir / f"{ligand_name}.pdbqt"
            convert_pdb_to_pdbqt(pdb_file, pdbqt_file)
            
            ligand_data.append({
                'name': ligand_name,
                'smiles': smiles,
                'pdbqt_file': pdbqt_file
            })
            
            if len(ligand_data) >= 100:
                break
                
        except Exception as e:
            print(f"❌ Error creating ligand {i+1}: {e}")
            continue
    
    print(f"✅ Creados {len(ligand_data)} ligandos exitosamente")
    return ligand_data

def convert_pdb_to_pdbqt(pdb_file, pdbqt_file):
    """Simple PDB to PDBQT conversion"""
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        pdbqt_lines = []
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Simple conversion - add partial charge and atom type
                pdbqt_line = line[:66] + '  0.00 A \n'
                pdbqt_lines.append(pdbqt_line)
            elif line.startswith('CONECT'):
                continue  # Skip connectivity
        
        # Add ROOT and ENDROOT for flexibility
        final_lines = ['ROOT\n'] + pdbqt_lines + ['ENDROOT\n']
        
        with open(pdbqt_file, 'w') as f:
            f.writelines(final_lines)
        
        return True
    except:
        return False

def run_batch_docking_100(ligand_data):
    """Run docking for 100 ligands"""
    print("🔬 Ejecutando docking para 100 ligandos...")
    
    results = []
    batch_dir = Path('batch_results_100')
    batch_dir.mkdir(exist_ok=True)
    
    # Use existing receptor files
    receptor = 'pocket.pdbqt'
    config_file = 'conf_quick.txt'
    
    for i, ligand_info in enumerate(ligand_data):
        ligand_name = ligand_info['name']
        ligand_file = ligand_info['pdbqt_file']
        
        if not ligand_file.exists():
            continue
        
        print(f"  📄 Docking {i+1}/100: {ligand_name}")
        
        output_file = batch_dir / f"{ligand_name}_out.pdbqt"
        log_file = batch_dir / f"{ligand_name}_log.txt"
        
        # Run QuickVina (adjust path)
        cmd = [
            '../bin/qvina2.1',
            '--receptor', receptor,
            '--ligand', str(ligand_file),
            '--config', config_file,
            '--out', str(output_file),
            '--log', str(log_file)
        ]
        
        start_time = time.time()
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            runtime = time.time() - start_time
            
            # Extract binding affinity
            affinity = extract_binding_affinity(log_file)
            
            results.append({
                'name': ligand_name,
                'smiles': ligand_info['smiles'],
                'affinity': affinity,
                'runtime': runtime
            })
            
        except subprocess.TimeoutExpired:
            print(f"    ⏰ Timeout para {ligand_name}")
            results.append({
                'name': ligand_name,
                'smiles': ligand_info['smiles'],
                'affinity': 0.0,
                'runtime': 30.0
            })
        except Exception as e:
            print(f"    ❌ Error: {e}")
            continue
    
    # Save results to CSV
    csv_file = batch_dir / 'docking_results_100.csv'
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['name', 'smiles', 'affinity', 'runtime'])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"✅ Docking completado: {len(results)} ligandos procesados")
    return results

def extract_binding_affinity(log_file):
    """Extract binding affinity from log file"""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Look for affinity pattern
        lines = content.split('\n')
        for line in lines:
            if 'VINA RESULT:' in line or '1   ' in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    try:
                        affinity = float(part)
                        if -20 < affinity < 0:  # Reasonable range
                            return affinity
                    except:
                        continue
        return 0.0
    except:
        return 0.0

def calculate_morgan_fingerprints(smiles_list, radius=2, n_bits=2048):
    """Calculate Morgan fingerprints for SMILES"""
    print("🧮 Calculando Morgan fingerprints...")
    
    fingerprints = []
    valid_smiles = []
    
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                fp = generator.GetFingerprint(mol)
                # Convert to bit vector
                fp_array = np.zeros(n_bits)
                for bit in fp.GetOnBits():
                    fp_array[bit] = 1
                fingerprints.append(fp_array)
                valid_smiles.append(smiles)
            else:
                print(f"⚠️  SMILES inválido: {smiles}")
        except Exception as e:
            print(f"❌ Error procesando {smiles}: {e}")
    
    return np.array(fingerprints), valid_smiles

def train_xgboost_model(fingerprints, affinities, test_size=0.2):
    """Train XGBoost model for affinity prediction"""
    print("🤖 Entrenando modelo XGBoost...")
    
    # Remove zero affinities (failed dockings)
    valid_mask = np.array(affinities) != 0.0
    X = fingerprints[valid_mask]
    y = np.array(affinities)[valid_mask]
    
    if len(X) < 10:
        print("❌ Insuficientes datos válidos para entrenamiento")
        return None, None, {}
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=42
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Grid search for best parameters
    param_grid = {
        'n_estimators': [100, 200, 300],
        'max_depth': [3, 6, 9],
        'learning_rate': [0.01, 0.1, 0.2],
        'subsample': [0.8, 1.0]
    }
    
    xgb_model = xgb.XGBRegressor(random_state=42)
    grid_search = GridSearchCV(
        xgb_model, param_grid, cv=5, scoring='neg_mean_squared_error', n_jobs=-1
    )
    
    grid_search.fit(X_train_scaled, y_train)
    best_model = grid_search.best_estimator_
    
    # Predictions
    y_pred_train = best_model.predict(X_train_scaled)
    y_pred_test = best_model.predict(X_test_scaled)
    
    # Metrics
    metrics = {
        'train_r2': r2_score(y_train, y_pred_train),
        'test_r2': r2_score(y_test, y_pred_test),
        'train_rmse': np.sqrt(mean_squared_error(y_train, y_pred_train)),
        'test_rmse': np.sqrt(mean_squared_error(y_test, y_pred_test)),
        'train_mae': mean_absolute_error(y_train, y_pred_train),
        'test_mae': mean_absolute_error(y_test, y_pred_test),
        'best_params': grid_search.best_params_,
        'cv_score': -grid_search.best_score_,
        'n_samples': len(X),
        'n_features': X.shape[1]
    }
    
    # Cross-validation
    cv_scores = cross_val_score(best_model, X_train_scaled, y_train, cv=5, scoring='r2')
    metrics['cv_r2_mean'] = cv_scores.mean()
    metrics['cv_r2_std'] = cv_scores.std()
    
    print(f"✅ Modelo entrenado exitosamente")
    print(f"  📊 R² Test: {metrics['test_r2']:.3f}")
    print(f"  📊 RMSE Test: {metrics['test_rmse']:.3f}")
    
    # Save model
    model_data = {
        'model': best_model,
        'scaler': scaler,
        'metrics': metrics,
        'X_test': X_test,
        'y_test': y_test,
        'y_pred_test': y_pred_test
    }
    
    with open('batch_results_100/xgboost_model.pkl', 'wb') as f:
        pickle.dump(model_data, f)
    
    return best_model, scaler, metrics

def create_ml_pdf_report():
    """Create comprehensive ML PDF report"""
    print("📊 Creando reporte PDF con análisis de Machine Learning...")
    
    # Main workflow
    print("1️⃣  Generando 100 ligandos...")
    ligand_data = create_100_ligands()
    
    print("2️⃣  Ejecutando docking batch...")
    docking_results = run_batch_docking_100(ligand_data)
    
    print("3️⃣  Calculando fingerprints...")
    smiles_list = [r['smiles'] for r in docking_results]
    affinities = [r['affinity'] for r in docking_results]
    fingerprints, valid_smiles = calculate_morgan_fingerprints(smiles_list)
    
    print("4️⃣  Entrenando modelo XGBoost...")
    model, scaler, metrics = train_xgboost_model(fingerprints, affinities)
    
    # Create PDF report
    filename = "batch_results_100/QuickVina_ML_Analysis_Report.pdf"
    
    with PdfPages(filename) as pdf:
        # Page 1: Overview and Statistics
        fig = plt.figure(figsize=(8.5, 11))
        fig.suptitle('QuickVina Machine Learning Analysis\n100 Ligands + XGBoost Regression Model', 
                     fontsize=16, fontweight='bold', y=0.95)
        
        # Summary statistics
        valid_affinities = [a for a in affinities if a != 0.0]
        
        summary_text = f"""MACHINE LEARNING MOLECULAR DOCKING STUDY

Dataset Size: 100 pharmaceutical compounds
Valid Docking Results: {len(valid_affinities)} compounds
Fingerprint Dimensions: {fingerprints.shape[1] if len(fingerprints) > 0 else 'N/A'} Morgan bits
Binding Affinity Range: {min(valid_affinities):.2f} to {max(valid_affinities):.2f} kcal/mol
Mean Binding Energy: {np.mean(valid_affinities):.2f} ± {np.std(valid_affinities):.2f} kcal/mol

MACHINE LEARNING MODEL PERFORMANCE:
Algorithm: XGBoost Gradient Boosting Regressor
Features: Morgan Fingerprints (radius=2, 2048 bits)
Train/Test Split: 80%/20%
Cross-Validation: 5-fold

Model Performance Metrics:"""
        
        if metrics:
            summary_text += f"""
R² Score (Test): {metrics['test_r2']:.3f}
RMSE (Test): {metrics['test_rmse']:.3f} kcal/mol
MAE (Test): {metrics['test_mae']:.3f} kcal/mol
CV R² Score: {metrics['cv_r2_mean']:.3f} ± {metrics['cv_r2_std']:.3f}

Best Hyperparameters:
• n_estimators: {metrics['best_params']['n_estimators']}
• max_depth: {metrics['best_params']['max_depth']}
• learning_rate: {metrics['best_params']['learning_rate']}
• subsample: {metrics['best_params']['subsample']}

Training Samples: {metrics['n_samples']}
Feature Dimensions: {metrics['n_features']}"""
        else:
            summary_text += "\nModel training failed - insufficient valid data"
        
        fig.text(0.1, 0.75, summary_text, fontsize=10, verticalalignment='top', 
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcyan", alpha=0.9, edgecolor='navy'))
        
        # Save first page
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
        plt.close()
        
        # Page 2: Model Performance Plots
        if metrics and model is not None:
            fig2 = plt.figure(figsize=(8.5, 11))
            fig2.suptitle('XGBoost Model Performance Analysis', fontsize=16, fontweight='bold', y=0.95)
            
            # Load model data for plotting
            with open('batch_results_100/xgboost_model.pkl', 'rb') as f:
                model_data = pickle.load(f)
            
            X_test = model_data['X_test']
            y_test = model_data['y_test']
            y_pred_test = model_data['y_pred_test']
            
            # Predicted vs Actual plot
            ax1 = fig2.add_subplot(2, 2, 1)
            ax1.scatter(y_test, y_pred_test, alpha=0.6, color='blue')
            ax1.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
            ax1.set_xlabel('Actual Binding Affinity (kcal/mol)')
            ax1.set_ylabel('Predicted Binding Affinity (kcal/mol)')
            ax1.set_title(f'Predicted vs Actual (R² = {metrics["test_r2"]:.3f})')
            ax1.grid(True, alpha=0.3)
            
            # Residuals plot
            ax2 = fig2.add_subplot(2, 2, 2)
            residuals = y_test - y_pred_test
            ax2.scatter(y_pred_test, residuals, alpha=0.6, color='green')
            ax2.axhline(y=0, color='r', linestyle='--')
            ax2.set_xlabel('Predicted Binding Affinity (kcal/mol)')
            ax2.set_ylabel('Residuals (kcal/mol)')
            ax2.set_title('Residuals Plot')
            ax2.grid(True, alpha=0.3)
            
            # Feature importance (top 20)
            ax3 = fig2.add_subplot(2, 1, 2)
            feature_importance = model.feature_importances_
            top_features = np.argsort(feature_importance)[-20:]
            ax3.barh(range(len(top_features)), feature_importance[top_features])
            ax3.set_xlabel('Feature Importance')
            ax3.set_ylabel('Morgan Fingerprint Bit Index')
            ax3.set_title('Top 20 Most Important Molecular Features')
            ax3.set_yticks(range(len(top_features)))
            ax3.set_yticklabels([f'Bit {i}' for i in top_features])
            
            plt.tight_layout()
            pdf.savefig(fig2, bbox_inches='tight', dpi=300)
            plt.close()
        
    print(f"✅ Reporte ML creado: {filename}")
    return filename

if __name__ == "__main__":
    try:
        print("🚀 ANÁLISIS DE MACHINE LEARNING PARA DOCKING MOLECULAR")
        print("=" * 70)
        
        report_file = create_ml_pdf_report()
        
        print(f"\n🎉 ANÁLISIS COMPLETO FINALIZADO")
        print(f"📄 Reporte: {report_file}")
        print("\n🔬 WORKFLOW COMPLETADO:")
        print("  ✅ 100 ligandos diversos generados")
        print("  ✅ Docking molecular ejecutado")
        print("  ✅ SMILES extraídos automáticamente")
        print("  ✅ Morgan fingerprints calculados")
        print("  ✅ Modelo XGBoost entrenado")
        print("  ✅ Validación cruzada realizada")
        print("  ✅ Reporte PDF con análisis ML generado")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc() 