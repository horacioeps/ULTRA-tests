#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Set RDKit environment
os.environ['RDBASE'] = '/opt/homebrew/opt/rdkit/share/RDKit'

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    from rdkit.Chem import rdFingerprintGenerator
    import xgboost as xgb
    from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
    from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
    from sklearn.preprocessing import StandardScaler
except ImportError as e:
    print(f"❌ Error importing: {e}")
    sys.exit(1)

# Comprehensive 100 drug-like SMILES database
DRUG_SMILES_100 = [
    # Analgesics & Anti-inflammatories (20)
    'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aspirin
    'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # Ibuprofen
    'CC(=O)NC1=CC=C(C=C1)O',  # Acetaminophen
    'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O',  # Morphine
    'CC1=CC=C(C=C1)C(=O)C2=CC=C(C=C2)O',  # Basic NSAID
    'CC(C)CC1=CC=C(C=C1)C(C)C(=O)NC2=CC=CC=C2',  # Ibuprofen derivative
    'CC(=O)NC1=CC=C(C=C1)OC',  # Acetaminophen analog
    'CC1=CC=C(C=C1)S(=O)(=O)NC2=CC=CC=C2',  # Sulfonamide analgesic
    'CC(C)(C)CC(=O)N1CCN(CC1)C2=CC=CC=C2',  # Synthetic analgesic
    'CN(C)CCN1C(=O)C2=CC=CC=C2C1=O',  # Phthalimide analgesic
    'CC1=CC=C(C=C1)C(=O)NCC2=CC=C(C=C2)O',  # Hydroxylated amide
    'CC(C)OC(=O)N1CCC(CC1)C(=O)O',  # Proline derivative
    'CN1C(=O)C(=C(C1=O)C)CC2=CC=CC=C2',  # Pyrrole derivative
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)Cl',  # Chlorinated amide
    'CC(C)N1CCN(CC1)C(=O)C2=CC=CC=C2',  # Piperazine amide
    'CN(C)C(=O)C1=CC=C(C=C1)OC2=CC=CC=C2',  # Diphenyl ether
    'CC1=CC=C(C=C1)S(=O)(=O)N2CCC(CC2)C(=O)O',  # Sulfonamide derivative
    'CN1CCN(CC1)C(=O)C2=CC=C(C=C2)F',  # Fluorinated piperazine
    'CC(C)(C)C1=CC=C(C=C1)C(=O)NC2=CC=CC=C2',  # Tert-butyl amide
    'CC1=CC=C(C=C1)C(=O)N2CCN(CC2)CC3=CC=CC=C3',  # Piperazine derivative
    
    # CNS Active Compounds (20)
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
    'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3',  # Diazepam
    'CC(C)NCC(C1=CC(=C(C=C1)O)CO)O',  # Salbutamol
    'CN(C)CCN1C2=CC=CC=C2SC3=CC=CC=C31',  # Tricyclic
    'CC1=CN=C(C=C1)NC(=O)C2=CC=CC=C2',  # Pyridine derivative
    'CN1C2=CC=CC=C2C(=O)N(C1=O)CC3=CC=CC=C3',  # Benzodiazepine-like
    'CC(C)N(CC1=CC=CC=C1)C(=O)C2=CC=CC=C2',  # Dopaminergic
    'CN(C)CCOC1=CC=C(C=C1)C(=O)NC2=CC=CC=C2',  # Ether-linked amide
    'CC1=CC2=C(C=C1)N(C(=O)N2)CC3=CC=CC=C3',  # Benzimidazolone
    'CN1CCN(CC1)C2=CC=C(C=C2)C(=O)NC3=CC=CC=C3',  # Piperazine CNS
    'CC(C)(C)OC(=O)N1CCC(CC1)C(=O)NC2=CC=CC=C2',  # Boc-protected
    'CN1C(=O)C2=C(C1=O)C=CC(=C2)N3CCN(CC3)C',  # Complex heterocycle
    'CC1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)N3CCOCC3',  # Morpholine
    'CN(C)CCN1C(=O)NC2=CC=CC=C21',  # Benzimidazolone
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)CF3',  # Trifluoromethyl
    'CN1CCC(CC1)NC(=O)C2=CC=C(C=C2)OC',  # Piperidine ether
    'CC(C)CC1=CC=C(C=C1)C(=O)N2CCCC2',  # Pyrrolidine amide
    'CN(C)C1=CC=C(C=C1)C(=O)NC2=CC=CC=C2F',  # Fluorinated aniline
    'CC1=CC=C(C=C1)S(=O)(=O)NC2=CC=C(C=C2)C',  # Sulfonamide CNS
    'CN1CCN(CC1)C(=O)C2=CC=C(C=C2)OCC3=CC=CC=C3',  # Benzyl ether,
    
    # Cardiovascular Drugs (20)
    'CC(C)C1=C(C(=C(N1CC(CC(=O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4',  # Atorvastatin
    'CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O',  # Warfarin
    'CN(C)C(=N)NC(=N)N',  # Metformin
    'C1=CC(=CC=C1CCN)O',  # Dopamine
    'CC(C)(C)NCC(C1=CC(=C(C=C1)O)CO)O',  # Terbutaline
    'CC1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)S(=O)(=O)N',  # CV sulfonamide
    'CN1CCN(CC1)C2=CC=C(C=C2)C(=O)NC3=CC=CC=C3',  # CV piperazine
    'CC(C)OC(=O)NC1=CC=C(C=C1)C(=O)O',  # Carboxylic acid derivative
    'CC1=CC=C(C=C1)OC(=O)C2=CC=C(C=C2)N',  # Aminobenzoate
    'CN(C)CCN1C(=O)C2=CC=CC=C2C1=S',  # Thioamide CV
    'CC(C)CC1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)OH',  # Hydroxylated CV
    'CC1=CC=C(C=C1)S(=O)(=O)N2CCN(CC2)C3=CC=CC=C3',  # Sulfonamide piperazine
    'CN1C2=CC=CC=C2C(=O)N(C1=O)C3=CC=CC=C3',  # Benzodiazepine CV
    'CC(C)(C)C1=CC=C(C=C1)C(=O)N2CCC(CC2)O',  # Hydroxylated piperidine
    'CN(C)C(=O)C1=CC=C(C=C1)N2CCN(CC2)C',  # Methylated piperazine
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)OC',  # Methoxy CV drug
    'CC(C)N1CCN(CC1)C(=O)C2=CC=CC=C2F',  # Fluorinated CV
    'CN1CCC(CC1)C(=O)NC2=CC=C(C=C2)Cl',  # Chlorinated CV
    'CC1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)S(=O)(=O)C',  # Methylsulfonyl
    'CN(C)CCO1C=CC=C1C(=O)NC2=CC=CC=C2',  # Furan derivative
    
    # Antibiotics & Antimicrobials (20)
    'CC1=CC2=C(C=C1C)N(C=N2)CC3=CN=CN3',  # Imidazole antifungal
    'CC(C)CC1=CC=C(C=C1)C(C)C(=O)NC2=CC=C(C=C2)O',  # Antibiotic
    'CC1=C(C=CC=C1Cl)NC(=O)CN2CCN(CC2)CC3=CC=CC=C3',  # Antimicrobial
    'C1=CC(=CC=C1C(=O)NN=CC2=CC=C(C=C2)[N+](=O)[O-])O',  # Nitrofuran
    'CN1CCN(CC1)C2=CC=C(C=C2)C(=O)NC3=CC=C(C=C3)N',  # Amino antibiotic
    'CC1=CC=C(C=C1)S(=O)(=O)NC2=CC=C(C=C2)C(=O)O',  # Sulfa drug
    'CN(C)C1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)NO2',  # Nitro antibiotic
    'CC(C)OC(=O)N1CCN(CC1)C2=CC=C(C=C2)Br',  # Brominated
    'CN1C2=CC=CC=C2C(=O)N(C1=O)CC3=CC=C(C=C3)I',  # Iodinated
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)S(=O)(=O)NH2',  # Sulfonamide antibiotic
    'CN(C)CCN1C(=O)C2=CC=C(C=C2)NC1=O',  # Cyclic antibiotic
    'CC(C)CC1=CC=C(C=C1)C(=O)N2CCC(CC2)N',  # Amino piperidine
    'CN1CCN(CC1)C(=O)C2=CC=C(C=C2)N(O)O',  # Hydroxylamine
    'CC1=CC=C(C=C1)OC(=O)C2=CC=C(C=C2)Cl',  # Chlorinated ester
    'CN(C)C(=O)C1=CC=C(C=C1)NC2=CC=CC=C2Br',  # Brominated amide
    'CC(C)N1CCN(CC1)C(=O)C2=CC=C(C=C2)F',  # Fluorinated antibiotic
    'CN1C(=O)C2=CC=C(C=C2)NC1=S',  # Thioamide antibiotic
    'CC1=CC=C(C=C1)S(=O)(=O)N2CCC(CC2)NH2',  # Amino sulfonamide
    'CN(C)CCO1C=CC=C1C(=O)N2CCN(CC2)C',  # Furan antibiotic
    'CC(C)(C)OC(=O)N1CCC(CC1)C(=O)NC2=CC=C(C=C2)Cl',  # Protected antibiotic
    
    # Diverse Pharmaceuticals (20)
    'CC1=CC2=C(C=C1)C(=O)N(C2=O)CC3=CC=CC=C3',  # Phthalimide
    'CN1C2=CC=CC=C2C(=O)N(C1=O)CC3=CC=C(C=C3)Cl',  # Chlorinated heterocycle
    'CC(C)CC1=CC=C(C=C1)C(=O)N2CCCC2C(=O)O',  # Proline complex
    'CN(C)CCN1C2=CC=CC=C2SC3=C1C=C(C=C3)CF3',  # Trifluoromethyl complex
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)N3CCOCC3',  # Morpholine complex
    'CN(C)C(=O)C1=CC=C(C=C1)N2CCN(CC2)CC3=CC=CC=C3',  # Complex piperazine
    'CC(C)(C)C1=CC=C(C=C1)S(=O)(=O)NC2=CC=CC=C2',  # Tert-butyl sulfonamide
    'CN1CCN(CC1)C(=O)C2=CC=C(C=C2)OC3=CC=CC=C3',  # Diphenyl ether complex
    'CC1=CC=C(C=C1)C(=O)N2CCC(CC2)N3C=CC=N3',  # Imidazole piperidine
    'CN(C)CCOC1=CC=C(C=C1)C(=O)NC2=CC=CC=C2F',  # Fluorinated ether
    'CC(C)(C)OC(=O)NC1=CC=C(C=C1)C(=O)N2CCCC2',  # Boc pyrrolidine
    'CC1=CC2=C(C=C1)N(C(=O)N2)CC3=CC=C(C=C3)Cl',  # Chlorinated benzimidazole
    'CN1C(=O)C2=C(C1=O)C=CC(=C2)N3CCN(CC3)CC4=CC=CC=C4',  # Ultra-complex
    'CC(C)CC1=CC=C(C=C1)C(=O)NC2=CC=C(C=C2)S(=O)(=O)N',  # Complex sulfonamide
    'CN(C)CCN1C(=O)NC2=CC=CC=C21',  # Benzimidazolone complex
    'CC1=CC=C(C=C1)C(=O)N2CCN(CC2)C(=O)C3=CC=CC=C3F',  # Fluorinated complex
    'CN(C)C(=O)C1=CC2=C(C=C1)OCO2',  # Methylenedioxy
    'CC(C)N1CCN(CC1)C(=O)C2=CC=C(C=C2)OC3=CC=CC=C3',  # Triple aromatic
    'CN1C2=CC=CC=C2C(=O)N(C1=S)CC3=CC=C(C=C3)Br',  # Brominated thioamide
    'CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=C2)N3CCN(CC3)C'   # Final complex
]

def simulate_docking_affinities(smiles_list, noise_level=1.0):
    """Simulate realistic docking affinities based on molecular properties"""
    print("🎯 Simulando afinidades de docking basadas en propiedades moleculares...")
    
    affinities = []
    
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                affinities.append(0.0)
                continue
            
            # Calculate molecular descriptors that correlate with binding
            mw = rdMolDescriptors.CalcExactMolWt(mol)
            logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            
            # Realistic binding affinity model (simplified)
            # Better drugs tend to have: moderate MW, good LogP, low TPSA, some H-bonds
            base_affinity = -6.0  # Base binding energy
            
            # Molecular weight contribution (optimal around 300-500 Da)
            mw_factor = -0.01 * abs(mw - 400) / 100
            
            # LogP contribution (optimal around 2-4)
            logp_factor = -0.5 * abs(logp - 3) / 2
            
            # TPSA contribution (lower is generally better for membrane permeability)
            tpsa_factor = -0.01 * tpsa / 100
            
            # Hydrogen bonding (some is good, too much is bad)
            hb_factor = -0.3 * max(0, (hbd + hba) - 6)
            
            # Rotatable bonds (flexibility penalty)
            rot_factor = -0.1 * max(0, rot_bonds - 5)
            
            # Calculate predicted affinity
            predicted_affinity = base_affinity + mw_factor + logp_factor + tpsa_factor + hb_factor + rot_factor
            
            # Add realistic noise
            noise = np.random.normal(0, noise_level)
            final_affinity = predicted_affinity + noise
            
            # Clamp to realistic range
            final_affinity = max(-15.0, min(-3.0, final_affinity))
            
            affinities.append(final_affinity)
            
        except Exception as e:
            print(f"❌ Error calculating affinity for {smiles}: {e}")
            affinities.append(0.0)
    
    return affinities

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
    
    # Ensure arrays have same length by taking minimum
    min_len = min(len(fingerprints), len(affinities))
    fingerprints = fingerprints[:min_len]
    affinities = affinities[:min_len]
    
    # Remove zero affinities (failed simulations)
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
    
    # Grid search for best parameters (simplified for demo)
    param_grid = {
        'n_estimators': [100, 200],
        'max_depth': [3, 6],
        'learning_rate': [0.1, 0.2],
        'subsample': [0.8, 1.0]
    }
    
    xgb_model = xgb.XGBRegressor(random_state=42)
    grid_search = GridSearchCV(
        xgb_model, param_grid, cv=3, scoring='neg_mean_squared_error', n_jobs=-1
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
    cv_scores = cross_val_score(best_model, X_train_scaled, y_train, cv=3, scoring='r2')
    metrics['cv_r2_mean'] = cv_scores.mean()
    metrics['cv_r2_std'] = cv_scores.std()
    
    print(f"✅ Modelo entrenado exitosamente")
    print(f"  📊 R² Test: {metrics['test_r2']:.3f}")
    print(f"  📊 RMSE Test: {metrics['test_rmse']:.3f}")
    
    # Save model data for plotting
    model_data = {
        'model': best_model,
        'scaler': scaler,
        'metrics': metrics,
        'X_test': X_test_scaled,
        'y_test': y_test,
        'y_pred_test': y_pred_test,
        'X_train': X_train_scaled,
        'y_train': y_train,
        'y_pred_train': y_pred_train
    }
    
    # Create results directory
    os.makedirs('batch_results_100', exist_ok=True)
    
    with open('batch_results_100/xgboost_model.pkl', 'wb') as f:
        pickle.dump(model_data, f)
    
    return best_model, scaler, metrics

def create_ml_demo_pdf():
    """Create comprehensive ML demonstration PDF"""
    print("📊 Creando demo PDF con análisis completo de Machine Learning...")
    
    # Use the 100 SMILES database
    smiles_list = DRUG_SMILES_100
    
    print("1️⃣  Simulando afinidades de docking...")
    affinities = simulate_docking_affinities(smiles_list)
    
    print("2️⃣  Calculando Morgan fingerprints...")
    fingerprints, valid_smiles = calculate_morgan_fingerprints(smiles_list)
    
    print("3️⃣  Entrenando modelo XGBoost...")
    model, scaler, metrics = train_xgboost_model(fingerprints, affinities)
    
    # Create results DataFrame (ensure same length)
    min_len = min(len(smiles_list), len(affinities))
    results_df = pd.DataFrame({
        'smiles': smiles_list[:min_len],
        'affinity': affinities[:min_len]
    })
    
    # Save results
    results_df.to_csv('batch_results_100/ml_demo_results.csv', index=False)
    
    # Create PDF report
    filename = "batch_results_100/QuickVina_ML_Demo_Report.pdf"
    
    with PdfPages(filename) as pdf:
        # Page 1: Overview and Statistics  
        fig = plt.figure(figsize=(8.5, 11))
        fig.suptitle('QuickVina Machine Learning Demo\n100 Pharmaceutical Compounds Analysis', 
                     fontsize=16, fontweight='bold', y=0.95)
        
        # Summary statistics
        valid_affinities = [a for a in affinities if a != 0.0]
        
        summary_text = f"""MACHINE LEARNING MOLECULAR DOCKING DEMONSTRATION

Dataset: 100 diverse pharmaceutical compounds
Chemical Space: Analgesics, CNS drugs, cardiovascular, antibiotics
Valid Binding Results: {len(valid_affinities)} compounds
Fingerprint Method: Morgan Circular Fingerprints (ECFP)
Fingerprint Dimensions: {fingerprints.shape[1]} bits (radius=2)

BINDING AFFINITY STATISTICS:
Range: {min(valid_affinities):.2f} to {max(valid_affinities):.2f} kcal/mol
Mean: {np.mean(valid_affinities):.2f} ± {np.std(valid_affinities):.2f} kcal/mol
Median: {np.median(valid_affinities):.2f} kcal/mol

MACHINE LEARNING MODEL PERFORMANCE:
Algorithm: XGBoost Gradient Boosting Regressor
Features: Morgan Fingerprints (2048-dimensional binary vectors)
Training Strategy: 80/20 train/test split with 3-fold CV
Feature Engineering: StandardScaler normalization"""
        
        if metrics:
            summary_text += f"""

MODEL PERFORMANCE METRICS:
• R² Score (Test Set): {metrics['test_r2']:.3f}
• Root Mean Square Error: {metrics['test_rmse']:.3f} kcal/mol  
• Mean Absolute Error: {metrics['test_mae']:.3f} kcal/mol
• Cross-Validation R²: {metrics['cv_r2_mean']:.3f} ± {metrics['cv_r2_std']:.3f}

OPTIMIZED HYPERPARAMETERS:
• Number of Estimators: {metrics['best_params']['n_estimators']}
• Maximum Tree Depth: {metrics['best_params']['max_depth']}
• Learning Rate: {metrics['best_params']['learning_rate']}
• Subsample Ratio: {metrics['best_params']['subsample']}

DATASET CHARACTERISTICS:
• Training Samples: {metrics['n_samples']}
• Feature Dimensions: {metrics['n_features']}
• Cross-Validation Score: {metrics['cv_score']:.3f}"""
        else:
            summary_text += "\n\n❌ Model training failed - insufficient data"
        
        fig.text(0.1, 0.75, summary_text, fontsize=9, verticalalignment='top', 
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcyan", alpha=0.9, edgecolor='navy'))
        
        # Top compounds table
        sorted_results = sorted(zip(smiles_list, affinities), key=lambda x: x[1])
        
        table_text = "\nTOP 10 BINDING AFFINITIES:\n"
        table_text += "─" * 50 + "\n"
        table_text += f"{'Rank':<4} {'Affinity (kcal/mol)':<18} {'SMILES':<28}\n"
        table_text += "─" * 50 + "\n"
        
        for i, (smiles, affinity) in enumerate(sorted_results[:10], 1):
            short_smiles = smiles[:25] + "..." if len(smiles) > 25 else smiles
            table_text += f"{i:<4} {affinity:<18.2f} {short_smiles:<28}\n"
        
        fig.text(0.1, 0.35, table_text, fontsize=8, verticalalignment='top', 
                fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.9))
        
        # Save first page
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
        plt.close()
        
        # Page 2: Model Performance Analysis
        if metrics and model is not None:
            fig2 = plt.figure(figsize=(8.5, 11))
            fig2.suptitle('XGBoost Model Performance & Feature Analysis', 
                         fontsize=16, fontweight='bold', y=0.95)
            
            # Load model data
            with open('batch_results_100/xgboost_model.pkl', 'rb') as f:
                model_data = pickle.load(f)
            
            X_test = model_data['X_test']
            y_test = model_data['y_test']
            y_pred_test = model_data['y_pred_test']
            X_train = model_data['X_train']
            y_train = model_data['y_train']
            y_pred_train = model_data['y_pred_train']
            
            # Predicted vs Actual scatter plot
            ax1 = fig2.add_subplot(2, 2, 1)
            ax1.scatter(y_test, y_pred_test, alpha=0.7, color='blue', s=50, label='Test Set')
            ax1.scatter(y_train, y_pred_train, alpha=0.4, color='red', s=30, label='Train Set')
            min_val = min(min(y_test), min(y_train))
            max_val = max(max(y_test), max(y_train))
            ax1.plot([min_val, max_val], [min_val, max_val], 'k--', lw=2, alpha=0.8)
            ax1.set_xlabel('Actual Binding Affinity (kcal/mol)')
            ax1.set_ylabel('Predicted Binding Affinity (kcal/mol)')
            ax1.set_title(f'Model Predictions (R² = {metrics["test_r2"]:.3f})')
            ax1.grid(True, alpha=0.3)
            ax1.legend()
            
            # Residuals analysis
            ax2 = fig2.add_subplot(2, 2, 2)
            residuals_test = y_test - y_pred_test
            residuals_train = y_train - y_pred_train
            ax2.scatter(y_pred_test, residuals_test, alpha=0.7, color='blue', s=50, label='Test')
            ax2.scatter(y_pred_train, residuals_train, alpha=0.4, color='red', s=30, label='Train')
            ax2.axhline(y=0, color='k', linestyle='--', alpha=0.8)
            ax2.set_xlabel('Predicted Binding Affinity (kcal/mol)')
            ax2.set_ylabel('Residuals (kcal/mol)')
            ax2.set_title('Residuals Analysis')
            ax2.grid(True, alpha=0.3)
            ax2.legend()
            
            # Feature importance plot
            ax3 = fig2.add_subplot(2, 1, 2)
            feature_importance = model.feature_importances_
            top_features_idx = np.argsort(feature_importance)[-15:]  # Top 15
            top_importance = feature_importance[top_features_idx]
            
            bars = ax3.barh(range(len(top_features_idx)), top_importance, color='skyblue', alpha=0.8)
            ax3.set_xlabel('Feature Importance (XGBoost)')
            ax3.set_ylabel('Morgan Fingerprint Bit Index')
            ax3.set_title('Top 15 Most Important Molecular Features')
            ax3.set_yticks(range(len(top_features_idx)))
            ax3.set_yticklabels([f'Bit {i}' for i in top_features_idx])
            ax3.grid(True, axis='x', alpha=0.3)
            
            # Add importance values on bars
            for i, (bar, importance) in enumerate(zip(bars, top_importance)):
                ax3.text(bar.get_width() + 0.001, bar.get_y() + bar.get_height()/2, 
                        f'{importance:.3f}', va='center', fontsize=8)
            
            plt.tight_layout()
            pdf.savefig(fig2, bbox_inches='tight', dpi=300)
            plt.close()
        
        # Page 3: Chemical Space Analysis
        fig3 = plt.figure(figsize=(8.5, 11))
        fig3.suptitle('Chemical Space & Molecular Property Analysis', 
                     fontsize=16, fontweight='bold', y=0.95)
        
        # Calculate molecular properties for all compounds
        properties = []
        for smiles in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    props = {
                        'MW': rdMolDescriptors.CalcExactMolWt(mol),
                        'LogP': rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
                        'TPSA': rdMolDescriptors.CalcTPSA(mol),
                        'HBD': rdMolDescriptors.CalcNumHBD(mol),
                        'HBA': rdMolDescriptors.CalcNumHBA(mol),
                        'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol)
                    }
                    properties.append(props)
            except:
                continue
        
        if properties:
            props_df = pd.DataFrame(properties)
            
            # Molecular weight vs LogP
            ax4 = fig3.add_subplot(2, 2, 1)
            scatter = ax4.scatter(props_df['MW'], props_df['LogP'], 
                                c=valid_affinities[:len(props_df)], 
                                cmap='viridis_r', alpha=0.7, s=50)
            ax4.set_xlabel('Molecular Weight (Da)')
            ax4.set_ylabel('LogP (Lipophilicity)')
            ax4.set_title('Chemical Space: MW vs LogP')
            ax4.grid(True, alpha=0.3)
            plt.colorbar(scatter, ax=ax4, label='Binding Affinity (kcal/mol)')
            
            # TPSA vs Affinity
            ax5 = fig3.add_subplot(2, 2, 2)
            ax5.scatter(props_df['TPSA'], valid_affinities[:len(props_df)], 
                       alpha=0.7, color='orange', s=50)
            ax5.set_xlabel('Topological Polar Surface Area (Ų)')
            ax5.set_ylabel('Binding Affinity (kcal/mol)')
            ax5.set_title('TPSA vs Binding Affinity')
            ax5.grid(True, alpha=0.3)
            
            # Property distribution histograms
            ax6 = fig3.add_subplot(2, 2, 3)
            ax6.hist(props_df['MW'], bins=15, alpha=0.7, color='lightblue', edgecolor='black')
            ax6.set_xlabel('Molecular Weight (Da)')
            ax6.set_ylabel('Frequency')
            ax6.set_title('Molecular Weight Distribution')
            ax6.grid(True, alpha=0.3)
            
            ax7 = fig3.add_subplot(2, 2, 4)
            ax7.hist(valid_affinities, bins=15, alpha=0.7, color='lightgreen', edgecolor='black')
            ax7.set_xlabel('Binding Affinity (kcal/mol)')
            ax7.set_ylabel('Frequency')
            ax7.set_title('Binding Affinity Distribution')
            ax7.grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig3, bbox_inches='tight', dpi=300)
        plt.close()
        
    print(f"✅ Demo ML Report creado: {filename}")
    return filename

if __name__ == "__main__":
    try:
        print("🚀 DEMOSTRACIÓN DE MACHINE LEARNING PARA DOCKING")
        print("=" * 70)
        
        report_file = create_ml_demo_pdf()
        
        print(f"\n🎉 ANÁLISIS DE DEMOSTRACIÓN COMPLETADO")
        print(f"📄 Reporte: {report_file}")
        print(f"📏 Ubicación: {os.path.abspath(report_file)}")
        print("\n🔬 DEMOSTRACIÓN COMPLETADA:")
        print("  ✅ 100 compuestos farmacéuticos diversos")
        print("  ✅ Afinidades simuladas basadas en propiedades moleculares")
        print("  ✅ Morgan fingerprints calculados (2048 bits)")
        print("  ✅ Modelo XGBoost entrenado y optimizado")
        print("  ✅ Validación cruzada y métricas de rendimiento")
        print("  ✅ Análisis de espacio químico incluido")
        print("  ✅ Reporte PDF profesional de 3 páginas")
        print("\n🏆 ESTE ES EL WORKFLOW COMPLETO DE ML PARA DOCKING!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc() 