#!/usr/bin/env python3

import pandas as pd
import numpy as np

# Read the results
df = pd.read_csv('batch_results/docking_results.csv')

# Sort by best affinity (most negative = strongest binding)
df_sorted = df.sort_values('Best_Affinity_kcal_mol')

# Add ranking
df_sorted['Rank'] = range(1, len(df_sorted) + 1)

# Calculate average runtime
avg_runtime = df_sorted['Runtime_seconds'].mean()
total_runtime = df_sorted['Runtime_seconds'].sum()

print("🧬 QUICKVINA BATCH DOCKING RESULTS")
print("=" * 60)
print(f"📊 Total Ligands Tested: {len(df_sorted)}")
print(f"⏱️  Total Runtime: {total_runtime:.1f} seconds")
print(f"⚡ Average per Ligand: {avg_runtime:.1f} seconds")
print()

# Create formatted table
print("🏆 RANKING BY BINDING AFFINITY")
print("-" * 60)
print(f"{'Rank':<4} {'Ligand Name':<18} {'Best Score':<12} {'2nd Best':<10} {'3rd Best':<10}")
print("-" * 60)

for _, row in df_sorted.iterrows():
    rank = int(row['Rank'])
    name = row['Ligand_Name']
    best = row['Best_Affinity_kcal_mol']
    second = row['Second_Best_kcal_mol']
    third = row['Third_Best_kcal_mol']
    
    # Add medal emojis for top 3
    medal = ""
    if rank == 1:
        medal = "🥇"
    elif rank == 2:
        medal = "🥈"
    elif rank == 3:
        medal = "🥉"
    
    print(f"{rank:<4} {name:<18} {best:<12.1f} {second:<10.1f} {third:<10.1f} {medal}")

print()
print("📈 ANALYSIS:")
print(f"🎯 Best Binder: {df_sorted.iloc[0]['Ligand_Name']} ({df_sorted.iloc[0]['Best_Affinity_kcal_mol']:.1f} kcal/mol)")
print(f"🔄 Worst Binder: {df_sorted.iloc[-1]['Ligand_Name']} ({df_sorted.iloc[-1]['Best_Affinity_kcal_mol']:.1f} kcal/mol)")
print(f"📊 Score Range: {df_sorted['Best_Affinity_kcal_mol'].min():.1f} to {df_sorted['Best_Affinity_kcal_mol'].max():.1f} kcal/mol")
print(f"⚖️  Average Score: {df_sorted['Best_Affinity_kcal_mol'].mean():.1f} kcal/mol")

# Save formatted results
df_sorted.to_csv('batch_results/ranked_results.csv', index=False)
print(f"\n💾 Results saved to: batch_results/ranked_results.csv") 