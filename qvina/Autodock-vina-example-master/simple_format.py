#!/usr/bin/env python3

import csv

# Read the results
results = []
with open('batch_results/docking_results.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        results.append({
            'name': row['Ligand_Name'],
            'best': float(row['Best_Affinity_kcal_mol']),
            'second': float(row['Second_Best_kcal_mol']),
            'third': float(row['Third_Best_kcal_mol']),
            'runtime': float(row['Runtime_seconds'])
        })

# Sort by best affinity (most negative = strongest binding)
results.sort(key=lambda x: x['best'])

# Calculate statistics
total_runtime = sum(r['runtime'] for r in results)
avg_runtime = total_runtime / len(results)
best_scores = [r['best'] for r in results]
min_score = min(best_scores)
max_score = max(best_scores)
avg_score = sum(best_scores) / len(best_scores)

print("🧬 QUICKVINA BATCH DOCKING RESULTS")
print("=" * 60)
print(f"📊 Total Ligands Tested: {len(results)}")
print(f"⏱️  Total Runtime: {total_runtime:.1f} seconds")
print(f"⚡ Average per Ligand: {avg_runtime:.1f} seconds")
print()

# Create formatted table
print("🏆 RANKING BY BINDING AFFINITY")
print("-" * 70)
print(f"{'Rank':<4} {'Ligand Name':<18} {'Best Score':<12} {'2nd Best':<10} {'3rd Best':<10} {'Medal'}")
print("-" * 70)

for i, result in enumerate(results, 1):
    name = result['name']
    best = result['best']
    second = result['second']
    third = result['third']
    
    # Add medal emojis for top 3
    medal = ""
    if i == 1:
        medal = "🥇"
    elif i == 2:
        medal = "🥈"
    elif i == 3:
        medal = "🥉"
    
    print(f"{i:<4} {name:<18} {best:<12.1f} {second:<10.1f} {third:<10.1f} {medal}")

print()
print("📈 ANALYSIS:")
print(f"🎯 Best Binder: {results[0]['name']} ({results[0]['best']:.1f} kcal/mol)")
print(f"🔄 Worst Binder: {results[-1]['name']} ({results[-1]['best']:.1f} kcal/mol)")
print(f"📊 Score Range: {min_score:.1f} to {max_score:.1f} kcal/mol")
print(f"⚖️  Average Score: {avg_score:.1f} kcal/mol")
print()
print("💡 NOTE: More negative scores indicate stronger binding affinity")
print("💾 Raw data available in: batch_results/docking_results.csv") 