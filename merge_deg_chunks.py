import pandas as pd
import os

# Set the directory where chunk files are saved
chunk_dir = 'OutputFiles/chunks'

# Get all CSV chunk files
chunk_files = sorted([
    os.path.join(chunk_dir, fname)
    for fname in os.listdir(chunk_dir)
    if fname.endswith('.csv')
])

# Safety check
if not chunk_files:
    raise FileNotFoundError("No chunk CSV files found in 'OutputFiles/chunks/'")

print(f"🔍 Found {len(chunk_files)} chunk files. Merging...")

# Merge all chunk files into one DataFrame
all_chunks = pd.concat([pd.read_csv(f) for f in chunk_files], ignore_index=True)

# Sort by p-value
all_chunks = all_chunks.sort_values(by='p_value')

# Filter for significant DEGs
deg_filtered = all_chunks[
    (all_chunks['p_value'] < 0.05) &
    (all_chunks['log2FC'].abs() > 1)
]

# Save merged results
output_dir = 'OutputFiles'
all_chunks.to_csv(os.path.join(output_dir, 'DEG_full_results.csv'), index=False)
deg_filtered.to_csv(os.path.join(output_dir, 'DEG_significant_filtered.csv'), index=False)

print("✅ Merging complete.")
print("📁 Full results saved to: OutputFiles/DEG_full_results.csv")
print("⭐ Significant DEGs saved to: OutputFiles/DEG_significant_filtered.csv")