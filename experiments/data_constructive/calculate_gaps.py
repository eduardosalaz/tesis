import pandas as pd

# Read the heuristics CSV file (assumes it's named 'heuristics.csv')
heuristics_df = pd.read_csv('small_constructive_data.csv')

# Read the gurobi results CSV file (assumes it's named 'gurobi_results.csv')
gurobi_df = pd.read_csv('results_gurobi_small.csv')

# Merge heuristics with gurobi results using instance_number and P columns
merged_df = heuristics_df.merge(gurobi_df, left_on=['instance_number', 'P'], right_on=['instance_id', 'p'], how='left')

# Calculate optimality gaps
def calculate_gap(heuristic_obj, gurobi_bound, status):
    # If status indicates error, return gap of 1.0 (100%)
    if pd.isna(status) or status == 'error':
        return 1.0
    if pd.isna(heuristic_obj) or pd.isna(gurobi_bound) or gurobi_bound == 0:
        return 1.0
    return (heuristic_obj - gurobi_bound) / gurobi_bound

merged_df['post_repair_gap'] = merged_df.apply(lambda row: calculate_gap(row['post_repair_obj'], row['best_bound'], row['status']), axis=1)
merged_df['post_ls_gap'] = merged_df.apply(lambda row: calculate_gap(row['post_ls_obj'], row['best_bound'], row['status']), axis=1)

# Check if heuristic gaps beat gurobi gap (smaller gap is better)
merged_df['post_repair_beats_gurobi'] = merged_df['post_repair_gap'] < merged_df['optimality_gap']
merged_df['post_ls_beats_gurobi'] = merged_df['post_ls_gap'] < merged_df['optimality_gap']

# Keep all original columns from heuristics table and add new columns
final_df = merged_df.copy()

# Rename gurobi columns to avoid confusion
final_df = final_df.rename(columns={
    'best_value': 'gurobi_best_value',
    'best_bound': 'gurobi_best_bound', 
    'optimality_gap': 'gurobi_gap',
    'time': 'gurobi_time'
})

# Drop the duplicate columns from the merge (instance_id and p from gurobi table)
final_df = final_df.drop(columns=['instance_id', 'p'], errors='ignore')

# Sort by instance_number, P, and method
final_df = final_df.sort_values(['instance_number', 'P', 'method'])

# Save to CSV
final_df.to_csv('combined_results_small.csv', index=False)

print(f"Processed {len(final_df)} rows and saved to combined_results.csv")
print(f"Methods found: {sorted(final_df['method'].unique())}")
print(f"Instance numbers: {sorted(final_df['instance_number'].unique())}")
print(f"P values: {sorted(final_df['P'].unique())}")