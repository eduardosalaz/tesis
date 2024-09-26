
import pandas as pd

# Load dataset
df = pd.read_csv('your_dataset_path.csv')

# List of bounds provided for each instance (ID 1 to 20)
bounds = [462006.983, 474508.581, 451829.403, 460361.988, 468430.222, 490028.106, 
          453977.334, 460849.226, 487419.813, 466940.063, 466771.252, 460985.561, 
          468312.438, 459248.658, 479009.770, 462309.765, 462066.303, 454192.089, 
          458314.077, 464814.975]

# Function to recalculate Gap% and update ValueOptim column
def recalculate_gap_safe(row):
    instance_id = row['ID'] - 1  # Convert to zero-indexed for accessing the bounds
    if 0 <= instance_id < len(bounds):
        bound = bounds[instance_id]
        row['ValueOptim'] = bound
        if row['ValueGRASP'] == 0:
            row['Gap%'] = 0
        else:
            row['Gap%'] = abs(row['ValueGRASP'] - bound) / abs(row['ValueGRASP'])
    return row

# Apply recalculation function
df_updated_safe = df.apply(recalculate_gap_safe, axis=1)

# Filter out rows where ValueGRASP is 0
df_filtered = df_updated_safe[df_updated_safe['ValueGRASP'] != 0]

# Group by the 125 combinations of alpha_a, alpha_l, and iters
summary_stats_combinations = df_filtered.groupby(['alpha_a', 'alpha_l', 'iters']).agg(
    TotalTimeGRASP_min=('TotalTimeGRASP', 'min'),
    TotalTimeGRASP_avg=('TotalTimeGRASP', 'mean'),
    TotalTimeGRASP_max=('TotalTimeGRASP', 'max'),
    Gap_min=('Gap%', 'min'),
    Gap_avg=('Gap%', 'mean'),
    Gap_max=('Gap%', 'max')
).reset_index()

# Save the results to a CSV file
summary_stats_combinations.to_csv('summary_stats_combinations.csv', index=False)

print("Summary statistics have been saved to 'summary_stats_combinations.csv'")
