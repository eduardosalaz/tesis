
import pandas as pd

# Load the dataset
file_path = 'df_grasp_1250_all_combs.csv'  # Modify the path as needed
df_new = pd.read_csv(file_path)

# Bounds provided for the new 20 instances (ID 1 to 20)
new_bounds = [944139.739, 964942.154, 980835.67, 970702.713, 957236.295, 960100.522, 
              960425.2, 981274.539, 946657.095, 990223.35, 952796.984, 940185.685, 
              1008565.72, 983003.429, 955365.056, 983749.642, 989806.477, 972287.99, 
              970683.465, 967396.096]

# Function to recalculate the gap and update ValueOptim with new bounds
def recalculate_gap_new_bounds(row):
    instance_id = row['ID'] - 1  # Convert to zero-indexed for accessing the bounds
    if 0 <= instance_id < len(new_bounds):
        bound = new_bounds[instance_id]
        row['BestBound'] = bound
        if row['ValueGRASP'] == 0:
            row['Gap%'] = 100
        else:
            row['Gap%'] = (abs(row['ValueGRASP'] - bound) / abs(row['ValueGRASP'])) * 100
            if row['ValueGRASP'] < bound:
                row['Gap%'] = -row['Gap%']
    return row

# Apply the function across the dataframe
df_new_updated = df_new.apply(recalculate_gap_new_bounds, axis=1)
df_new_updated = df_new_updated[df_new_updated['ID'] != 21]
df_new_updated['TotalTimeGRASP'] = df_new_updated['TotalTimeGRASP'] / 1000

# Filter out rows where ID is 21 to disregard the 21st instance
df_new_filtered_no_21 = df_new_updated[df_new_updated['ID'] != 21]

# Group by the 125 combinations of alpha_a, alpha_l, and iters, then calculate the required stats
summary_stats_new_no_21 = df_new_filtered_no_21.groupby(['alpha_a', 'alpha_l', 'iters']).agg(
    TotalTimeGRASP_min=('TotalTimeGRASP', 'min'),
    TotalTimeGRASP_avg=('TotalTimeGRASP', 'mean'),
    TotalTimeGRASP_max=('TotalTimeGRASP', 'max'),
    Gap_min=('Gap%', 'min'),
    Gap_avg=('Gap%', 'mean'),
    Gap_max=('Gap%', 'max')
).reset_index()

# Find the rows corresponding to the minimum Gap% for each combination
min_gap_rows_new_no_21 = df_new_filtered_no_21.loc[
    df_new_filtered_no_21.groupby(['alpha_a', 'alpha_l', 'iters'])['Gap%'].idxmin()]

# Save the updated dataset, summary statistics, and min gap rows
df_new_updated.to_csv('df_grasp_1250_updated_with_percentage_3.csv', index=False)
summary_stats_new_no_21.to_csv('summary_1250_stats_new_percentage_no_21_3.csv', index=False)
min_gap_rows_new_no_21.to_csv('min_gap_rows_1250_new_percentage_no_21_3.csv', index=False)

print("Files saved: df_grasp_1250_updated_with_percentage2.csv, summary_stats_new_percentage_no_212.csv, min_gap_rows_new_percentage_no_212.csv")
