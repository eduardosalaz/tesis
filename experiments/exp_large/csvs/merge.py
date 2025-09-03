import pandas as pd
import glob
import os

def merge_csvs_to_single_file(input_folder_path, output_file_path, file_pattern="*.csv"):
    """
    Merge multiple CSV files into a single CSV file.
    
    Args:
        input_folder_path: Path to folder containing CSV files
        output_file_path: Path for the merged output CSV file
        file_pattern: Pattern to match CSV files (default: "*.csv")
    """
    
    # Get list of all CSV files in the folder
    csv_files = glob.glob(os.path.join(input_folder_path, file_pattern))
    
    if not csv_files:
        print(f"No CSV files found in {input_folder_path}")
        return
    
    print(f"Found {len(csv_files)} CSV files to merge:")
    for file in csv_files:
        print(f"  - {os.path.basename(file)}")
    
    # Read and combine all CSV files
    dataframes = []
    
    for file in csv_files:
        try:
            df = pd.read_csv(file)
            print(f"Loaded {file}: {len(df)} rows, {len(df.columns)} columns")
            dataframes.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    if not dataframes:
        print("No valid CSV files were loaded.")
        return
    
    # Merge all dataframes
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    print(f"\nMerged dataset: {len(merged_df)} total rows, {len(merged_df.columns)} columns")
    
    # Write to output file
    merged_df.to_csv(output_file_path, index=False)
    print(f"Merged data saved to: {output_file_path}")
    
    return merged_df

# Usage example:
if __name__ == "__main__":
    # Specify your input folder and output file
    input_folder = "."  # Change this to your folder path
    output_file = "merged_data.csv"          # Change this to your desired output filename
    
    # Merge the files
    merged_data = merge_csvs_to_single_file(input_folder, output_file)
    
    # Optional: Display basic info about the merged dataset
    if merged_data is not None:
        print("\nDataset Info:")
        print(f"Shape: {merged_data.shape}")
        print(f"Columns: {list(merged_data.columns)}")
        print(f"\nFirst few rows:")
        print(merged_data.head())