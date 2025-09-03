import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def analyze_single_dataset(df, dataset_name):
    """
    Analyze initialization methods for a single dataset by P value and success status
    """
    
    print(f"\n{'='*60}")
    print(f"ANALYZING {dataset_name.upper()} DATASET")
    print(f"{'='*60}")
    
    # Filter for successful runs only
    successful_data = df[df['status'] == 'success'].copy()
    
    print(f"Total records: {len(df)}")
    print(f"Successful records: {len(successful_data)}")
    print(f"Success rate: {len(successful_data)/len(df)*100:.1f}%")
    print()
    
    print("P values available:")
    print(sorted(df['P'].unique()))
    print()
    
    print("Methods available:")
    print(successful_data['method'].unique())
    print()
    
    return successful_data

def create_method_comparison_single(data, dataset_name, metrics=['init_time', 'prop_split_bus', 'total_time', 'post_ls_obj']):
    """
    Create comprehensive comparison of methods for a single dataset
    """
    
    results = {}
    
    print(f"\n=== DETAILED ANALYSIS FOR {dataset_name.upper()} DATASET ===")
    
    for p_value in sorted(data['P'].unique()):
        p_data = data[data['P'] == p_value]
        
        print(f"\n--- P Value = {p_value} ---")
        print(f"Total instances: {len(p_data)}")
        print(f"Instances per method:")
        print(p_data['method'].value_counts())
        print()
        
        # Summary statistics by method
        summary = p_data.groupby('method')[metrics].agg(['mean', 'std', 'count']).round(4)
        results[p_value] = summary
        
        print("Summary Statistics:")
        print(summary)
        print()
        
        # Find best method for each metric
        best_methods = {}
        for metric in metrics:
            method_means = p_data.groupby('method')[metric].mean()
            if metric in ['init_time', 'total_time']:
                # Lower is better for time metrics
                best_method = method_means.idxmin()
                best_value = method_means.min()
                worst_method = method_means.idxmax()
                worst_value = method_means.max()
            elif metric == 'post_ls_obj':
                # Lower is typically better for objective values
                best_method = method_means.idxmin()
                best_value = method_means.min()
                worst_method = method_means.idxmax()
                worst_value = method_means.max()
            else:  # prop_split_bus - showing both perspectives
                best_method_min = method_means.idxmin()
                best_method_max = method_means.idxmax()
                best_methods[f'{metric}_min'] = (best_method_min, method_means.min())
                best_methods[f'{metric}_max'] = (best_method_max, method_means.max())
                continue
                
            best_methods[metric] = {
                'best': (best_method, best_value),
                'worst': (worst_method, worst_value)
            }
        
        print(f"Best/Worst performing methods for P={p_value}:")
        for metric, values in best_methods.items():
            if isinstance(values, tuple):  # For prop_split_bus min/max
                method, value = values
                print(f"  {metric}: {method} ({value:.4f})")
            else:
                best_method, best_value = values['best']
                worst_method, worst_value = values['worst']
                print(f"  {metric} - Best: {best_method} ({best_value:.4f}), Worst: {worst_method} ({worst_value:.4f})")
        print()
    
    return results

def create_visualizations_single(data, dataset_name, save_plots=True):
    """
    Create visualizations for a single dataset
    """
    
    metrics = ['init_time', 'prop_split_bus', 'total_time', 'post_ls_obj']
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Overall comparison across all P values
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()
    
    for i, metric in enumerate(metrics):
        ax = axes[i]
        
        # Create box plots comparing methods across P values
        sns.boxplot(data=data, x='method', y=metric, hue='P', ax=ax)
        ax.set_title(f'{metric} by Method ({dataset_name} dataset)')
        ax.tick_params(axis='x', rotation=45)
        
        if metric in ['init_time', 'total_time']:
            ax.set_ylabel(f'{metric} (seconds)')
        elif metric == 'post_ls_obj':
            ax.set_ylabel('Objective Value')
            ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    
    plt.suptitle(f'Method Comparison - {dataset_name.upper()} Dataset', fontsize=16)
    plt.tight_layout()
    if save_plots:
        plt.savefig(f'method_comparison_{dataset_name}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create separate plots for each P value
    for p_value in sorted(data['P'].unique()):
        p_data = data[data['P'] == p_value]
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.flatten()
        
        for i, metric in enumerate(metrics):
            ax = axes[i]
            sns.barplot(data=p_data, x='method', y=metric, ax=ax, ci='sd')
            ax.set_title(f'{metric} by Method')
            ax.tick_params(axis='x', rotation=45)
            
            # Add value labels on bars
            for container in ax.containers:
                ax.bar_label(container, fmt='%.3f', rotation=90, fontsize=8)
        
        plt.suptitle(f'{dataset_name.upper()} Dataset - P Value = {p_value}', fontsize=16)
        plt.tight_layout()
        if save_plots:
            plt.savefig(f'method_comparison_{dataset_name}_P{p_value}.png', dpi=300, bbox_inches='tight')
        plt.show()

def statistical_analysis_single(data, dataset_name):
    """
    Perform statistical tests for a single dataset
    """
    
    print(f"\n=== STATISTICAL ANALYSIS - {dataset_name.upper()} DATASET ===")
    
    metrics = ['init_time', 'prop_split_bus', 'total_time', 'post_ls_obj']
    
    for p_value in sorted(data['P'].unique()):
        p_data = data[data['P'] == p_value]
        
        print(f"\n--- P Value = {p_value} ---")
        
        for metric in metrics:
            print(f"\n{metric}:")
            
            # Group data by method
            method_groups = [group[metric].values for name, group in p_data.groupby('method')]
            method_names = [name for name, group in p_data.groupby('method')]
            
            if len(method_groups) > 2:
                # Perform ANOVA if more than 2 groups
                try:
                    f_stat, p_val = stats.f_oneway(*method_groups)
                    print(f"  ANOVA F-statistic: {f_stat:.4f}, p-value: {p_val:.4f}")
                    
                    if p_val < 0.05:
                        print("  *** Significant differences detected between methods! ***")
                        
                        # Perform pairwise t-tests for significant results
                        print("  Pairwise comparisons:")
                        for i, method1 in enumerate(method_names):
                            for j, method2 in enumerate(method_names[i+1:], i+1):
                                t_stat, t_p = stats.ttest_ind(method_groups[i], method_groups[j])
                                significance = "***" if t_p < 0.001 else "**" if t_p < 0.01 else "*" if t_p < 0.05 else ""
                                print(f"    {method1} vs {method2}: p={t_p:.4f} {significance}")
                    else:
                        print("  No significant differences detected.")
                except Exception as e:
                    print(f"  ANOVA failed: {e}")

def generate_recommendations_single(data, dataset_name):
    """
    Generate recommendations for best initialization method for a single dataset
    """
    
    print(f"\n{'='*60}")
    print(f"RECOMMENDATIONS FOR {dataset_name.upper()} DATASET")
    print(f"{'='*60}")
    
    recommendations = {}
    
    # You can adjust these weights based on your priorities
    weights = {
        'init_time': -0.3,      # Negative because lower is better
        'total_time': -0.4,     # Negative because lower is better  
        'post_ls_obj': -0.2,    # Negative because lower is better
        'prop_split_bus': 0.1   # Adjust based on whether higher or lower split is better
    }
    
    for p_value in sorted(data['P'].unique()):
        p_data = data[data['P'] == p_value]
        
        print(f"\n--- P Value = {p_value} ---")
        
        # Calculate composite scores
        method_scores = {}
        detailed_scores = {}
        
        for method in p_data['method'].unique():
            method_data = p_data[p_data['method'] == method]
            score = 0
            method_details = {}
            
            for metric, weight in weights.items():
                if metric in method_data.columns:
                    metric_mean = method_data[metric].mean()
                    # Normalize to 0-1 scale
                    metric_min = p_data[metric].min()
                    metric_max = p_data[metric].max()
                    if metric_max != metric_min:
                        norm_value = (metric_mean - metric_min) / (metric_max - metric_min)
                    else:
                        norm_value = 0
                    
                    weighted_score = weight * norm_value
                    score += weighted_score
                    method_details[metric] = {
                        'raw_value': metric_mean,
                        'normalized': norm_value,
                        'weighted': weighted_score
                    }
            
            method_scores[method] = score
            detailed_scores[method] = method_details
        
        # Sort methods by score
        sorted_methods = sorted(method_scores.items(), key=lambda x: x[1], reverse=True)
        best_method = sorted_methods[0][0]
        recommendations[p_value] = best_method
        
        print(f"Recommended method: {best_method}")
        print("\nMethod rankings (higher score is better):")
        for i, (method, score) in enumerate(sorted_methods, 1):
            print(f"  {i}. {method}: {score:.4f}")
        
        # Show detailed breakdown for top method
        print(f"\nDetailed breakdown for {best_method}:")
        for metric, details in detailed_scores[best_method].items():
            print(f"  {metric}: {details['raw_value']:.4f} (weighted contribution: {details['weighted']:.4f})")
    
    return recommendations

def export_results_single(data, results, recommendations, dataset_name):
    """
    Export results to Excel file
    """
    
    with pd.ExcelWriter(f'{dataset_name}_analysis_results.xlsx') as writer:
        # Export raw data
        data.to_excel(writer, sheet_name='Raw_Data', index=False)
        
        # Export summary statistics for each P value
        for p_value, summary in results.items():
            sheet_name = f'Summary_P{p_value}'.replace('.', '_')
            summary.to_excel(writer, sheet_name=sheet_name)
        
        # Export recommendations
        rec_df = pd.DataFrame.from_dict(recommendations, orient='index', columns=['Recommended_Method'])
        rec_df.index.name = 'P_Value'
        rec_df.to_excel(writer, sheet_name='Recommendations')
    
    print(f"Results exported to {dataset_name}_analysis_results.xlsx")

# Main function to analyze a single dataset
def analyze_dataset(file_path, dataset_name):
    """
    Complete analysis pipeline for a single dataset
    """
    
    # Load dataset
    print(f"Loading {dataset_name} dataset from {file_path}...")
    df = pd.read_csv(file_path)
    
    # Analyze the data
    successful_data = analyze_single_dataset(df, dataset_name)
    
    # Create detailed comparison
    results = create_method_comparison_single(successful_data, dataset_name)
    
    # Create visualizations
    create_visualizations_single(successful_data, dataset_name)
    
    # Perform statistical analysis
    statistical_analysis_single(successful_data, dataset_name)
    
    # Generate recommendations
    recommendations = generate_recommendations_single(successful_data, dataset_name)
    
    # Export results
    export_results_single(successful_data, results, recommendations, dataset_name)
    
    return successful_data, results, recommendations

# Usage - Run each dataset independently
if __name__ == "__main__":
    
    # Analyze small dataset
    print("STARTING ANALYSIS...")
    
    small_data, small_results, small_recommendations = analyze_dataset(
        "small_constructive_data.csv", "small"  # Replace with your actual file path
    )
    
    medium_data, medium_results, medium_recommendations = analyze_dataset(
        "medium_constructive_data.csv", "medium"  # Replace with your actual file path
    )
    
    large_data, large_results, large_recommendations = analyze_dataset(
        "large_constructive_data.csv", "large"  # Replace with your actual file path
    )
    
    # Final summary across all datasets
    print(f"\n{'='*80}")
    print("FINAL SUMMARY ACROSS ALL DATASETS")
    print(f"{'='*80}")
    
    all_recommendations = {
        'small': small_recommendations,
        'medium': medium_recommendations,
        'large': large_recommendations
    }
    
    for dataset, recs in all_recommendations.items():
        print(f"\n{dataset.upper()} dataset recommendations:")
        for p_value, method in recs.items():
            print(f"  P={p_value}: {method}")
    
    print(f"\nAnalysis complete! Check the individual Excel files and plots for detailed results.")