import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_and_prepare_data(file_path, dataset_name):
    """
    Load dataset and prepare for analysis
    """
    print(f"\nLoading {dataset_name} dataset from {file_path}...")
    df = pd.read_csv(file_path)
    
    # Filter for successful runs only
    df_success = df[df['status'] == 'success'].copy()
    
    # Calculate average gap if not present (using post_ls_gap as final gap)
    if 'final_gap' not in df_success.columns:
        df_success['final_gap'] = df_success['post_ls_gap'] * 100
    
    print(f"Dataset shape: {df_success.shape}")
    print(f"Success rate: {len(df_success)/len(df)*100:.1f}%")
    print(f"Methods: {df_success['method'].unique()}")
    print(f"P values: {sorted(df_success['P'].unique())}")
    print(f"Instance numbers per instance_id: {df_success.groupby('instance_id_x')['instance_number'].nunique().unique()}")
    
    return df_success

def analyze_by_instance_aggregation(data, dataset_name):
    """
    Analyze data properly accounting for repeated instances
    Each instance is solved 3 times (instance_number 1, 2, 3)
    """
    print(f"\n{'='*60}")
    print(f"INSTANCE-LEVEL ANALYSIS - {dataset_name.upper()}")
    print(f"{'='*60}")
    
    # Key metrics to analyze
    metrics = {
        'init_time': 'Initialization Time (s)',
        'prop_split_bus': 'Proportion Split BUs',
        'final_gap': 'Optimality Gap',
        'total_time': 'Total Time (s)',
        'post_ls_obj': 'Final Objective Value'
    }
    
    results = {}
    
    for p_value in sorted(data['P'].unique()):
        print(f"\n--- P = {p_value} ---")
        p_data = data[data['P'] == p_value]
        
        # Group by instance_id and method, then average across the 3 runs
        aggregated = p_data.groupby(['instance_id_x', 'method'])[list(metrics.keys())].mean().reset_index()
        
        # Now group by method for statistics
        method_stats = aggregated.groupby('method')[list(metrics.keys())].agg(['mean', 'std', 'count'])
        
        results[p_value] = {
            'raw_data': p_data,
            'aggregated_data': aggregated,
            'statistics': method_stats
        }
        
        print(f"Methods compared: {aggregated['method'].unique()}")
        print(f"Instances per method: {aggregated.groupby('method').size().to_dict()}")
        print("\nAverage Performance by Method:")
        print(method_stats.round(4))
    
    return results

def create_clean_visualizations(results, dataset_name, save_plots=True):
    """
    Create clean visualizations with just the essential 4 plots
    """
    # Set style
    sns.set_style("whitegrid")
    colors = sns.color_palette("husl", n_colors=6)
    
    # For each P value, create the 4 essential plots
    for p_value, p_results in results.items():
        aggregated = p_results['aggregated_data']
        
        # Create 2x2 subplot with the 4 key metrics
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'{dataset_name.upper()} Dataset - P = {p_value}', fontsize=14, fontweight='bold')
        
        # 1. Initialization Time
        ax1 = axes[0, 0]
        sns.boxplot(data=aggregated, x='method', y='init_time', ax=ax1, palette=colors)
        ax1.set_title('Initialization Time', fontsize=11, fontweight='bold')
        ax1.set_xlabel('')
        ax1.set_ylabel('Time (seconds)')
        ax1.tick_params(axis='x', rotation=45)
        
        # 2. Optimality Gap
        ax2 = axes[0, 1]
        sns.boxplot(data=aggregated, x='method', y='final_gap', ax=ax2, palette=colors)
        ax2.set_title('Optimality Gap', fontsize=11, fontweight='bold')
        ax2.set_xlabel('')
        ax2.set_ylabel('Gap (proportion)')
        ax2.tick_params(axis='x', rotation=45)
        
        # 3. Proportion of Split BUs
        ax3 = axes[1, 0]
        sns.boxplot(data=aggregated, x='method', y='prop_split_bus', ax=ax3, palette=colors)
        ax3.set_title('Proportion of Split BUs', fontsize=11, fontweight='bold')
        ax3.set_xlabel('')
        ax3.set_ylabel('Proportion')
        ax3.tick_params(axis='x', rotation=45)
        
        # 4. Total Time
        ax4 = axes[1, 1]
        sns.boxplot(data=aggregated, x='method', y='total_time', ax=ax4, palette=colors)
        ax4.set_title('Total Solution Time', fontsize=11, fontweight='bold')
        ax4.set_xlabel('')
        ax4.set_ylabel('Time (seconds)')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        if save_plots:
            plt.savefig(f'{dataset_name}_P{p_value}_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Create separate time-quality trade-off plot
        create_tradeoff_plot(aggregated, p_value, dataset_name, save_plots)

def create_tradeoff_plot(aggregated, p_value, dataset_name, save_plots=True):
    """
    Create a simple time vs quality trade-off plot
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    # Plot each method
    for method in aggregated['method'].unique():
        method_data = aggregated[aggregated['method'] == method]
        ax.scatter(method_data['init_time'], method_data['final_gap'], 
                  label=method, alpha=0.7, s=80)
    
    ax.set_xlabel('Initialization Time (seconds)', fontsize=11)
    ax.set_ylabel('Optimality Gap', fontsize=11)
    ax.set_title(f'Time-Quality Trade-off - {dataset_name.upper()} (P = {p_value})', 
                fontsize=12, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_plots:
        plt.savefig(f'{dataset_name}_P{p_value}_tradeoff.png', dpi=300, bbox_inches='tight')
    plt.show()

def perform_statistical_tests(results, dataset_name):
    """
    Perform robust statistical tests accounting for repeated measures
    """
    print(f"\n{'='*60}")
    print(f"STATISTICAL ANALYSIS - {dataset_name.upper()}")
    print(f"{'='*60}")
    
    metrics_to_test = ['init_time', 'final_gap', 'prop_split_bus']
    
    for p_value, p_results in results.items():
        print(f"\n--- P = {p_value} ---")
        aggregated = p_results['aggregated_data']
        
        methods = aggregated['method'].unique()
        
        if len(methods) < 2:
            print("Not enough methods to compare")
            continue
        
        for metric in metrics_to_test:
            print(f"\n{metric.replace('_', ' ').upper()}:")
            
            # Prepare data for each method
            method_data = {}
            for method in methods:
                method_data[method] = aggregated[aggregated['method'] == method][metric].values
            
            # Check if we have enough data
            if all(len(data) > 1 for data in method_data.values()):
                # Perform Kruskal-Wallis test (non-parametric alternative to ANOVA)
                groups = list(method_data.values())
                h_stat, p_val = stats.kruskal(*groups)
                
                print(f"  Kruskal-Wallis H-statistic: {h_stat:.4f}, p-value: {p_val:.4f}")
                
                if p_val < 0.05:
                    print("  *** Significant differences detected! ***")
                    
                    # Perform pairwise Mann-Whitney U tests
                    print("\n  Pairwise comparisons (Mann-Whitney U):")
                    method_names = list(method_data.keys())
                    
                    for i, method1 in enumerate(method_names):
                        for j, method2 in enumerate(method_names[i+1:], i+1):
                            u_stat, p_val_pair = stats.mannwhitneyu(
                                method_data[method1], 
                                method_data[method2],
                                alternative='two-sided'
                            )
                            
                            # Apply Bonferroni correction
                            n_comparisons = len(method_names) * (len(method_names) - 1) / 2
                            adjusted_p = p_val_pair * n_comparisons
                            
                            significance = ""
                            if adjusted_p < 0.001:
                                significance = "***"
                            elif adjusted_p < 0.01:
                                significance = "**"
                            elif adjusted_p < 0.05:
                                significance = "*"
                            
                            print(f"    {method1} vs {method2}:")
                            print(f"      p-value (adjusted): {adjusted_p:.4f} {significance}")
                else:
                    print("  No significant differences detected")

def calculate_pareto_frontier(results, dataset_name, save_plots=True):
    """
    Calculate and visualize Pareto frontier for time vs quality trade-off
    """
    print(f"\n{'='*60}")
    print(f"PARETO FRONTIER ANALYSIS - {dataset_name.upper()}")
    print(f"{'='*60}")
    
    for p_value, p_results in results.items():
        print(f"\n--- P = {p_value} ---")
        
        stats = p_results['statistics']
        
        # Extract mean values
        methods = stats.index
        init_times = stats['init_time']['mean'].values
        gaps = stats['final_gap']['mean'].values
        
        # Find Pareto optimal methods
        pareto_methods = []
        for i, method in enumerate(methods):
            is_dominated = False
            for j in range(len(methods)):
                if i != j:
                    # Check if method j dominates method i
                    if init_times[j] <= init_times[i] and gaps[j] <= gaps[i]:
                        if init_times[j] < init_times[i] or gaps[j] < gaps[i]:
                            is_dominated = True
                            break
            if not is_dominated:
                pareto_methods.append(method)
        
        print(f"Pareto optimal methods: {pareto_methods}")
        print(f"Dominated methods: {[m for m in methods if m not in pareto_methods]}")
        
        # Create Pareto frontier plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 7))
        
        # Plot all methods with clear distinction
        for i, method in enumerate(methods):
            if method in pareto_methods:
                # Pareto optimal - filled circles
                ax.scatter(init_times[i], gaps[i], label=f'{method} (Pareto)', 
                          s=150, marker='o', alpha=0.8, edgecolors='black', linewidth=2)
            else:
                # Dominated - triangles
                ax.scatter(init_times[i], gaps[i], label=f'{method} (dominated)', 
                          s=100, marker='^', alpha=0.5)
            
            # Annotate
            ax.annotate(method, (init_times[i], gaps[i]), 
                       xytext=(5, 5), textcoords='offset points', fontsize=10,
                       fontweight='bold' if method in pareto_methods else 'normal')
        
        # Connect Pareto optimal points
        if len(pareto_methods) > 1:
            pareto_points = [(init_times[i], gaps[i]) for i, m in enumerate(methods) 
                           if m in pareto_methods]
            pareto_points.sort(key=lambda x: x[0])
            
            x_pareto = [p[0] for p in pareto_points]
            y_pareto = [p[1] for p in pareto_points]
            ax.plot(x_pareto, y_pareto, 'g--', alpha=0.5, linewidth=1.5)
        
        ax.set_xlabel('Mean Initialization Time (seconds)', fontsize=12)
        ax.set_ylabel('Mean Optimality Gap', fontsize=12)
        ax.set_title(f'Pareto Frontier - {dataset_name.upper()} (P = {p_value})', 
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', framealpha=0.9)
        
        # Add explanation
        ax.text(0.02, 0.98, 
               'Lower is better for both axes\n○ = Pareto optimal\n△ = Dominated',
               transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        if save_plots:
            plt.savefig(f'{dataset_name}_P{p_value}_pareto.png', dpi=300, bbox_inches='tight')
        plt.show()

def generate_recommendations(results, dataset_name):
    """
    Generate context-aware recommendations
    """
    print(f"\n{'='*60}")
    print(f"RECOMMENDATIONS - {dataset_name.upper()}")
    print(f"{'='*60}")
    
    recommendations = {}
    
    for p_value, p_results in results.items():
        print(f"\n--- P = {p_value} ---")
        
        stats = p_results['statistics']
        
        # Different recommendation scenarios
        scenarios = {
            'speed_critical': {
                'weights': {'init_time': -0.7, 'final_gap': -0.2, 'prop_split_bus': -0.1},
                'description': 'When initialization speed is critical'
            },
            'quality_critical': {
                'weights': {'init_time': -0.1, 'final_gap': -0.7, 'prop_split_bus': -0.2},
                'description': 'When solution quality is most important'
            },
            'balanced': {
                'weights': {'init_time': -0.4, 'final_gap': -0.4, 'prop_split_bus': -0.2},
                'description': 'Balanced trade-off between speed and quality'
            }
        }
        
        recommendations[p_value] = {}
        
        for scenario_name, scenario in scenarios.items():
            print(f"\n  {scenario['description']}:")
            
            # Calculate scores
            method_scores = {}
            for method in stats.index:
                score = 0
                for metric, weight in scenario['weights'].items():
                    if metric in ['init_time', 'final_gap', 'prop_split_bus']:
                        # Normalize (0-1, where 0 is best)
                        value = stats.loc[method, (metric, 'mean')]
                        min_val = stats[(metric, 'mean')].min()
                        max_val = stats[(metric, 'mean')].max()
                        
                        if max_val != min_val:
                            normalized = (value - min_val) / (max_val - min_val)
                        else:
                            normalized = 0
                        
                        score += weight * normalized
                
                method_scores[method] = score
            
            # Sort by score (higher is better since weights are negative)
            sorted_methods = sorted(method_scores.items(), key=lambda x: x[1], reverse=True)
            best_method = sorted_methods[0][0]
            
            recommendations[p_value][scenario_name] = best_method
            
            print(f"    Recommended: {best_method}")
            print(f"    Ranking:")
            for rank, (method, score) in enumerate(sorted_methods, 1):
                print(f"      {rank}. {method}: score={score:.4f}")
    
    return recommendations

def export_comprehensive_results(results, recommendations, dataset_name):
    """
    Export all results to Excel with multiple sheets
    """
    with pd.ExcelWriter(f'{dataset_name}_analysis_results.xlsx') as writer:
        
        # Sheet 1: Summary statistics for all P values
        summary_all = []
        for p_value, p_results in results.items():
            stats = p_results['statistics']
            stats_flat = stats.stack(level=0).reset_index()
            stats_flat.columns = ['method', 'metric', 'mean', 'std', 'count']
            stats_flat['P'] = p_value
            summary_all.append(stats_flat)
        
        if summary_all:
            pd.concat(summary_all).to_excel(writer, sheet_name='Summary_Statistics', index=False)
        
        # Sheet 2-N: Detailed results for each P value
        for p_value, p_results in results.items():
            sheet_name = f'P_{p_value}_Details'
            p_results['aggregated_data'].to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Sheet N+1: Recommendations
        rec_data = []
        for p_value, scenarios in recommendations.items():
            for scenario, method in scenarios.items():
                rec_data.append({
                    'P_value': p_value,
                    'Scenario': scenario,
                    'Recommended_Method': method
                })
        
        pd.DataFrame(rec_data).to_excel(writer, sheet_name='Recommendations', index=False)
        
        # Sheet N+2: Statistical test results (create summary)
        test_summary = []
        for p_value, p_results in results.items():
            aggregated = p_results['aggregated_data']
            for metric in ['init_time', 'final_gap', 'prop_split_bus']:
                methods = aggregated['method'].unique()
                if len(methods) > 1:
                    groups = [aggregated[aggregated['method'] == m][metric].values for m in methods]
                    if all(len(g) > 1 for g in groups):
                        h_stat, p_val = stats.kruskal(*groups)
                        test_summary.append({
                            'P_value': p_value,
                            'Metric': metric,
                            'Test': 'Kruskal-Wallis',
                            'Statistic': h_stat,
                            'p_value': p_val,
                            'Significant': p_val < 0.05
                        })
        
        if test_summary:
            pd.DataFrame(test_summary).to_excel(writer, sheet_name='Statistical_Tests', index=False)
    
    print(f"\nResults exported to {dataset_name}_analysis_results.xlsx")

# Main analysis function
def analyze_dataset_comprehensive(file_path, dataset_name):
    """
    Complete analysis pipeline for a dataset
    """
    # Load and prepare data
    data = load_and_prepare_data(file_path, dataset_name)
    
    # Analyze with proper instance aggregation
    results = analyze_by_instance_aggregation(data, dataset_name)
    
    # Create clean visualizations
    create_clean_visualizations(results, dataset_name)
    
    # Statistical testing
    perform_statistical_tests(results, dataset_name)
    
    # Pareto frontier analysis
    calculate_pareto_frontier(results, dataset_name)
    
    # Generate recommendations
    recommendations = generate_recommendations(results, dataset_name)
    
    # Export results
    export_comprehensive_results(results, recommendations, dataset_name)
    
    return results, recommendations

# Example usage
if __name__ == "__main__":
    # Analyze each dataset separately
    datasets = [
        ('combined_results_small.csv', 'small'),
        ('combined_results_medium.csv', 'medium'),
        ('combined_results_large.csv', 'large')
    ]
    
    all_results = {}
    all_recommendations = {}
    
    for file_path, dataset_name in datasets:
        try:
            print(f"\n{'='*80}")
            print(f"ANALYZING {dataset_name.upper()} DATASET")
            print(f"{'='*80}")
            
            results, recommendations = analyze_dataset_comprehensive(file_path, dataset_name)
            all_results[dataset_name] = results
            all_recommendations[dataset_name] = recommendations
            
        except FileNotFoundError:
            print(f"Warning: {file_path} not found. Skipping {dataset_name} dataset.")
        except Exception as e:
            print(f"Error analyzing {dataset_name}: {e}")
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")