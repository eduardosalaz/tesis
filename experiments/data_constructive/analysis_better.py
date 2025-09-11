import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import kruskal, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

def enhanced_statistical_analysis(data, dataset_name):
    """
    Enhanced statistical analysis with focus on method rankings and tradeoffs
    """
    
    print(f"\n{'='*80}")
    print(f"ENHANCED STATISTICAL ANALYSIS - {dataset_name.upper()} DATASET")
    print(f"{'='*80}")
    
    key_metrics = ['init_time', 'total_time', 'post_ls_gap']
    results_summary = {}
    
    for p_value in sorted(data['P'].unique()):
        p_data = data[data['P'] == p_value]
        
        print(f"\n{'-'*60}")
        print(f"P Value = {p_value}")
        print(f"{'-'*60}")
        
        p_results = {}
        
        for metric in key_metrics:
            print(f"\n📊 ANALYSIS FOR {metric.upper()}:")
            print("=" * 50)
            
            # Get method performance
            method_stats = p_data.groupby('method')[metric].agg(['count', 'mean', 'std', 'median', 'min', 'max']).round(4)
            method_stats = method_stats.sort_values('mean')  # Sort by mean performance
            
            print("\n📈 Method Performance (sorted by mean):")
            print(method_stats)
            
            # Statistical significance testing
            method_groups = [group[metric].values for name, group in p_data.groupby('method')]
            method_names = list(p_data.groupby('method').groups.keys())
            
            if len(method_groups) > 2:
                # Use Kruskal-Wallis (non-parametric) as it's more robust
                h_stat, p_val = kruskal(*method_groups)
                print(f"\n🔬 Kruskal-Wallis Test: H = {h_stat:.4f}, p-value = {p_val:.4f}")
                
                if p_val < 0.05:
                    print("✅ SIGNIFICANT DIFFERENCES DETECTED!")
                    
                    # Pairwise comparisons with effect size
                    print("\n🔍 Pairwise Comparisons (Mann-Whitney U):")
                    pairwise_results = []
                    
                    for i, method1 in enumerate(method_names):
                        group1 = [group[metric].values for name, group in p_data.groupby('method') if name == method1][0]
                        for j, method2 in enumerate(method_names[i+1:], i+1):
                            group2 = [group[metric].values for name, group in p_data.groupby('method') if name == method2][0]
                            
                            u_stat, u_p = mannwhitneyu(group1, group2, alternative='two-sided')
                            
                            # Calculate effect size (Cohen's d approximation)
                            mean1, mean2 = np.mean(group1), np.mean(group2)
                            std1, std2 = np.std(group1, ddof=1), np.std(group2, ddof=1)
                            pooled_std = np.sqrt(((len(group1)-1)*std1**2 + (len(group2)-1)*std2**2) / (len(group1)+len(group2)-2))
                            cohens_d = abs(mean1 - mean2) / pooled_std if pooled_std > 0 else 0
                            
                            # Effect size interpretation
                            if cohens_d < 0.2:
                                effect = "negligible"
                            elif cohens_d < 0.5:
                                effect = "small"
                            elif cohens_d < 0.8:
                                effect = "medium"
                            else:
                                effect = "large"
                            
                            significance = "***" if u_p < 0.001 else "**" if u_p < 0.01 else "*" if u_p < 0.05 else ""
                            
                            better_method = method1 if mean1 < mean2 else method2  # Assuming lower is better for all metrics
                            improvement = abs((mean1 - mean2) / max(mean1, mean2)) * 100
                            
                            print(f"   {method1} vs {method2}: p={u_p:.4f} {significance}")
                            print(f"      → Better method: {better_method} ({improvement:.1f}% improvement)")
                            print(f"      → Effect size: {cohens_d:.3f} ({effect})")
                            
                            pairwise_results.append({
                                'method1': method1, 'method2': method2, 'p_value': u_p,
                                'better_method': better_method, 'improvement_pct': improvement,
                                'effect_size': cohens_d, 'effect_magnitude': effect
                            })
                    
                    # Find the best method overall
                    best_method = method_stats.index[0]  # First in sorted list (lowest mean)
                    best_mean = method_stats.loc[best_method, 'mean']
                    
                    print(f"\n🏆 BEST METHOD: {best_method} (mean = {best_mean:.4f})")
                    
                    # Show how much better the best method is compared to others
                    print(f"\n📊 Performance compared to best method ({best_method}):")
                    for method in method_stats.index[1:]:  # Skip the best method
                        other_mean = method_stats.loc[method, 'mean']
                        worse_by = ((other_mean - best_mean) / best_mean) * 100
                        print(f"   {method}: {worse_by:.1f}% worse")
                
                else:
                    print("❌ No significant differences detected")
                    best_method = method_stats.index[0]
                    print(f"🏆 Best method (numerically): {best_method}")
            
            p_results[metric] = {
                'best_method': best_method if 'best_method' in locals() else method_stats.index[0],
                'method_stats': method_stats,
                'significant': p_val < 0.05 if 'p_val' in locals() else False,
                'p_value': p_val if 'p_val' in locals() else None
            }
        
        results_summary[p_value] = p_results
    
    return results_summary

def tradeoff_analysis(data, dataset_name):
    """
    Analyze tradeoffs between speed and quality metrics
    """
    
    print(f"\n{'='*80}")
    print(f"TRADEOFF ANALYSIS - {dataset_name.upper()} DATASET")
    print(f"{'='*80}")
    
    for p_value in sorted(data['P'].unique()):
        p_data = data[data['P'] == p_value]
        
        print(f"\n{'-'*60}")
        print(f"P Value = {p_value} - Speed vs Quality Tradeoffs")
        print(f"{'-'*60}")
        
        # Calculate method averages
        method_summary = p_data.groupby('method').agg({
            'init_time': 'mean',
            'total_time': 'mean',
            'post_ls_gap': 'mean'
        }).round(4)
        
        # Normalize metrics for comparison (0-1 scale)
        normalized = method_summary.copy()
        for col in normalized.columns:
            min_val, max_val = normalized[col].min(), normalized[col].max()
            if max_val != min_val:
                normalized[col] = (normalized[col] - min_val) / (max_val - min_val)
        
        # Calculate efficiency scores (lower time + lower gap = better)
        method_summary['speed_score'] = 1 - normalized['total_time']  # Higher = faster
        method_summary['quality_score'] = 1 - normalized['post_ls_gap']  # Higher = better quality
        method_summary['efficiency_score'] = (method_summary['speed_score'] + method_summary['quality_score']) / 2
        
        # Sort by efficiency
        method_summary = method_summary.sort_values('efficiency_score', ascending=False)
        
        print("\n🎯 METHOD EFFICIENCY RANKING:")
        print("=" * 50)
        print(f"{'Method':<20} {'Total Time':<12} {'Gap %':<10} {'Speed Score':<12} {'Quality Score':<13} {'Efficiency':<10}")
        print("-" * 85)
        
        for method in method_summary.index:
            row = method_summary.loc[method]
            print(f"{method:<20} {row['total_time']:<12.4f} {row['post_ls_gap']:<10.4f} "
                  f"{row['speed_score']:<12.4f} {row['quality_score']:<13.4f} {row['efficiency_score']:<10.4f}")
        
        # Identify specialized methods
        print(f"\n🚀 SPEED SPECIALIST: {method_summary['speed_score'].idxmax()}")
        print(f"   → Fastest total time: {method_summary.loc[method_summary['speed_score'].idxmax(), 'total_time']:.4f} seconds")
        
        print(f"\n🎯 QUALITY SPECIALIST: {method_summary['quality_score'].idxmax()}")
        print(f"   → Best optimality gap: {method_summary.loc[method_summary['quality_score'].idxmax(), 'post_ls_gap']:.4f}%")
        
        print(f"\n⚖️ BEST OVERALL: {method_summary['efficiency_score'].idxmax()}")
        print(f"   → Balanced performance with efficiency score: {method_summary['efficiency_score'].max():.4f}")
        
        # Tradeoff insights
        fastest_method = method_summary['total_time'].idxmin()
        best_quality_method = method_summary['post_ls_gap'].idxmin()
        
        if fastest_method != best_quality_method:
            fastest_time = method_summary.loc[fastest_method, 'total_time']
            fastest_gap = method_summary.loc[fastest_method, 'post_ls_gap']
            quality_time = method_summary.loc[best_quality_method, 'total_time']
            quality_gap = method_summary.loc[best_quality_method, 'post_ls_gap']
            
            time_tradeoff = ((quality_time - fastest_time) / fastest_time) * 100
            gap_tradeoff = ((fastest_gap - quality_gap) / quality_gap) * 100
            
            print(f"\n⚡ SPEED vs QUALITY TRADEOFF:")
            print(f"   → {fastest_method} is {time_tradeoff:.1f}% faster than {best_quality_method}")
            print(f"   → But {best_quality_method} has {gap_tradeoff:.1f}% better optimality gap")

def create_comprehensive_comparison_plots(data, dataset_name):
    """
    Create enhanced visualizations showing speed vs quality tradeoffs
    """
    
    fig = plt.figure(figsize=(18, 6))
    
    for i, p_value in enumerate(sorted(data['P'].unique())):
        p_data = data[data['P'] == p_value]
        
        # Speed vs Quality scatter plot
        plt.subplot(1, 3, i+1)
        method_means = p_data.groupby('method').agg({
            'total_time': 'mean',
            'post_ls_gap': 'mean'
        })
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(method_means)))
        for j, (method, row) in enumerate(method_means.iterrows()):
            plt.scatter(row['total_time'], row['post_ls_gap'], 
                       s=200, c=[colors[j]], alpha=0.7, label=method)
            plt.annotate(method, (row['total_time'], row['post_ls_gap']), 
                        xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        plt.xlabel('Total Time (seconds)')
        plt.ylabel('Optimality Gap (%)')
        plt.title(f'Speed vs Quality - P={p_value}')
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{dataset_name}_speed_vs_quality.png', dpi=300, bbox_inches='tight')
    plt.show()

def method_recommendation_engine(data, dataset_name):
    """
    Intelligent method recommendation based on use case priorities
    """
    
    print(f"\n{'='*80}")
    print(f"METHOD RECOMMENDATION ENGINE - {dataset_name.upper()} DATASET")
    print(f"{'='*80}")
    
    use_cases = {
        'speed_critical': {'total_time': 0.8, 'post_ls_gap': 0.2},
        'quality_critical': {'total_time': 0.2, 'post_ls_gap': 0.8},
        'balanced': {'total_time': 0.5, 'post_ls_gap': 0.5},
        'initialization_speed': {'init_time': 0.9, 'post_ls_gap': 0.1}
    }
    
    for p_value in sorted(data['P'].unique()):
        p_data = data[data['P'] == p_value]
        
        print(f"\n{'-'*60}")
        print(f"RECOMMENDATIONS FOR P = {p_value}")
        print(f"{'-'*60}")
        
        # Calculate method averages
        method_stats = p_data.groupby('method').agg({
            'init_time': 'mean',
            'total_time': 'mean',
            'post_ls_gap': 'mean'
        }).round(4)
        
        # Normalize for scoring (lower is better for all metrics)
        normalized = method_stats.copy()
        for col in normalized.columns:
            min_val, max_val = normalized[col].min(), normalized[col].max()
            if max_val != min_val:
                # Invert so that lower values get higher scores
                normalized[col] = 1 - (normalized[col] - min_val) / (max_val - min_val)
        
        print(f"{'Use Case':<20} {'Recommended Method':<20} {'Score':<8} {'Reasoning':<50}")
        print("-" * 100)
        
        for use_case, weights in use_cases.items():
            scores = {}
            for method in method_stats.index:
                score = 0
                for metric, weight in weights.items():
                    if metric in normalized.columns:
                        score += weight * normalized.loc[method, metric]
                scores[method] = score
            
            best_method = max(scores, key=scores.get)
            best_score = scores[best_method]
            
            # Generate reasoning
            reasoning = ""
            if use_case == 'speed_critical':
                time_val = method_stats.loc[best_method, 'total_time']
                reasoning = f"Fastest total time ({time_val:.3f}s)"
            elif use_case == 'quality_critical':
                gap_val = method_stats.loc[best_method, 'post_ls_gap']
                reasoning = f"Best optimality gap ({gap_val:.3f}%)"
            elif use_case == 'balanced':
                reasoning = f"Best speed-quality balance"
            elif use_case == 'initialization_speed':
                init_val = method_stats.loc[best_method, 'init_time']
                reasoning = f"Fastest initialization ({init_val:.3f}s)"
            
            print(f"{use_case:<20} {best_method:<20} {best_score:<8.3f} {reasoning:<50}")

# Updated main analysis function
def comprehensive_method_analysis(file_path, dataset_name):
    """
    Complete enhanced analysis pipeline
    """
    
    print(f"🚀 LOADING {dataset_name.upper()} DATASET...")
    df = pd.read_csv(file_path)
    
    # Filter for successful runs
    successful_data = df[df['status'] == 'success'].copy()
    successful_data['post_ls_gap'] = successful_data['post_ls_gap'] * 100  # Convert to percentage
    
    print(f"✅ Data loaded: {len(successful_data)} successful runs out of {len(df)} total")
    print(f"📊 Methods: {list(successful_data['method'].unique())}")
    print(f"🎯 P values: {sorted(successful_data['P'].unique())}")
    
    # Run enhanced analyses
    statistical_results = enhanced_statistical_analysis(successful_data, dataset_name)
    tradeoff_analysis(successful_data, dataset_name)
    create_comprehensive_comparison_plots(successful_data, dataset_name)
    method_recommendation_engine(successful_data, dataset_name)
    
    return successful_data, statistical_results

# Usage example
if __name__ == "__main__":
    # Replace with your actual file paths
    datasets = {
        'small': 'combined_results_small.csv',
        'medium': 'combined_results_medium.csv',
        'large': 'combined_results_large.csv'
    }
    
    all_results = {}
    
    for dataset_name, file_path in datasets.items():
        print(f"\n{'#'*100}")
        print(f"PROCESSING {dataset_name.upper()} DATASET")
        print(f"{'#'*100}")
        
        data, results = comprehensive_method_analysis(file_path, dataset_name)
        all_results[dataset_name] = {'data': data, 'results': results}
        
        print(f"\n✅ {dataset_name.upper()} ANALYSIS COMPLETE!")
    
    print(f"\n{'#'*100}")
    print("🎉 ALL ANALYSES COMPLETE!")
    print(f"{'#'*100}")