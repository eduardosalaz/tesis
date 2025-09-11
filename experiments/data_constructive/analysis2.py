import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

class HeuristicMethodAnalyzer:
    """
    Advanced statistical analyzer for heuristic initialization methods
    """
    
    def __init__(self, file_path, dataset_name):
        self.file_path = file_path
        self.dataset_name = dataset_name
        self.data = None
        self.aggregated_data = None
        self.results = {}
        
    def load_data(self):
        """Load and preprocess the data"""
        print(f"\nLoading {self.dataset_name} dataset...")
        self.data = pd.read_csv(self.file_path)
        
        # Filter successful runs
        self.data = self.data[self.data['status'] == 'success'].copy()
        
        # Add final gap if not present
        if 'final_gap' not in self.data.columns:
            self.data['final_gap'] = self.data['post_ls_gap'] * 100
        
        print(f"Loaded {len(self.data)} successful runs")
        print(f"Methods: {self.data['method'].unique()}")
        print(f"P values: {sorted(self.data['P'].unique())}")
        
        # Aggregate by instance (average over 3 runs)
        self._aggregate_instances()
        
    def _aggregate_instances(self):
        """Aggregate data by instance (average over 3 runs)"""
        agg_columns = ['init_time', 'prop_split_bus', 'final_gap', 'total_time', 'post_ls_obj']
        
        self.aggregated_data = self.data.groupby(
            ['instance_id_x', 'P', 'method']
        )[agg_columns].mean().reset_index()
        
        print(f"Aggregated to {len(self.aggregated_data)} instance-method combinations")
        
    def compute_efficiency_metrics(self):
        """
        Compute efficiency metrics combining time and quality
        """
        print(f"\n{'='*60}")
        print(f"EFFICIENCY METRICS - {self.dataset_name.upper()}")
        print(f"{'='*60}")
        
        for p_value in sorted(self.aggregated_data['P'].unique()):
            print(f"\n--- P = {p_value} ---")
            
            p_data = self.aggregated_data[self.aggregated_data['P'] == p_value].copy()
            
            # Compute efficiency scores
            # Efficiency = Quality / Time (higher is better)
            # Quality = 1 / (1 + gap) to avoid division by zero
            p_data['quality_score'] = 1 / (1 + p_data['final_gap'])
            p_data['time_efficiency'] = p_data['quality_score'] / p_data['init_time']
            p_data['overall_efficiency'] = p_data['quality_score'] / p_data['total_time']
            
            # Store results
            self.results[f'efficiency_P{p_value}'] = p_data
            
            # Print summary
            efficiency_summary = p_data.groupby('method')[
                ['time_efficiency', 'overall_efficiency']
            ].agg(['mean', 'std'])
            
            print("\nEfficiency Metrics (higher is better):")
            print(efficiency_summary.round(4))
            
            # Rank methods
            rankings = {}
            for metric in ['time_efficiency', 'overall_efficiency']:
                ranked = p_data.groupby('method')[metric].mean().sort_values(ascending=False)
                rankings[metric] = list(ranked.index)
            
            print("\nMethod Rankings by Efficiency:")
            for metric, ranking in rankings.items():
                print(f"  {metric}: {' > '.join(ranking)}")
                
    def analyze_consistency(self):
        """
        Analyze consistency/reliability of methods across instances
        """
        print(f"\n{'='*60}")
        print(f"CONSISTENCY ANALYSIS - {self.dataset_name.upper()}")
        print(f"{'='*60}")
        
        consistency_results = {}
        
        for p_value in sorted(self.aggregated_data['P'].unique()):
            print(f"\n--- P = {p_value} ---")
            
            p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
            
            consistency_metrics = []
            
            for method in p_data['method'].unique():
                method_data = p_data[p_data['method'] == method]
                
                # Calculate coefficient of variation (CV) for each metric
                # Lower CV indicates more consistent performance
                metrics_cv = {}
                for metric in ['init_time', 'final_gap', 'prop_split_bus']:
                    mean_val = method_data[metric].mean()
                    std_val = method_data[metric].std()
                    cv = (std_val / mean_val * 100) if mean_val != 0 else 0
                    metrics_cv[f'{metric}_cv'] = cv
                
                # Calculate inter-quartile range ratio (IQR/median)
                metrics_iqr = {}
                for metric in ['init_time', 'final_gap', 'prop_split_bus']:
                    q1 = method_data[metric].quantile(0.25)
                    q3 = method_data[metric].quantile(0.75)
                    median = method_data[metric].median()
                    iqr_ratio = ((q3 - q1) / median * 100) if median != 0 else 0
                    metrics_iqr[f'{metric}_iqr_ratio'] = iqr_ratio
                
                consistency_metrics.append({
                    'method': method,
                    **metrics_cv,
                    **metrics_iqr,
                    'overall_consistency': np.mean(list(metrics_cv.values()))
                })
            
            consistency_df = pd.DataFrame(consistency_metrics)
            consistency_df = consistency_df.sort_values('overall_consistency')
            
            consistency_results[p_value] = consistency_df
            
            print("\nConsistency Metrics (lower values = more consistent):")
            print(consistency_df.set_index('method').round(2))
            
            print(f"\nMost consistent method: {consistency_df.iloc[0]['method']}")
            print(f"Least consistent method: {consistency_df.iloc[-1]['method']}")
        
        self.results['consistency'] = consistency_results
        
    def perform_dominance_analysis(self):
        """
        Perform pairwise dominance analysis
        """
        print(f"\n{'='*60}")
        print(f"DOMINANCE ANALYSIS - {self.dataset_name.upper()}")
        print(f"{'='*60}")
        
        for p_value in sorted(self.aggregated_data['P'].unique()):
            print(f"\n--- P = {p_value} ---")
            
            p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
            methods = p_data['method'].unique()
            
            # Create dominance matrix
            dominance_matrix = pd.DataFrame(
                index=methods, 
                columns=methods,
                data=0
            )
            
            # For each pair of methods, count instances where one dominates
            instances = p_data['instance_id_x'].unique()
            
            for instance in instances:
                inst_data = p_data[p_data['instance_id_x'] == instance]
                
                for method1 in methods:
                    m1_data = inst_data[inst_data['method'] == method1]
                    if len(m1_data) == 0:
                        continue
                        
                    for method2 in methods:
                        if method1 == method2:
                            continue
                            
                        m2_data = inst_data[inst_data['method'] == method2]
                        if len(m2_data) == 0:
                            continue
                        
                        # Check if method1 dominates method2 on this instance
                        # (better gap AND better time)
                        if (m1_data['final_gap'].values[0] < m2_data['final_gap'].values[0] and
                            m1_data['init_time'].values[0] < m2_data['init_time'].values[0]):
                            dominance_matrix.loc[method1, method2] += 1
            
            # Calculate dominance percentage
            n_instances = len(instances)
            dominance_pct = (dominance_matrix / n_instances * 100).round(1)
            
            print("\nDominance Matrix (% of instances where row dominates column):")
            print(dominance_pct)
            
            # Calculate dominance scores
            dominance_scores = dominance_pct.sum(axis=1) - dominance_pct.sum(axis=0)
            dominance_scores = dominance_scores.sort_values(ascending=False)
            
            print("\nNet Dominance Scores (higher = more dominant):")
            for method, score in dominance_scores.items():
                print(f"  {method}: {score:.1f}")
                
    def analyze_performance_profiles(self):
        """
        Create performance profiles (Dolan-Moré profiles)
        """
        print(f"\n{'='*60}")
        print(f"PERFORMANCE PROFILES - {self.dataset_name.upper()}")
        print(f"{'='*60}")
        
        for p_value in sorted(self.aggregated_data['P'].unique()):
            p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
            
            # Create performance profiles for different metrics
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            
            metrics = ['init_time', 'final_gap', 'total_time']
            
            for idx, metric in enumerate(metrics):
                ax = axes[idx]
                
                # Pivot data for easier computation
                pivot = p_data.pivot_table(
                    index='instance_id_x', 
                    columns='method', 
                    values=metric
                )
                
                # Calculate performance ratios
                min_per_instance = pivot.min(axis=1)
                
                # Plot performance profile for each method
                for method in pivot.columns:
                    ratios = pivot[method] / min_per_instance
                    ratios = ratios.dropna().sort_values()
                    
                    # Calculate cumulative distribution
                    y = np.arange(1, len(ratios) + 1) / len(ratios)
                    
                    ax.plot(ratios, y, label=method, linewidth=2)
                
                ax.set_xlabel('Performance Ratio τ')
                ax.set_ylabel('P(r ≤ τ)')
                ax.set_title(f'{metric.replace("_", " ").title()}')
                ax.legend()
                ax.grid(True, alpha=0.3)
                ax.set_xlim(1, None)
            
            plt.suptitle(f'Performance Profiles - {self.dataset_name.upper()} (P={p_value})', 
                        fontsize=14)
            plt.tight_layout()
            plt.savefig(f'{self.dataset_name}_P{p_value}_performance_profiles.png', 
                       dpi=300, bbox_inches='tight')
            plt.show()
            
    def analyze_correlations(self):
        """
        Analyze correlations between different metrics
        """
        print(f"\n{'='*60}")
        print(f"CORRELATION ANALYSIS - {self.dataset_name.upper()}")
        print(f"{'='*60}")
        
        for p_value in sorted(self.aggregated_data['P'].unique()):
            print(f"\n--- P = {p_value} ---")
            
            p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
            
            # Calculate correlations for each method
            for method in p_data['method'].unique():
                method_data = p_data[p_data['method'] == method]
                
                if len(method_data) > 3:  # Need enough data for correlation
                    corr_matrix = method_data[
                        ['init_time', 'prop_split_bus', 'final_gap', 'total_time']
                    ].corr()
                    
                    print(f"\n{method} Correlations:")
                    
                    # Find significant correlations
                    for i in range(len(corr_matrix.columns)):
                        for j in range(i+1, len(corr_matrix.columns)):
                            col1 = corr_matrix.columns[i]
                            col2 = corr_matrix.columns[j]
                            corr_val = corr_matrix.iloc[i, j]
                            
                            if abs(corr_val) > 0.5:  # Strong correlation
                                print(f"  {col1} vs {col2}: {corr_val:.3f}")
                                
    def multivariate_analysis(self):
        """
        Perform multivariate analysis using PCA
        """
        print(f"\n{'='*60}")
        print(f"MULTIVARIATE ANALYSIS (PCA) - {self.dataset_name.upper()}")
        print(f"{'='*60}")
        
        for p_value in sorted(self.aggregated_data['P'].unique()):
            print(f"\n--- P = {p_value} ---")
            
            p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
            
            # Prepare data for PCA
            features = ['init_time', 'prop_split_bus', 'final_gap', 'total_time']
            X = p_data[features].values
            
            # Standardize features
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            
            # Perform PCA
            pca = PCA()
            X_pca = pca.fit_transform(X_scaled)
            
            # Print explained variance
            print("\nPCA Explained Variance Ratio:")
            for i, var_exp in enumerate(pca.explained_variance_ratio_):
                print(f"  PC{i+1}: {var_exp:.3f} ({var_exp*100:.1f}%)")
            
            print(f"\nCumulative variance explained by first 2 PCs: "
                  f"{pca.explained_variance_ratio_[:2].sum()*100:.1f}%")
            
            # Feature loadings
            loadings = pd.DataFrame(
                pca.components_[:2].T,
                columns=['PC1', 'PC2'],
                index=features
            )
            
            print("\nFeature Loadings on Principal Components:")
            print(loadings.round(3))
            
            # Visualize PCA
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
            
            # Scatter plot in PC space
            colors = {method: i for i, method in enumerate(p_data['method'].unique())}
            
            for method in p_data['method'].unique():
                mask = p_data['method'] == method
                ax1.scatter(X_pca[mask, 0], X_pca[mask, 1], 
                          label=method, alpha=0.6, s=50)
            
            ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
            ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
            ax1.set_title(f'PCA Scatter - {self.dataset_name.upper()} (P={p_value})')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Loading plot
            for i, feature in enumerate(features):
                ax2.arrow(0, 0, loadings.iloc[i, 0], loadings.iloc[i, 1],
                         head_width=0.05, head_length=0.05, fc='blue', ec='blue')
                ax2.text(loadings.iloc[i, 0]*1.1, loadings.iloc[i, 1]*1.1, 
                        feature, fontsize=10)
            
            ax2.set_xlabel('PC1')
            ax2.set_ylabel('PC2')
            ax2.set_title('PCA Loading Plot')
            ax2.set_xlim(-1, 1)
            ax2.set_ylim(-1, 1)
            ax2.grid(True, alpha=0.3)
            ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
            ax2.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
            
            plt.tight_layout()
            plt.savefig(f'{self.dataset_name}_P{p_value}_pca_analysis.png', 
                       dpi=300, bbox_inches='tight')
            plt.show()
            
    def bootstrap_confidence_intervals(self, n_bootstrap=1000):
        """
        Compute bootstrap confidence intervals for method performance
        """
        print(f"\n{'='*60}")
        print(f"BOOTSTRAP CONFIDENCE INTERVALS - {self.dataset_name.upper()}")
        print(f"{'='*60}")
        
        np.random.seed(42)  # For reproducibility
        
        for p_value in sorted(self.aggregated_data['P'].unique()):
            print(f"\n--- P = {p_value} ---")
            
            p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
            
            bootstrap_results = {}
            
            for method in p_data['method'].unique():
                method_data = p_data[p_data['method'] == method]
                
                bootstrap_results[method] = {}
                
                for metric in ['init_time', 'final_gap']:
                    values = method_data[metric].values
                    
                    # Perform bootstrap
                    bootstrap_means = []
                    for _ in range(n_bootstrap):
                        sample = np.random.choice(values, size=len(values), replace=True)
                        bootstrap_means.append(np.mean(sample))
                    
                    # Calculate confidence intervals
                    ci_lower = np.percentile(bootstrap_means, 2.5)
                    ci_upper = np.percentile(bootstrap_means, 97.5)
                    mean_estimate = np.mean(bootstrap_means)
                    
                    bootstrap_results[method][metric] = {
                        'mean': mean_estimate,
                        'ci_lower': ci_lower,
                        'ci_upper': ci_upper
                    }
            
            # Display results
            print("\n95% Bootstrap Confidence Intervals:")
            
            for metric in ['init_time', 'final_gap']:
                print(f"\n{metric.replace('_', ' ').title()}:")
                for method in sorted(bootstrap_results.keys()):
                    res = bootstrap_results[method][metric]
                    print(f"  {method}: {res['mean']:.4f} [{res['ci_lower']:.4f}, {res['ci_upper']:.4f}]")
                    
            self.results[f'bootstrap_P{p_value}'] = bootstrap_results
            
    def generate_comprehensive_report(self):
        """
        Generate a comprehensive report with all analyses
        """
        report_file = f'{self.dataset_name}_comprehensive_statistical_report.txt'
        
        with open(report_file, 'w') as f:
            f.write(f"{'='*80}\n")
            f.write(f"COMPREHENSIVE STATISTICAL ANALYSIS REPORT\n")
            f.write(f"Dataset: {self.dataset_name.upper()}\n")
            f.write(f"{'='*80}\n\n")
            
            # Summary statistics
            f.write("SUMMARY STATISTICS\n")
            f.write("-" * 40 + "\n")
            for p_value in sorted(self.aggregated_data['P'].unique()):
                p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
                
                f.write(f"\nP = {p_value}:\n")
                summary = p_data.groupby('method')[
                    ['init_time', 'final_gap', 'prop_split_bus']
                ].agg(['mean', 'std', 'min', 'max'])
                
                f.write(summary.to_string())
                f.write("\n")
            
            # Key findings
            f.write("\n" + "="*80 + "\n")
            f.write("KEY FINDINGS AND RECOMMENDATIONS\n")
            f.write("-" * 40 + "\n")
            
            for p_value in sorted(self.aggregated_data['P'].unique()):
                f.write(f"\nP = {p_value}:\n")
                
                p_data = self.aggregated_data[self.aggregated_data['P'] == p_value]
                
                # Find best methods for different criteria
                best_time = p_data.groupby('method')['init_time'].mean().idxmin()
                best_gap = p_data.groupby('method')['final_gap'].mean().idxmin()
                
                f.write(f"  Fastest initialization: {best_time}\n")
                f.write(f"  Best solution quality: {best_gap}\n")
                
                # Check if same method is best for both
                if best_time == best_gap:
                    f.write(f"  *** {best_time} dominates on both criteria ***\n")
                else:
                    f.write(f"  Trade-off exists between {best_time} (speed) and {best_gap} (quality)\n")
            
            f.write("\n" + "="*80 + "\n")
            f.write("Report generated successfully\n")
        
        print(f"\nComprehensive report saved to {report_file}")
        
    def run_all_analyses(self):
        """
        Run all statistical analyses
        """
        self.load_data()
        self.compute_efficiency_metrics()
        self.analyze_consistency()
        self.perform_dominance_analysis()
        self.analyze_performance_profiles()
        self.analyze_correlations()
        self.multivariate_analysis()
        self.bootstrap_confidence_intervals()
        self.generate_comprehensive_report()
        
        return self.results

# Main execution
if __name__ == "__main__":
    # Analyze each dataset
    datasets = [
        ('combined_results_small.csv', 'small'),
        ('combined_results_medium.csv', 'medium'),
        ('combined_results_larg.csv', 'large')
    ]
    
    all_analyzers = {}
    
    for file_path, dataset_name in datasets:
        try:
            print(f"\n{'='*80}")
            print(f"ADVANCED ANALYSIS FOR {dataset_name.upper()} DATASET")
            print(f"{'='*80}")
            
            analyzer = HeuristicMethodAnalyzer(file_path, dataset_name)
            results = analyzer.run_all_analyses()
            all_analyzers[dataset_name] = analyzer
            
        except FileNotFoundError:
            print(f"Warning: {file_path} not found. Skipping.")
        except Exception as e:
            print(f"Error analyzing {dataset_name}: {e}")
    
    print(f"\n{'='*80}")
    print("ALL ANALYSES COMPLETE")
    print(f"{'='*80}")