
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import ttest_rel

# Load the data
df = pd.read_csv('path_to_your_output.csv')

# Descriptive Statistics
descriptive_stats = df.groupby(' threads')[' milliseconds'].agg(['mean', 'std'])
print(descriptive_stats)

# Repeated Measures ANOVA
df_long = df.pivot(index='instance', columns=' threads', values=' milliseconds').reset_index().melt(id_vars=['instance'], 
                                                                                                value_name='Runtime', 
                                                                                                var_name='Threads')
anova = AnovaRM(df_long, 'Runtime', 'instance', within=['Threads'])
anova_results = anova.fit()
print(anova_results.summary())

# Tukey's HSD Post-hoc Test
posthoc = pairwise_tukeyhsd(df_long['Runtime'], df_long['Threads'], alpha=0.05)
print(posthoc.summary())

# Effect Size
def cohen_d(group1, group2):
    n1, n2 = len(group1), len(group2)
    s1, s2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    return (np.mean(group1) - np.mean(group2)) / pooled_std

effect_sizes = {}
for thread in df[' threads'].unique():
    if thread != 1:
        d = cohen_d(df[df[' threads'] == 1][' milliseconds'], df[df[' threads'] == thread][' milliseconds'])
        effect_sizes[thread] = d
print("\nEffect Sizes:", effect_sizes)

# Paired T-Tests
t_test_results = {}
thread_counts = sorted(df[' threads'].unique())
for i in range(len(thread_counts) - 1):
    data1 = df[df[' threads'] == thread_counts[i]][' milliseconds']
    data2 = df[df[' threads'] == thread_counts[i + 1]][' milliseconds']
    t_stat, p_val = ttest_rel(data1, data2)
    t_test_results[(thread_counts[i], thread_counts[i + 1])] = (t_stat, p_val)
print("\nPaired T-Tests Results:", t_test_results)

# Box Plots
plt.figure(figsize=(12, 8), dpi=150)
boxplot = df.boxplot(column=' milliseconds', by=' threads', grid=False, patch_artist=True, showfliers=False)
plt.title('Distribution of Runtimes for Different Thread Counts')
plt.suptitle('')  # Suppress the default title
plt.ylabel('Runtime (in milliseconds)')
plt.xlabel('Number of Threads')
plt.tight_layout()
plt.savefig("boxplots_output.png")
plt.show()
