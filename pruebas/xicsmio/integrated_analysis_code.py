
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import f_oneway, stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

# Load data
data_small = pd.read_csv("/path_to_small_dataset.csv")
data_medium = pd.read_csv("/path_to_medium_dataset.csv")

# Descriptive Statistics
def descriptive_statistics(data):
    return data.groupby(['alpha_a', 'alpha_l', 'iters']).agg(['mean', 'std'])

# One-way ANOVA
def one_way_anova(data):
    return (
        f_oneway(*[data['ValueGRASP'][data['iters'] == i] for i in data['iters'].unique()]),
        f_oneway(*[data['TotalTimeGRASP'][data['iters'] == i] for i in data['iters'].unique()])
    )

# Factorial ANOVA
def factorial_anova(data):
    model_value = ols('ValueGRASP ~ C(alpha_a) * C(alpha_l) * C(iters)', data).fit()
    model_time = ols('TotalTimeGRASP ~ C(alpha_a) * C(alpha_l) * C(iters)', data).fit()
    return (
        anova_lm(model_value, typ=2),
        anova_lm(model_time, typ=2)
    )

# Tukey HSD
def tukey_hsd(data):
    return (
        pairwise_tukeyhsd(data['ValueGRASP'], data['iters']).summary(),
        pairwise_tukeyhsd(data['TotalTimeGRASP'], data['iters']).summary()
    )

# Tukey-Kramer Test
def tukey_kramer(data, variable):
    unique_values = sorted(data['iters'].unique())
    value_combinations = [(i, j) for idx, i in enumerate(unique_values) for j in unique_values[idx+1:]]
    rejected_list = []

    for combo in value_combinations:
        group1 = data[data['iters'] == combo[0]][variable]
        group2 = data[data['iters'] == combo[1]][variable]
        t_stat, p_val = stats.ttest_ind(group1, group2)
        
        if p_val < 0.05:  # Using 0.05 as the significance level
            rejected_list.append((combo, True, p_val))
        else:
            rejected_list.append((combo, False, p_val))

    return pd.DataFrame(rejected_list, columns=['Group Combo', 'Reject Null?', 'p-value'])

# Main Effects Plots
def plot_main_effects(dataset, title_suffix):
    grouped_means = dataset.groupby(['alpha_a', 'alpha_l', 'iters']).agg({'ValueGRASP': 'mean'}).reset_index()
    fig, ax = plt.subplots(3, 1, figsize=(15, 20), dpi=150)
    grouped_means.groupby('alpha_a')['ValueGRASP'].mean().plot(kind='line', marker='o', ax=ax[0])
    grouped_means.groupby('alpha_l')['ValueGRASP'].mean().plot(kind='line', marker='o', ax=ax[1])
    grouped_means.groupby('iters')['ValueGRASP'].mean().plot(kind='line', marker='o', ax=ax[2])
    plt.tight_layout()
    plt.show()

# Pareto Front Plot
def plot_pareto_front(dataset, title_suffix):
    pareto_data = dataset.groupby(['iters']).agg({
        'TotalTimeGRASP': 'mean',
        'ValueGRASP': 'mean'
    }).reset_index()
    pareto_data_sorted = pareto_data.sort_values(by='ValueGRASP')
    pareto_front = [pareto_data_sorted.iloc[0]]
    for i in range(1, len(pareto_data_sorted)):
        if pareto_data_sorted.iloc[i]['TotalTimeGRASP'] < pareto_front[-1]['TotalTimeGRASP']:
            pareto_front.append(pareto_data_sorted.iloc[i])
    plt.figure(figsize=(12, 8), dpi=150)
    plt.scatter(pareto_data['TotalTimeGRASP'], pareto_data['ValueGRASP'], c='blue', label='Data Points')
    plt.plot([point['TotalTimeGRASP'] for point in pareto_front], [point['ValueGRASP'] for point in pareto_front], c='red', label='Pareto Front')
    plt.xlabel('Mean Runtime (TotalTimeGRASP)')
    plt.ylabel('Mean Objective Value (ValueGRASP)')
    plt.title(f'Pareto Front for Runtime vs. Objective Value {title_suffix}')
    plt.legend()
    plt.grid(True)
    plt.show()

# Note: After loading the datasets, you can call the above functions as required to perform the respective analyses and plots.
