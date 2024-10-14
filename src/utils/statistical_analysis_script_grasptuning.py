
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# Load the dataset
df_625 = pd.read_csv('df_grasp_625_allcombs_1to10.csv')

# One-Way ANOVA for runtime (TotalTimeGRASP) and objective function (ValueGRASP)
iters_groups_runtime = [group['TotalTimeGRASP'].values for name, group in df_625.groupby('iters')]
alpha_a_groups_runtime = [group['TotalTimeGRASP'].values for name, group in df_625.groupby('alpha_a')]
alpha_l_groups_runtime = [group['TotalTimeGRASP'].values for name, group in df_625.groupby('alpha_l')]

iters_groups_obj = [group['ValueGRASP'].values for name, group in df_625.groupby('iters')]
alpha_a_groups_obj = [group['ValueGRASP'].values for name, group in df_625.groupby('alpha_a')]
alpha_l_groups_obj = [group['ValueGRASP'].values for name, group in df_625.groupby('alpha_l')]

anova_runtime_iters = stats.f_oneway(*iters_groups_runtime)
anova_runtime_alpha_a = stats.f_oneway(*alpha_a_groups_runtime)
anova_runtime_alpha_l = stats.f_oneway(*alpha_l_groups_runtime)

anova_obj_iters = stats.f_oneway(*iters_groups_obj)
anova_obj_alpha_a = stats.f_oneway(*alpha_a_groups_obj)
anova_obj_alpha_l = stats.f_oneway(*alpha_l_groups_obj)

# Factorial ANOVA for interactions
formula_runtime = 'TotalTimeGRASP ~ C(iters) + C(alpha_a) + C(alpha_l) + C(iters):C(alpha_a) + C(iters):C(alpha_l) + C(alpha_a):C(alpha_l)'
model_runtime = ols(formula_runtime, data=df_625).fit()
anova_table_runtime = sm.stats.anova_lm(model_runtime, typ=2)

formula_obj = 'ValueGRASP ~ C(iters) + C(alpha_a) + C(alpha_l) + C(iters):C(alpha_a) + C(iters):C(alpha_l) + C(alpha_a):C(alpha_l)'
model_obj = ols(formula_obj, data=df_625).fit()
anova_table_obj = sm.stats.anova_lm(model_obj, typ=2)

# Tukey HSD tests
tukey_runtime_iters = pairwise_tukeyhsd(endog=df_625['TotalTimeGRASP'], groups=df_625['iters'], alpha=0.05)
tukey_obj_iters = pairwise_tukeyhsd(endog=df_625['ValueGRASP'], groups=df_625['iters'], alpha=0.05)

tukey_runtime_alpha_a = pairwise_tukeyhsd(endog=df_625['TotalTimeGRASP'], groups=df_625['alpha_a'], alpha=0.05)
tukey_obj_alpha_a = pairwise_tukeyhsd(endog=df_625['ValueGRASP'], groups=df_625['alpha_a'], alpha=0.05)

tukey_runtime_alpha_l = pairwise_tukeyhsd(endog=df_625['TotalTimeGRASP'], groups=df_625['alpha_l'], alpha=0.05)
tukey_obj_alpha_l = pairwise_tukeyhsd(endog=df_625['ValueGRASP'], groups=df_625['alpha_l'], alpha=0.05)

# Aggregations for min, mean, max, and std dev for alpha_a, alpha_l, and iters
agg_functions = {
    'TotalTimeGRASP': ['min', 'mean', 'max', 'std'],
    'ValueGRASP': ['min', 'mean', 'max', 'std'],
    'Gap%': ['min', 'mean', 'max', 'std']
}

alpha_a_stats = df_625.groupby('alpha_a').agg(agg_functions)
alpha_l_stats = df_625.groupby('alpha_l').agg(agg_functions)
iters_stats = df_625.groupby('iters').agg(agg_functions)

# Display results
print('ANOVA Results (Runtime - Iters):', anova_runtime_iters)
print('ANOVA Results (Runtime - Alpha_a):', anova_runtime_alpha_a)
print('ANOVA Results (Runtime - Alpha_l):', anova_runtime_alpha_l)

print('ANOVA Results (Objective - Iters):', anova_obj_iters)
print('ANOVA Results (Objective - Alpha_a):', anova_obj_alpha_a)
print('ANOVA Results (Objective - Alpha_l):', anova_obj_alpha_l)

print('Factorial ANOVA (Runtime):\n', anova_table_runtime)
print('Factorial ANOVA (Objective):\n', anova_table_obj)

print(tukey_runtime_iters)
print(tukey_obj_iters)
print(tukey_runtime_alpha_a)
print(tukey_obj_alpha_a)
print(tukey_runtime_alpha_l)
print(tukey_obj_alpha_l)

print('Alpha_a Stats:\n', alpha_a_stats)
print('Alpha_l Stats:\n', alpha_l_stats)
print('Iters Stats:\n', iters_stats)
