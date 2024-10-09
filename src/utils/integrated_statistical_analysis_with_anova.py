
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Assuming the data is loaded in a dataframe called 'df'
mean_runtimes = df.groupby('threads')['milliseconds'].mean()

# Polynomial Regression (Degree 2)
coeffs_poly = np.polyfit(mean_runtimes.index, mean_runtimes.values, 2)
poly_eq = np.poly1d(coeffs_poly)

# Logarithmic Fit: y = a * log(x) + b
log_coeffs = np.polyfit(np.log(mean_runtimes.index), mean_runtimes.values, 1)
log_eq = lambda x: log_coeffs[0] * np.log(x) + log_coeffs[1]

# Exponential Decay: y = a * exp(-b * x) + c
def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

# Power Law: y = a * x^(-b)
def power_law(x, a, b):
    return a * np.power(x, -b)

# Convert thread counts to float to avoid issues with negative powers in power law
thread_counts_float = mean_runtimes.index.astype(float)

# Fitting models
popt_exp, _ = curve_fit(exp_decay, thread_counts_float, mean_runtimes.values, p0=(40000, 0.1, 1000), maxfev=2000)
popt_power, _ = curve_fit(power_law, thread_counts_float, mean_runtimes.values, p0=(40000, 1), maxfev=2000)

# Error calculation (RSS)
def rss(y_true, y_pred):
    return np.sum((y_true - y_pred) ** 2)

rss_poly = rss(mean_runtimes.values, poly_eq(mean_runtimes.index))
rss_log = rss(mean_runtimes.values, log_eq(mean_runtimes.index))
rss_exp = rss(mean_runtimes.values, exp_decay(thread_counts_float, *popt_exp))
rss_power = rss(mean_runtimes.values, power_law(thread_counts_float, *popt_power))

# Plotting the curves with their respective errors
plt.figure(figsize=(12, 8))
plt.scatter(thread_counts_float, mean_runtimes.values, label="Mean Runtimes", color='blue')

# Polynomial Fit
plt.plot(mean_runtimes.index, poly_eq(mean_runtimes.index), label=f"Polynomial Fit (Degree 2), RSS: {rss_poly:.2f}", color='purple')

# Logarithmic Fit
plt.plot(mean_runtimes.index, log_eq(mean_runtimes.index), label=f"Logarithmic Fit, RSS: {rss_log:.2f}", color='orange')

# Exponential Decay Fit
plt.plot(thread_counts_float, exp_decay(thread_counts_float, *popt_exp), label=f"Exponential Decay Fit, RSS: {rss_exp:.2f}", color='red')

# Power Law Fit
plt.plot(thread_counts_float, power_law(thread_counts_float, *popt_power), label=f"Power Law Fit, RSS: {rss_power:.2f}", color='green')

plt.title("Mean Runtime vs Thread Count with Various Curve Fits and RSS")
plt.xlabel("Number of Threads")
plt.ylabel("Mean Runtime (ms)")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Saving the plot
plt.savefig("runtime_curve_fits_with_errors.png")
plt.show()

# Print RSS values for reference
print("RSS values for each fit:")
print(f"Polynomial Fit (Degree 2): {rss_poly}")
print(f"Logarithmic Fit: {rss_log}")
print(f"Exponential Decay Fit: {rss_exp}")
print(f"Power Law Fit: {rss_power}")

# One-way ANOVA on Objective Values ('value')
anova_model = ols('value ~ C(threads)', data=df).fit()
anova_table = sm.stats.anova_lm(anova_model, typ=2)

# Print ANOVA results
print("\nOne-way ANOVA on Objective Values (value):")
print(anova_table)
