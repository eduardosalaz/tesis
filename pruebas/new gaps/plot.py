import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams
import matplotlib as mpl

# Add Computer Modern font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman'] + plt.rcParams['font.serif']
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# Read the CSV file
df = pd.read_csv('output_grasp_625_threads.csv')

# Create the figure and axis objects
fig, ax = plt.subplots(figsize=(12, 8))

# Set a clean style
plt.style.use('seaborn-v0_8-whitegrid')

# Create the boxplot
boxplot = ax.boxplot([df[df['threads'] == i]['milliseconds'] for i in range(1, 9)],
                     patch_artist=True,
                     boxprops=dict(facecolor='#50C878', color='black'),
                     medianprops=dict(color='black'),
                     whiskerprops=dict(color='black'),
                     capprops=dict(color='black'),
                     flierprops=dict(color='black', markeredgecolor='black'))

# Customize the plot
ax.set_title(r'\textbf{Runtimes of the GRASP Metaheuristic with threads for Size 1 instances}', fontsize=16, pad=20)
ax.set_xlabel(r'\textbf{Number of threads}', fontsize=14)
ax.set_ylabel(r'\textbf{Runtime (milliseconds)}', fontsize=14)

# Increase font size for x-axis tick labels
ax.set_xticks(range(1, 9))
ax.set_xticklabels(range(1, 9), fontsize=12)

# Increase font size for y-axis tick labels
ax.tick_params(axis='y', labelsize=12)

# Adjust the y-axis to start from a lower value
y_min = df['milliseconds'].min() * 0.9  # Start at 90% of the minimum value
y_max = df['milliseconds'].max() * 1.1  # End at 110% of the maximum value
ax.set_ylim(bottom=y_min, top=y_max)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Adjust layout
plt.tight_layout()

# Save the plot in different formats with higher DPI
#plt.savefig('runtime_vs_threads_boxplot_larger_ticks_size2.svg', format='svg', dpi=500, bbox_inches='tight')
#plt.savefig('runtime_vs_threads_boxplot_larger_ticks_size2.pdf', format='pdf', dpi=500, bbox_inches='tight')
plt.savefig('runtime_vs_threads_boxplot_625.eps', format='eps', dpi=500, bbox_inches='tight')

# Display the plot (optional, comment out if not needed)
plt.show()