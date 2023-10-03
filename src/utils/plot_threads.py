import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data
df = pd.read_csv('output.csv')

# Use the stylish theme
plt.style.use('ggplot')

# Define distinct colors for differentiation using the Set3 palette
colors = plt.cm.Set3(np.linspace(0, 1, len(df['instance'].unique())))

# Set up the figure with high resolution
plt.figure(figsize=(15, 10), dpi=300)

# Plot each instance's data
for color, (instance, group_data) in zip(colors, df.groupby('instance')):
    plt.plot(group_data[' threads'], group_data[' milliseconds'], 
             label=instance, color=color, linewidth=1.0, marker='o', markersize=4)

# Customize the plot
plt.xlabel('Number of Threads')
plt.ylabel('Runtime (in milliseconds)')
plt.title('Runtime vs Number of Threads for Different Instances')
plt.legend(loc="upper right", bbox_to_anchor=(1.3, 1), fontsize='small')
plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.tight_layout()

# Save the plot as a PDF
plt.savefig("output_plot.pdf", format="pdf")

# Display the plot (optional if you just want the PDF)
plt.show()
