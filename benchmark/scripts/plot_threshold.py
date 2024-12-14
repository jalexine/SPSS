import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import warnings

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pandas")
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")

# Get the file path and optional output file name from the command-line arguments
if len(sys.argv) < 2:
    print("Usage: python script.py <file_path> [output_file]")
    sys.exit(1)

file_path = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) > 2 else file_path.replace('.csv', '_plot.png')

# Extract dataset name for the title
dataset_name = os.path.basename(file_path).replace(".csv", "")
dataset_name = dataset_name.replace("stats_", "").rsplit("_", 1)[0]

# Read the CSV data
data = pd.read_csv(file_path)


# Use column indices instead of names to match the expected structure
threshold_col = data.columns[2]
SPSS_Size_col = data.columns[5]
k_col = data.columns[1]

# Select only the necessary columns
data = data[[threshold_col, SPSS_Size_col, k_col]].copy()
data.columns = ["threshold", "SPSS_Size", "k"]

# Filter for k = 21, 31, 41
filtered_data = data[data['k'].isin([21, 31, 41])].copy()

# Convert columns to appropriate types
filtered_data.loc[:, "threshold"] = pd.to_numeric(filtered_data["threshold"], errors='coerce')
filtered_data.loc[:, "SPSS_Size"] = pd.to_numeric(filtered_data["SPSS_Size"], errors='coerce')

# Round SPSS_Size values for annotations
filtered_data.loc[:, "SPSS_Size"] = filtered_data["SPSS_Size"].round(-3)

# Plot the graph with curves for each k
plt.figure(figsize=(12, 8))  # Larger figure for better readability
sns.lineplot(
    data=filtered_data,
    x="threshold",
    y="SPSS_Size",
    hue="k",
    marker="o",
    linewidth=2,
    alpha=0.8,  # Add transparency to lines
    palette=["#FFC0CB", "#FF69B4",  "#8A2BE2"]  # Use pink, pale pink, and purple for curves
)

# Set y-axis range to expand further beyond the data range
y_min = filtered_data["SPSS_Size"].min() - 150000
y_max = filtered_data["SPSS_Size"].max() + 150000
plt.ylim(y_min, y_max)

# Set x-axis ticks to display 1, 2, 3, 4, 5 explicitly
plt.xticks([2, 3, 4, 5], fontsize=14)

# Format y-axis to avoid scientific notation
plt.ticklabel_format(style='plain', axis='y')

# Add labels and a title
plt.title(dataset_name, fontsize=18)
plt.xlabel("Solidity Threshold (t)", fontsize=16)
plt.ylabel("Total SPSS Size (|SPSS(K)|)", fontsize=16)
plt.legend(title="Value of k", title_fontsize=14, fontsize=12, loc="upper right")
plt.grid(visible=True, linestyle="--", alpha=0.6)
plt.tight_layout()

# Save the plot to a file
plt.savefig(output_file)
print(f"Plot saved as {output_file}")
