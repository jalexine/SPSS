import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import warnings

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pandas")
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")

# Get file paths and optional output file name from the command-line arguments
if len(sys.argv) < 3:
    print("Usage: python script.py <file_path1> <file_path2> [output_file]")
    sys.exit(1)

file_path1 = sys.argv[1]
file_path2 = sys.argv[2]
output_file = sys.argv[3] if len(sys.argv) > 3 else "comparison_plot.png"

# Function to process data
def process_data(file_path):
    data = pd.read_csv(file_path)
    threshold_col = data.columns[2]
    SPSS_Count_col = data.columns[6]
    k_col = data.columns[1]

    # Select and clean data
    data = data[[threshold_col, SPSS_Count_col, k_col]].copy()
    data.columns = ["threshold", "SPSS_Count", "k"]
    data = data[data['k'].isin([21, 31, 41])].copy()
    data["threshold"] = pd.to_numeric(data["threshold"], errors='coerce')
    data["SPSS_Count"] = pd.to_numeric(data["SPSS_Count"], errors='coerce')
    return data

# Process both CSV files
data1 = process_data(file_path1)
data2 = process_data(file_path2)

# Extract dataset names for titles
dataset_name1 = os.path.basename(file_path1).replace(".csv", "").replace("stats_", "").rsplit("_", 1)[0]
dataset_name2 = os.path.basename(file_path2).replace(".csv", "").replace("stats_", "").rsplit("_", 1)[0]

# Plot the two graphs side by side
fig, axes = plt.subplots(1, 2, figsize=(18, 8), sharey=True)
# Plot for the first dataset
sns.lineplot(
    data=data1,
    x="threshold",
    y="SPSS_Count",
    hue="k",
    marker="o",
    linewidth=2,
    alpha=0.8,
    palette=["#FFC0CB", "#FF69B4", "#8A2BE2"],
    ax=axes[0]
)
axes[0].set_title(dataset_name1, fontsize=18)
axes[0].set_xlabel("Solidity Threshold (t)", fontsize=16)
axes[0].set_ylabel("SPSS(K) Character Count", fontsize=16)
axes[0].set_xticks([2, 3, 4, 5])
axes[0].grid(visible=True, linestyle="--", alpha=0.6)

# Move legend inside the first plot
axes[0].legend(title="Value of k", title_fontsize=14, fontsize=12, loc="upper right")

# Plot for the second dataset
sns.lineplot(
    data=data2,
    x="threshold",
    y="SPSS_Count",
    hue="k",
    marker="o",
    linewidth=2,
    alpha=0.8,
    palette=["#FFC0CB", "#FF69B4", "#8A2BE2"],
    ax=axes[1]
)
axes[1].set_title(dataset_name2, fontsize=18)
axes[1].set_xlabel("Solidity Threshold (t)", fontsize=16)
axes[1].set_ylabel("")  # Avoid redundant y-axis label
axes[1].set_xticks([2, 3, 4, 5])
axes[1].grid(visible=True, linestyle="--", alpha=0.6)

# Move legend inside the second plot
axes[1].legend(title="Value of k", title_fontsize=14, fontsize=12, loc="upper right")

# Adjust layout to accommodate the global legend
plt.tight_layout(rect=[0, 0, 1, 0.9])  # Leave space at the top for the legend

plt.savefig(output_file)
print(f"Comparison plot saved as {output_file}")
