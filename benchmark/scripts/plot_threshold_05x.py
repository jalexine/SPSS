import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import warnings

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pandas")
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")

# Get file path and optional output file name from the command-line arguments
if len(sys.argv) < 2:
    print("Usage: python script.py <file_path> [output_file]")
    sys.exit(1)

file_path = sys.argv[1]
output_file = sys.argv[2] if len(sys.argv) > 2 else "comparison_spss_plot.png"

# Function to process data for SPSS Count and SPSS Size
def process_data(file_path):
    data = pd.read_csv(file_path)
    threshold_col = data.columns[2]  # Solidity threshold (3rd column)
    SPSS_Size_col = data.columns[5]  # SPSS Size (6th column)
    SPSS_Count_col = data.columns[6]  # SPSS Count (7th column)
    k_col = data.columns[1]  # k value (2nd column)

    # Select and clean data
    data = data[[threshold_col, SPSS_Size_col, SPSS_Count_col, k_col]].copy()
    data.columns = ["threshold", "SPSS_Size", "SPSS_Count", "k"]
    data = data[data['k'].isin([21, 31, 41])].copy()
    data["threshold"] = pd.to_numeric(data["threshold"], errors='coerce')
    data["SPSS_Size"] = pd.to_numeric(data["SPSS_Size"], errors='coerce')
    data["SPSS_Count"] = pd.to_numeric(data["SPSS_Count"], errors='coerce')
    return data

# Process the input file
data = process_data(file_path)

# Extract dataset name for the title
dataset_name = os.path.basename(file_path).replace(".csv", "").replace("stats_", "").rsplit("_", 1)[0]

# Plot the two graphs side by side
fig, axes = plt.subplots(1, 2, figsize=(18, 8), sharey=False, constrained_layout=True)

# Plot SPSS Count (Left Panel)
sns.lineplot(
    data=data,
    x="threshold",
    y="SPSS_Count",
    hue="k",
    marker="o",
    linewidth=2,
    alpha=0.8,
    palette=["#FFC0CB", "#FF69B4", "#8A2BE2"],
    ax=axes[0]
)
axes[0].set_title(f"{dataset_name}", fontsize=18)
axes[0].set_xlabel("Solidity Threshold (t)", fontsize=16)
axes[0].set_ylabel("SPSS Count", fontsize=16)
axes[0].set_xticks([2, 3, 4, 5])
axes[0].grid(visible=True, linestyle="--", alpha=0.6)
axes[0].legend(title="Value of k", title_fontsize=14, fontsize=12, loc="upper right")

# Plot SPSS Size (Right Panel)
sns.lineplot(
    data=data,
    x="threshold",
    y="SPSS_Size",
    hue="k",
    marker="o",
    linewidth=2,
    alpha=0.8,
    palette=["#FFC0CB", "#FF69B4", "#8A2BE2"],
    ax=axes[1]
)
axes[1].set_title(f"{dataset_name}", fontsize=18)
axes[1].set_xlabel("Solidity Threshold (t)", fontsize=16)
axes[1].set_ylabel("SPSS Size", fontsize=16)
axes[1].set_xticks([2, 3, 4, 5])
axes[1].grid(visible=True, linestyle="--", alpha=0.6)
axes[1].legend(title="Value of k", title_fontsize=14, fontsize=12, loc="upper right")

# Save the figure
plt.savefig(output_file)
print(f"Comparison plot saved as {output_file}")

