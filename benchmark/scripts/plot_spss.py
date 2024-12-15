import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

# Function to generate the plot
def plot_spss_two(input_csv1, input_csv2, output_png):
    plt.figure(figsize=(10, 6))

    # Colors for the two datasets
    colors = ['#FF69B4','#800080']  # Violet, Pink
    input_csvs = [input_csv1, input_csv2]

    y_min = 0
    y_max_values = []

    for idx, input_csv in enumerate(input_csvs):
        # Read the CSV file
        data = pd.read_csv(input_csv)

        # Filter for t = 3 (adjust the filter as necessary)
        data_t2 = data[data['t'] == 2]

        # Select the columns (k = column 2, SPSS size = column 6)
        data_t2 = data_t2[['k', 'SPSS(K)']].copy()

        # Rename columns for better readability
        data_t2.columns = ['k', 'SPSS_Size']

        # Plot the curve
        label = os.path.basename(input_csv)
        label = label.replace("_unitig.csv", "").replace("_simplitig.csv", "").replace("stats_", "")

        plt.plot(
            data_t2['k'],
            data_t2['SPSS_Size'],
            label=label,
            marker='o',
            linewidth=2,
            color=colors[idx]
        )

        # Collect y_max for dynamic scaling
        y_max_values.append(data_t2['SPSS_Size'].max())

    # Dynamically calculate y-axis limits and ticks for small or large values
    y_max = max(y_max_values)
    if y_max < 155000:  # Automatic zoom for values < 100k
        y_step = 15000  # Smaller step size for better visibility
        y_min = int(np.floor(min(y_max_values) / y_step)) * y_step
        y_max = int(np.ceil(y_max / y_step)) * y_step  # Round up to nearest step
        y_ticks = np.arange(y_min, y_max + y_step, y_step)
        if len(y_ticks) < 5:  # S'assure qu'il y a au moins 5 ticks
            y_ticks = np.linspace(y_min, y_max, 5)
    else:
        y_step = max(10000, (y_max - y_min) // 5)
        y_max = int(np.ceil(y_max / y_step)) * y_step
        y_ticks = np.arange(y_min, y_max + y_step, y_step)

    # Set x-axis ticks and limits
    plt.xticks([21, 31, 41, 51, 61, 71], fontsize=12)
    plt.xlim(19, 73)  # Set bounds for the horizontal axis

    # Set y-axis ticks and limits
    plt.yticks(y_ticks, fontsize=12)
    plt.ylim(y_min, y_max)  # Automatically zoom for small values

    # Display the true values on the y-axis without scientific notation
    plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))  # Format with commas

    # Add title, labels, and legend
    plt.xlabel("k", fontsize=14)
    plt.ylabel("Total SPSS Size", fontsize=14)
    plt.legend(title="Datasets", loc='upper left', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # Save the plot to the specified file
    plt.savefig(output_png)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate SPSS plots from two CSV files.")
    parser.add_argument("input_csv1", help="Path to the first input CSV file.")
    parser.add_argument("input_csv2", help="Path to the second input CSV file.")
    parser.add_argument("output_png", help="Path to the output PNG file.")

    args = parser.parse_args()
    plot_spss_two(args.input_csv1, args.input_csv2, args.output_png)
