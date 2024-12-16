import os
import re
import pandas as pd
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Parse log files and extract relevant data.")
parser.add_argument("input_dir", type=str, help="Path to the directory containing log files.")
parser.add_argument("output_file", type=str, help="Path to save the filtered CSV output.")
args = parser.parse_args()

# Path to the directory containing logs
log_dir = args.input_dir

# Initialize data storage
data = []

# Parse each log file in the input directory
for log_file in os.listdir(log_dir):
    if log_file.endswith(".log"):
        file_path = os.path.join(log_dir, log_file)
        with open(file_path, "r") as f:
            content = f.read()

            # Extract main information from the log file content
            dataset_match = re.search(r"Running query: dataset=(.+?), mode=(.+?), k=(\d+), t=(\d+)", content)
            if dataset_match:
                dataset = dataset_match.group(1)
                mode = dataset_match.group(2)
                k = int(dataset_match.group(3))
                t = int(dataset_match.group(4))

                # Extract execution time (Elapsed time, m:ss)
                time_match = re.search(r"Elapsed \(wall clock\) time.*?:\s+(\d+):(\d+\.\d+)", content)
                if time_match:
                    minutes = int(time_match.group(1))
                    seconds = float(time_match.group(2))
                    elapsed_time_formatted = f"{minutes:02}:{seconds:.2f}"
                else:
                    elapsed_time_formatted = "00:00.00"

                # Extract peak memory usage (Maximum resident set size in kbytes)
                mem_match = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", content)
                if mem_match:
                    memory = int(mem_match.group(1)) / 1024  # Convert to MB
                else:
                    memory = None

                # Append collected data to the list
                data.append({
                    "Dataset": dataset,
                    "k": k,
                    "t": t,
                    "Mode": mode,
                    "Elapsed Time (mm:ss)": elapsed_time_formatted,
                    "Peak Mem. (MB)": round(memory, 2) if memory else None
                })

# Convert the collected data into a DataFrame
df = pd.DataFrame(data)

# Apply filters:
# - Include k=21 or k=31 when t=2
# - Include k=31 or k=41 when t=3
df_filtered = df[
    ((df["k"].isin([21, 31])) & (df["t"] == 2)) |
    ((df["k"].isin([31, 41])) & (df["t"] == 3))
]

# Sort the data by Dataset, k, t, and Mode
df_filtered = df_filtered.sort_values(by=["Dataset", "k", "t", "Mode"])

# Save the filtered data into a clean CSV file
df_filtered.to_csv(args.output_file, index=False)

# Display a preview of the filtered data
print(df_filtered)
