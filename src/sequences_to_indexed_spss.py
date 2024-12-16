from collections import defaultdict
from fmi import FmIndex
import sys
from timer import Timer
import csv
import argparse
import random
import pickle
import os


def canonical_kmer(kmer):
    """
    Determines the canonical form of a k-mer by selecting the
    lexicographically smaller between the k-mer and its reverse complement

    Args:
        kmer (str): The k-mer sequence.

    Returns:
        str: The canonical k-mer.
    """
    # Create translation table
    complement_table = str.maketrans("ACGT", "TGCA")
    # Compute reverse complement
    rev_comp = kmer.translate(complement_table)[::-1]
    # Return the lexicographically smaller k-mer
    return kmer if kmer <= rev_comp else rev_comp


def count_kmers(fasta_file, k):
    """
    Counts canonical k-mers directly from a FASTA file,
    processing sequences in chunks.

    Args:
        fasta_file (str): The path to the FASTA file.
        k (int): The size of the k-mers.

    Returns:
        dict: A dictionary with canonical k-mers as keys
        and their counts as values.
    """
    kmer_counts = defaultdict(int)
    # faster to use str.maketrans than creating a fonction
    complement_table = str.maketrans("ACGT", "TGCA")

    with open(fasta_file, "r") as f:
        current_sequence = []

        for line in f:
            if line.startswith(">"):  # Header, process accumulated sequence
                if current_sequence:
                    sequence = "".join(current_sequence)
                    for i in range(len(sequence) - k + 1):
                        kmer = sequence[i:i + k]
                        rev_kmer = kmer.translate(complement_table)[::-1]
                        canonical = kmer if kmer <= rev_kmer else rev_kmer
                        kmer_counts[canonical] += 1
                    current_sequence = []
            else:
                current_sequence.append(line.strip())

        # Process the final sequence
        if current_sequence:
            sequence = "".join(current_sequence)
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i + k]
                rev_kmer = kmer.translate(complement_table)[::-1]
                canonical = kmer if kmer <= rev_kmer else rev_kmer
                kmer_counts[canonical] += 1

    return kmer_counts


def filter_kmers(kmer_counts, threshold):
    """
    Filters out k-mers that do not meet the specified solidity threshold.

    Args:
        kmer_counts (dict): A dictionary of k-mer counts.
        threshold (int): The minimum count a k-mer must have to be retained.

    Returns:
        set: A set of k-mers that meet or exceed the solidity threshold.
    """
    total_kmers = sum(kmer_counts.values())
    if threshold > total_kmers:
        print(
            f"error : solidity threshold ({threshold}) exceeds the total "
            f"number of k-mers ({total_kmers})."
        )
        sys.exit(1)
    return {kmer for kmer, count in kmer_counts.items() if count >= threshold}


def generate_simplitigs(kmers, k):
    """
    Generates maximal simplitigs from a set of solid k-mers.

    Args:
        kmers (set): A set of solid k-mers.
        k (int): The size of the k-mers.

    Returns:
        list: A list of maximal simplitigs.
    """
    maximal_simplitigs = []

    def extend_simplitig(K, simplitig, direction):
        """
        Extends a simplitig in a given direction (forward or backward).

        Args:
            K (set): The set of remaining k-mers.
            simplitig (str): The current simplitig.
            direction (str): Direction of extension ('forwards' or 'backwards')

        Returns:
            tuple: Updated set of k-mers and the extended simplitig.
        """
        extending = True
        while extending:
            extending = False
            if direction == 'forwards':
                q = simplitig[-(k - 1):]  # Get the last k-1 bases
                extend_kmer_func = lambda x: q + x
            elif direction == 'backwards':
                q = simplitig[:k - 1]  # Get the first k-1 bases
                extend_kmer_func = lambda x: x + q

            # Try extending with each possible base ('A', 'C', 'G', 'T')
            for x in ['A', 'C', 'G', 'T']:
                kmer = extend_kmer_func(x)
                # Convert to canonical form
                canonical_kmer_form = canonical_kmer(kmer)
                if canonical_kmer_form in K:
                    extending = True
                    simplitig = simplitig + x if direction == 'forwards' else x + simplitig
                    K.remove(canonical_kmer_form)  # Remove canonical k-mer
                    break
        return K, simplitig

    # Iterate through the k-mers and generate simplitigs
    while len(kmers) > 0:
        # Pick a k-mer to seed the simplitig
        seed_kmer = random.choice(list(kmers))
        # Remove the canonical form of the seed k-mer
        kmers.remove(canonical_kmer(seed_kmer))

        # Initialize the simplitig with the seed k-mer
        simplitig = canonical_kmer(seed_kmer)
        kmers, simplitig = extend_simplitig(kmers, simplitig, 'backwards')
        kmers, simplitig = extend_simplitig(kmers, simplitig, 'forwards')

        maximal_simplitigs.append(simplitig)

    return maximal_simplitigs


def generate_unitigs(kmers, k):
    """
    Generates unitigs dynamically by traversing the graph.

    Args:
        kmers (set): A set of solid k-mers.
        k (int): The size of the k-mers.

    Returns:
        list: A list of unitigs.
    """
    visited_kmers = set()
    unitigs = []

    def extend_unitig(current_kmer, direction):
        """
        Extends a unitig in a given direction.

        Args:
            current_kmer (str): The current k-mer.
            direction (str): direction of extension ('forward' or 'backward')

        Returns:
            str: The extended unitig.
        """
        unitig = current_kmer
        while True:
            if direction == 'forward':
                suffix = unitig[-(k - 1):]
                next_kmer = None
                for base in 'ACGT':
                    candidate_kmer = suffix + base
                    if candidate_kmer in kmers and candidate_kmer not in visited_kmers:
                        next_kmer = candidate_kmer
                        break
                if next_kmer:
                    visited_kmers.add(next_kmer)
                    unitig += next_kmer[-1]
                else:
                    break
            elif direction == 'backward':
                prefix = unitig[:k - 1]
                next_kmer = None
                for base in 'ACGT':
                    candidate_kmer = base + prefix
                    if candidate_kmer in kmers and candidate_kmer not in visited_kmers:
                        next_kmer = candidate_kmer
                        break
                if next_kmer:
                    visited_kmers.add(next_kmer)
                    unitig = next_kmer[0] + unitig
                else:
                    break
        return unitig

    # Iterate through the k-mers and generate unitigs
    for kmer in kmers:
        if kmer not in visited_kmers:
            visited_kmers.add(kmer)
            # Extend in both directions
            unitig = extend_unitig(kmer, 'backward')
            unitig = extend_unitig(unitig, 'forward')
            unitigs.append(unitig)

    return unitigs


def generate_spss(kmers, k, mode='simplitig'):
    """
    Creates the Solid Prefix Sequence Set (SPSS) by concatenating
    simplitigs or unitigs based on the selected mode.

    Args:
        kmers (set): The set of solid k-mers.
        k (int): The size of the k-mers.
        mode (str, optional): Mode of SPSS construction
        ('simplitig' or 'unitig'). Defaults to 'simplitig'
    Returns:
        str: The concatenated SPSS string.
    """
    # Check the mode to decide which method to use for SPSS construction
    if mode == 'simplitig':
        simplitigs = generate_simplitigs(kmers, k)
    elif mode == 'unitig':
        simplitigs = generate_unitigs(kmers, k)
    # '#' separates individual sequences, and '$' marks the end
    # for compatibility with BWT and indexing.
    return "#".join(simplitigs) + "$"


def prepare_output_dir(base_name, mode, k, t, extension="dump", parent_dir="benchmark"):
    """
    Prepares the output directory for saving the FM-index file.

    Args:
        base_name (str): The base name derived from the FASTA file.
        mode (str): The mode used for SPSS construction
            ('simplitig' or 'unitig').
        k (int): The size of the k-mers.
        t (int): The solidity threshold for k-mers.
        extension (str, optional): The file extension.
            Defaults to "dump".
        parent_dir (str, optional): parent directory where outputs are stored
            Defaults to "benchmark".

    Returns:
        str: The full path to the output file.
    """
    # Create the output directory if it doesn't exist
    output_dir = os.path.abspath(os.path.join(parent_dir, base_name))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Format the file name with relevant parameters
    file_name = f"fmi_{base_name}_{mode}_k{k}_t{t}.{extension}"
    # Return the full path to the file
    return os.path.join(output_dir, file_name)


def save_fm_index(fm_index, fasta_file, mode, k, t, output_file=None):
    """
    Serializes and saves the FM-index to a file.
    Automatically creates folders if needed.

    Args:
        fm_index (FmIndex): The FM-index to save.
        fasta_file (str): The path to the FASTA file.
        mode (str): mode for SPSS construction ('simplitig' or 'unitig').
        k (int): The size of the k-mers.
        t (int): The solidity threshold for k-mers.
        output_file (str, optional): The path to save the FM-index.
            If not provided or is a directory, a default file name is generated

    Returns:
        str: The path where the FM-index was saved.
    """
    # Extract the base name from the FASTA file
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]

    # Generate default output file name if none or it's a directory
    if not output_file or os.path.isdir(output_file):
        output_file = prepare_output_dir(base_name, mode, k, t, extension="dump")

    # Ensure the parent directory exists
    parent_dir = os.path.dirname(output_file)
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

    # Save the FM-index to the specified file
    try:
        with open(output_file, 'wb') as f:
            pickle.dump(fm_index, f)
        print(f"FM-index saved to {output_file}")
    except Exception as e:
        print(f"Error saving FM-index: {e}")

    return output_file


def save_benchmark_results(fasta_file, mode, stats, k, threshold):
    """
    Saves benchmark results to a CSV file in the stats directory
    and prints them in their original format.

    Args:
        fasta_file (str): The input FASTA file.
        mode (str): mode used for SPSS construction ('simplitig' or 'unitig')
        stats (dict): A dictionary containing benchmarking statistics.
        k (int): The size of the k-mers.
        threshold (int): The solidity threshold for k-mers.

    Returns:
        None
    """
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    # Unique directory for stats
    parent_dir = os.path.abspath("benchmark/stats")

    if not os.path.exists(parent_dir):
        # Create the "stats" folder if it doesn't exist
        os.makedirs(parent_dir)

    # CSV file name based only on the mode and dataset
    csv_file_name = os.path.join(parent_dir, f"stats_{base_name}_{mode}.csv")

    # Prepare the processed stat
    processed_stats = {
        'Dataset': base_name,
        'k': k,
        't': threshold,
        'TIME_SELECTING_KMERS (sec)': float(stats['TIME_SELECTING_KMERS'].split()[0]),
        'TIME_SPSS_CONSTRUCTION (sec)': float(stats['TIME_SPSS_CONSTRUCTION'].split()[0]),
        'SPSS(K)': int(stats['SPSS(K)'].split()[0]),
        'SPSS(K) count': int(stats['SPSS(K) count'].split()[0]),
        'TIME_BUILD_FMI (sec)': float(stats['TIME_BUILD_FMI'].split()[0]),
    }

    fieldnames = list(processed_stats.keys())

    # Write to the CSV file
    with open(csv_file_name, mode='a', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        # If the file is empty, write the header
        if csv_file.tell() == 0:
            writer.writeheader()
        writer.writerow(processed_stats)

    print(f"Benchmark results saved to {csv_file_name}")


def test_fm_index(fm_index, spss, filtered_kmers, sample_size=100):
    """
    Tests the FM-index by comparing results from its `contains` method
    with Python's `in` operator on the SPSS string.

    Args:
        fm_index (FmIndex): The FM-index to validate.
        spss (str): The SPSS string for brute-force validation.
        filtered_kmers (set): The set of solid k-mers to test.
        sample_size (int, optional): Number of k-mers to sample for testing.
            Defaults to 100.

    Returns:
        bool: True if all tested k-mers match between FM-index and `in`, False otherwise.
    """
    # Randomly sample up to `sample_size` k-mers from the filtered k-mers
    sampled_kmers = random.sample(list(filtered_kmers), min(sample_size, len(filtered_kmers)))

    # Validate each k-mer using FM-index and Python's `in` operator
    for kmer in sampled_kmers:
        fm_contains = fm_index.contains(kmer)  # Query the FM-index
        # Use Python's `in` operator on SPSS
        python_contains = (kmer in spss)
        # If there's a mismatch, report it and fail the test
        if fm_contains != python_contains:
            print(f"Mismatch detected for k-mer '{kmer}':")
            print(f"  FM-index result: {fm_contains}")
            print(f"  Python `in` result: {python_contains}")
            return False

    # Return True if all sampled k-mers match
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate and serialize the FM-index of SPSS.",
        usage="python src/sequences_to_indexed_spss.py -i <input.fasta> -k <kmer_size> -t <threshold> [-q <output_file>] [-m {simplitig,unitig}]"
    )
    parser.add_argument("-i", required=True, help="Input FASTA file containing genomic sequences.")
    parser.add_argument("-k", type=int, required=True, help="Size of k-mers.")
    parser.add_argument("-t", type=int, required=True, help="Solidity threshold for k-mers.")
    parser.add_argument("-o", required=False, help="Output file for serialized FM-index")
    parser.add_argument("-stats", required=False, help="Output file for benchmark stats (optional).")
    parser.add_argument("-m", choices=['simplitig', 'unitig'], default='simplitig', help="Mode of SPSS construction ('simplitig' default or 'unitig').")

    # Capture missing arguments error
    try:
        args = parser.parse_args()
    except SystemExit:
        print("\n\n \033[95m♡ pls use the -h or --help option for usage details \033[0m")
        sys.exit(1)

    fasta_file = args.i
    k = args.k
    threshold = args.t
    mode = args.m
    output_file = args.o  # Optional - can be passed or left empty
    stats_file = args.stats  # csv optionnal

    stats = {}

    # Step 1: Count k-mers
    with Timer() as total_time:
        kmer_counts = count_kmers(fasta_file, k)
    stats['TIME_SELECTING_KMERS'] = f"{round(total_time.t, 2)} seconds"

    # Step 2: Filter solid k-mers
    filtered_kmers = filter_kmers(kmer_counts, threshold)

    # Step 3: Generate SPSS
    with Timer() as total_time:
        spss = generate_spss(filtered_kmers, k, mode=mode)
    spss_size = len(spss)
    num_spss = spss.count("#") + 1
    stats['TIME_SPSS_CONSTRUCTION'] = f"{round(total_time.t, 2)} seconds"

    stats['SPSS(K)'] = f"{spss_size} characters"
    stats['SPSS(K) count'] = f"{num_spss} sequences"

    # Step 4: Build and validate the FM-index
    with Timer() as total_time:
        fm_index = FmIndex(spss)
    stats['TIME_BUILD_FMI'] = f"{round(total_time.t, 2)} seconds"

    # Test FM-index
    if test_fm_index(fm_index, spss, filtered_kmers):
        print("FM-index validation passed!")
    else:
        print("Oh noooo FM-index validation failed!")
        sys.exit(1)

    # Step 5: Save the FM-index if specified
    if output_file:
        save_fm_index(fm_index, fasta_file, mode, k, threshold, output_file=output_file)
    else:
        print("FM-index not saved because -q was not specified.")

    # Step 6: Save benchmark stats if specified
    if stats_file:
        save_benchmark_results(fasta_file, mode, stats, k, threshold)
    else:
        print("Benchmark stats not saved because -stats was not specified.")

    # Print stats summary to console
    print("\nBenchmark Summary:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
    print("-" * 50)