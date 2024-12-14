from collections import defaultdict
from fmi import FmIndex
import time
import sys
from timer import Timer
import csv
import argparse
import random
import pickle
import os


def reverse_complement(kmer):
    """
    Compute the reverse complement of a k-mer using translation tables.

    Args:
        kmer (str): The k-mer sequence.

    Returns:
        str: The reverse complement of the input k-mer.
    """
    return kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical_kmer(kmer):
    """
    Determines the canonical form of a k-mer by selecting the lexicographically smaller between the k-mer and its reverse complement.

    Args:
        kmer (str): The k-mer sequence

    Returns:
        str: The canonical k-mer.
    """
    rev_comp = reverse_complement(kmer)
    return kmer if kmer <= rev_comp else rev_comp


def count_kmers(fasta_file, k):
    """
    Counts canonical k-mers directly from a FASTA file, processing sequences in chunks.

    Args:
        fasta_file (str): The path to the FASTA file.
        k (int): The size of the k-mers.

    Returns:
        dict: A dictionary with canonical k-mers as keys and their counts as values.
    """
    kmer_counts = defaultdict(int)
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
        print(f"error : solidity threshold ({threshold}) exceeds the total number of k-mers ({total_kmers}).")
        sys.exit(1)
    return {kmer for kmer, count in kmer_counts.items() if count >= threshold}

def remove_reverse_complements(kmers):
    """
    Removes redundant reverse complements by keeping only the canonical k-mers.

    Args:
        kmers (set): A set of k-mers (strings).

    Returns:
        set: A cleaned set containing only canonical k-mers.
    """
    cleaned_kmers = set() 
    for kmer in kmers:
        canonical = canonical_kmer(kmer)  # Canonical k-mer
        cleaned_kmers.add(canonical)  # add only canonical k-mer
    return cleaned_kmers

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
            direction (str): Direction of extension ('forwards' or 'backwards').

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
                canonical_kmer_form = canonical_kmer(kmer)
                if canonical_kmer_form in K:
                    extending = True
                    simplitig = simplitig + x if direction == 'forwards' else x + simplitig
                    K.remove(canonical_kmer_form)  # remove canonical kmer
                    break
        return K, simplitig

    kmers = remove_reverse_complements(kmers)

    # Iterate through the k-mers and generate simplitigs
    while len(kmers) > 0:
        seed_kmer = random.choice(list(kmers)) # Pick a k-mer to seed the simplitig
        kmers.remove(seed_kmer)  # Remove the seed k-mer from the set

        # Initialize the simplitig with the seed k-mer
        simplitig = seed_kmer
        kmers, simplitig = extend_simplitig(kmers, simplitig, 'backwards')  # Extend backwards

        kmers, simplitig = extend_simplitig(kmers, simplitig, 'forwards')  # Extend forwards

        maximal_simplitigs.append(simplitig)  

    return maximal_simplitigs


def generate_unitigs(kmers):
    """
    Generates unitigs from a set of k-mers.

    A unitig is a maximal sequence formed by traversing the k-mer graph, 
    where each prefix leads to exactly one suffix and the path is extended 
    while there is only one valid neighbor.

    Args:
        kmers (set): The set of k-mers.

    Returns:
        list: A list of unitigs.
    """
    graph = defaultdict(list)
    for kmer in kmers:
        prefix = kmer[:-1]  # Prefix of the k-mer (first k-1 bases)
        suffix = kmer[1:]   # Suffix of the k-mer (last k-1 bases)
        graph[prefix].append(suffix)

    unitigs = []
    visited_edges = set()

    # Traverse the graph to find unitigs
    for node in graph:
        for neighbor in graph[node]:
            edge = (node, neighbor)
            if edge in visited_edges:
                continue

            visited_edges.add(edge)
            unitig = [node, neighbor]

            # Extend the unitig while there is exactly one valid neighbor
            current = neighbor
            while current in graph and len(graph[current]) == 1:
                next_node = graph[current][0]
                next_edge = (current, next_node)
                if next_edge in visited_edges:
                    break
                visited_edges.add(next_edge)
                unitig.append(next_node)
                current = next_node

            # Join unitig and add to list
            unitigs.append(''.join([unitig[0]] + [u[-1] for u in unitig[1:]]))

    return unitigs


def generate_spss(kmers, k, mode='simplitig'):
    """
    Creates the Solid Prefix Sequence Set (SPSS) by concatenating simplitigs or unitigs based on the selected mode.

    Args:
        kmers (set): The set of solid k-mers.
        k (int): The size of the k-mers.
        mode (str, optional): Mode of SPSS construction ('simplitig' or 'unitig'). Defaults to 'simplitig'.

    Returns:
        str: The concatenated SPSS string.
    """
    if mode == 'simplitig':
        simplitigs = generate_simplitigs(kmers, k)
    elif mode == 'unitig':
        simplitigs = generate_unitigs(kmers)
    return "#".join(simplitigs) + "$"


def prepare_output_dir(base_name, mode, k, t, extension="dump", parent_dir="benchmark"):
    """
    Prepares the output directory for saving the FM-index file.

    Args:
        base_name (str): The base name derived from the FASTA file.
        mode (str): The mode used for SPSS construction ('simplitig' or 'unitig').
        k (int): The size of the k-mers.
        t (int): The solidity threshold for k-mers.
        extension (str, optional): The file extension. Defaults to "dump".
        parent_dir (str, optional): The parent directory where outputs are stored. Defaults to "benchmark".

    Returns:
        str: The full path to the output file.
    """
    output_dir = os.path.abspath(os.path.join(parent_dir, base_name))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_name = f"fmi_{base_name}_{mode}_k{k}_t{t}.{extension}"
    return os.path.join(output_dir, file_name)


def save_fm_index(fm_index, fasta_file, mode, k, t, output_file=None):
    """
    Serializes and saves the FM-index to a file. Automatically creates folders if needed.

    Args:
        fm_index (FmIndex): The FM-index to save.
        fasta_file (str): The path to the FASTA file.
        mode (str): The mode used for SPSS construction ('simplitig' or 'unitig').
        k (int): The size of the k-mers.
        t (int): The solidity threshold for k-mers.
        output_file (str, optional): The path to save the FM-index. If not provided or is a directory, a default file name is generated.

    Returns:
        str: The path where the FM-index was saved.
    """
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]

    if not output_file or os.path.isdir(output_file):
        output_file = prepare_output_dir(base_name, mode, k, t, extension="dump")

    try:
        with open(output_file, 'wb') as f:
            pickle.dump(fm_index, f)
        print(f"FM-index saved to {output_file}")
    except Exception as e:
        print(f"Error saving FM-index: {e}")

    return output_file


def save_benchmark_results(fasta_file, mode, stats):
    """
    Saves benchmark results to a CSV file in the stats directory and prints them in their original format.

    Args:
        fasta_file (str): The input FASTA file.
        mode (str): The mode used for SPSS construction ('simplitig' or 'unitig').
        stats (dict): A dictionary containing benchmarking statistics.

    Returns:
        None
    """
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    parent_dir = os.path.abspath("benchmark/stats")  # Unique directory for stats

    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)  # Create the "stats" folder if it doesn't exist

    # CSV file name based only on the mode and dataset
    csv_file_name = os.path.join(parent_dir, f"stats_{base_name}_{mode}.csv")

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

        if csv_file.tell() == 0:  # If the file is empty, write the header
            writer.writeheader()
        writer.writerow(processed_stats)

    print(f"Benchmark results saved to {csv_file_name}")


def test_fm_index(fm_index, spss, filtered_kmers, sample_size=100):
    """
    Tests the correctness of the FM-index by sampling k-mers and comparing the results with a brute-force search.

    Args:
        fm_index (FmIndex): The FM-index to test.
        spss (str): The SPSS to test against.
        filtered_kmers (set): The set of filtered k-mers.
        sample_size (int, optional): The number of k-mers to sample for testing. Defaults to 100.

    Returns:
        bool: True if the FM-index passes the test, False otherwise.
    """
    sampled_kmers = random.sample(list(filtered_kmers), min(sample_size, len(filtered_kmers)))

    for kmer in sampled_kmers:
        fm_contains = fm_index.contains(kmer)
        brute_contains = (kmer in spss)

        if fm_contains != brute_contains:
            print(f"test fm-index : its so oveeeeer '{kmer}'")
            print(f"FM-index: {fm_contains}, Brute-force search: {brute_contains}")
            return False

    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate and serialize the FM-index of SPSS.",
        usage="python src/sequences_to_indexed_spss.py -i <input.fasta> -k <kmer_size> -t <threshold> [-o <output_file>] [-m {simplitig,unitig}]"
    )    
    parser.add_argument("-i", required=True, help="Input FASTA file containing genomic sequences.")
    parser.add_argument("-k", type=int, required=True, help="Size of k-mers.")
    parser.add_argument("-t", type=int, required=True, help="Solidity threshold for k-mers.")
    parser.add_argument("-o", required=False, help="Output file for serialized FM-index.")  # Make it optional
    parser.add_argument("-m", choices=['simplitig', 'unitig'], default='simplitig', help="Mode of SPSS construction ('simplitig' default or 'unitig').")

    # Capture missing arguments error
    try:
        args = parser.parse_args()
    except SystemExit:
        print("\n\n \033[95m♡ Please use the -h or --help option for usage details \033[0m")
        sys.exit(1)

    fasta_file = args.i
    k = args.k
    threshold = args.t
    mode = args.m
    output_file = args.o  # Optional - can be passed or left empty

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
    # You can use the test if you need it.. 
    # if not test_fm_index(fm_index, spss, filtered_kmers):
        # print("Oh no! FM-index validation failed.")
        # exit(1)

    # Step 5: Serialize the FM-index
    if output_file and output_file.endswith(".dump"):
        dump_file_name = save_fm_index(fm_index, fasta_file, mode, k, threshold, output_file=output_file)
    else:
        dump_file_name = save_fm_index(fm_index, fasta_file, mode, k, threshold, output_file=None)
        save_benchmark_results(fasta_file, mode, stats)

    for key, value in stats.items():
        print(f"  {key}: {value}")
    print("-" * 50)  # Prints a line with 50 dashes
