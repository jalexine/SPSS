from collections import defaultdict
from fmi import FmIndex
import time
from timer import Timer
import csv
import argparse
import random
import pickle
import os


def reverse_complement(kmer):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[base] for base in reversed(kmer))

def canonical_kmer(kmer):
    rev_comp = reverse_complement(kmer)
    return min(kmer, rev_comp)

def count_kmers(fasta_file, k):
    kmer_counts = defaultdict(int)
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'): 
                for i in range(len(line) - k + 1):
                    kmer = line[i:i+k]
                    kmer = canonical_kmer(kmer)  
                    kmer_counts[kmer] += 1
    return kmer_counts


def filter_kmers(kmer_counts, threshold):
    return {kmer for kmer, count in kmer_counts.items() if count >= threshold}
from collections import defaultdict

def generate_simplitigs(kmers, k):
    """
    Generates maximal simplitigs from a set of solid k-mers.

    Args:
        kmers (set): The set of solid k-mers.
        k (int): The size of the k-mers.

    Returns:
        list: A list of maximal simplitigs.
    """
    maximal_simplitigs = []
    visited_kmers = set()

    def extend_simplitig(K, simplitig, direction):
        """
        Extends a simplitig in the specified direction (forwards or backwards).

        Args:
            K (set): The set of remaining k-mers.
            simplitig (str): The current simplitig.
            direction (str): The direction to extend ('forwards' or 'backwards').

        Returns:
            tuple: Updated (K, simplitig).
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
                if kmer in K:
                    extending = True
                    simplitig = simplitig + x if direction == 'forwards' else x + simplitig
                    K.remove(kmer)
                    break
        return K, simplitig

    # Iterate through the k-mers and generate simplitigs
    while len(kmers) > 0:
        seed_kmer = next(iter(kmers))  # Pick a k-mer to seed the simplitig
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
    if mode == 'simplitig':
        simplitigs = generate_simplitigs(kmers, k)
    elif mode == 'unitig':
        simplitigs = generate_unitigs(kmers)
    return "#".join(simplitigs) + "$"

def save_fm_index(fm_index, fasta_file, mode, output_file=None):
    """
    Serializes and saves the FM-index to a file. Automatically creates folders if needed.

    Args:
        fm_index (FmIndex): The FM-index to save.
        fasta_file (str): The path to the FASTA file.
        mode (str): The mode used for SPSS construction ('simplitig' or 'unitig').
        output_file (str, optional): The path to save the FM-index. If not provided, a default file name is generated.

    Returns:
        str: The path where the FM-index was saved.
    """
    parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

    # Dynamically generate the subfolder name based on the FASTA file name
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    fasta_folder = f"{base_name}_reads"  # Customize this as needed

    # Create the folder "benchmark/<dynamic_folder_name>" if it doesn't exist
    reads_folder = os.path.join(parent_dir, 'benchmark', fasta_folder)
    if not os.path.exists(reads_folder):
        os.makedirs(reads_folder)

    # If no output file is provided, generate the default name
    if not output_file:
        output_file = os.path.join(reads_folder, f"{base_name}_{mode}_{k}_{threshold}.dump")
    
    # If output_file is a relative path, prepend it with the benchmark folder
    elif os.path.isabs(output_file): 
        pass
    else:
        output_file = os.path.join(reads_folder, output_file)

    # Save the FM-index to the generated file path
    try:
        with open(output_file, 'wb') as f:
            pickle.dump(fm_index, f)
        print(f"FM-index saved to {output_file}")
    except Exception as e:
        print(f"Error saving FM-index: {e}") 
    return output_file




def save_benchmark_results(fasta_file, mode, stats):
    """
    Saves benchmark results to a CSV file and prints them to the console.

    Args:
        fasta_file (str): The input FASTA file.
        mode (str): The mode used for SPSS construction ('simplitig' or 'unitig').
        stats (dict): A dictionary containing benchmarking statistics.

    Returns:
        None
    """
    parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]

    csv_file_name = os.path.join(parent_dir, 'benchmark', f"stats_{base_name}.{mode}.csv")

    with open(csv_file_name, mode='a', newline='') as csv_file:
        fieldnames = ['Dataset', 'TIME_SELECTING_KMERS', 'TIME_SPSS_CONSTRUCTION', 'SPSS(K)', 'SPSS(K) count', 'TIME_BUILD_FMI']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        if csv_file.tell() == 0:
            writer.writeheader()
        writer.writerow(stats)
    for key, value in stats.items():
        print(f"{key}: {value}")


def test_fm_index(fm_index, spss, filtered_kmers, sample_size=100):
    """
    Tests the correctness of the FM-index by sampling k-mers and comparing the results with a brute-force search.

    Args:
        fm_index (FmIndex): The FM-index to test.
        spss (str): The SPSS to test against.
        filtered_kmers (set): The set of filtered k-mers.
        sample_size (int): The number of k-mers to sample for testing.

    Returns:
        bool: True if the FM-index passes the test, False otherwise.
    """
    sampled_kmers = random.sample(list(filtered_kmers), min(sample_size, len(filtered_kmers)))

    for kmer in sampled_kmers:
        fm_contains = fm_index.contains(kmer)

        brute_contains = (kmer in spss)

        if fm_contains != brute_contains:
            print(f"test fm-index : its so oveeeeer '{kmer}'")
            print(f"FM-index: {fm_contains}, Recherche brute: {brute_contains}")
            return False

    return True



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and serialize the FM-index of SPSS.")
    parser.add_argument("-i", required=True, help="Input FASTA file containing genomic sequences.")
    parser.add_argument("-k", type=int, required=True, help="Size of k-mers.")
    parser.add_argument("-t", type=int, required=True, help="Solidity threshold for k-mers.")
    parser.add_argument("-o", required=False, help="Output file for serialized FM-index.")  # Make it optional
    parser.add_argument("-m", choices=['simplitig', 'unitig'], default='simplitig', help="Mode of spss construction ('simplitig' or 'unitig').")

    args = parser.parse_args()

    fasta_file = args.i
    k = args.k
    threshold = args.t
    mode = args.m
    output_file = args.o  # Optional - can be passed or left empty

    stats = {}
    print(f"\nProceeding with {os.path.basename(fasta_file)}:")

    # Step 1: Count k-mers
    with Timer() as total_time:
        kmer_counts = count_kmers(fasta_file, k)
    stats['TIME_SELECTING_KMERS'] = f"{round(total_time.t, 2)} seconds"
    
    # Step 2: Filter solid k-mers
    filtered_kmers = filter_kmers(kmer_counts, threshold)

    # Step 3: Filter build spss
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
    #you can use the test if you need it.. 
    #if not test_fm_index(fm_index, spss, filtered_kmers):
        #print("oh no FM-index validation failed.")
        #exit(1)

    # Step 5: Serialize the FM-index
    dump_file_name = save_fm_index(fm_index, fasta_file, mode, output_file)

    # Optional: Save results for benchmarking
    save_benchmark_results(fasta_file, mode, stats)
    print("-" * 50)  # Prints a line with 50 dashes


