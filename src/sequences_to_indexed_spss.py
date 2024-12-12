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

def prepare_output_dir(base_name, mode, k, t, extension="dump", parent_dir="benchmark"):
    """
    Prepare directory for outputs
    """
    output_dir = os.path.abspath(os.path.join(parent_dir, base_name))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    file_name = f"fmi_{base_name}_{mode}_k{k}_t{t}.{extension}"
    return os.path.join(output_dir, file_name)


def save_fm_index(fm_index, fasta_file, mode, k, t, output_file=None):
    """
    Serializes and saves the FM-index to a file. Automatically creates folders if needed.
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
    """
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    parent_dir = os.path.abspath("benchmark/stats")  # Répertoire unique pour les stats

    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)  # Crée le dossier "stats" s'il n'existe pas

    # Nom du fichier CSV basé uniquement sur le mode et le dataset
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

    # Écriture dans le fichier CSV
    with open(csv_file_name, mode='a', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        if csv_file.tell() == 0:  # Si le fichier est vide, écris l'en-tête
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
    if output_file and output_file.endswith(".dump"):
        dump_file_name = save_fm_index(fm_index, fasta_file, mode, k, threshold, output_file=output_file)
    else:
        dump_file_name = save_fm_index(fm_index, fasta_file, mode, k, threshold, output_file=None)
        save_benchmark_results(fasta_file, mode, stats)

    for key, value in stats.items():
        print(f"  {key}: {value}")
    print("-" * 50)  # Prints a line with 50 dashes


