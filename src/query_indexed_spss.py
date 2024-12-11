"""
Given a set of queries, check using a FM-index if the query is in a set of kmers represented by their SPSS
"""

import argparse
import fmi
from sequences_to_indexed_spss import reverse_complement


def shared_kmers(
    sequence: str,
    my_fmi: fmi.FmIndex,
    k: int,
) -> float:
    """Returns the number of shared kmers between a sequence and a fmi

    Args:
        sequence (str): sequence
        fmi (fmi.Fmi): fmi
        k (int): k

    Returns:
        float: number of shared kmers between a sequence and a fmi, rounded to 3 decimals
    """
    nb_kmers = len(sequence) - k + 1
    nb_shared_kmers = 0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        if my_fmi.contains(kmer):
            nb_shared_kmers += 1
        else:
            kmer = reverse_complement(kmer)
            if my_fmi.contains(kmer):
                nb_shared_kmers += 1
    return round(nb_shared_kmers / nb_kmers, 3)


def parse_input_fasta_file(
    fasta_file_name: str, k: int, my_fmi: fmi.FmIndex, output_file_name: str
):
    """
    For each sequence in the fasta file, calculate the number of shared kmers with the fmi
    """
    with open(output_file_name, encoding="utf-8", mode="w") as output_stream, open(
        fasta_file_name, encoding="utf-8"
    ) as fasta_file:
        while True:
            line = fasta_file.readline()
            if not line:
                break
            if not line[0] == ">":
                print(f"Fasta format error, line {line}")
            comment = line.strip(">").strip()
            line = fasta_file.readline()
            line = line.strip()
            ratio = shared_kmers(line, my_fmi, k)
            output_stream.write(f"{comment}\t{ratio}\n")


def main():
    """
    Main function
    QueryIndexedSPSS.py -q query_file_name.fa -i fmi_index_file_name -k kmer_size -o output_file [-h]
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("-q", "--query", help="query file name", required=True)
    parser.add_argument("-i", "--index", help="index file name", required=True)
    parser.add_argument("-k", "--kmer", help="kmer size", type=int, required=True)
    parser.add_argument("-o", "--output", help="output file name", required=True)
    args = parser.parse_args()
    fmi_index = fmi.load_fm_index(args.index)
    parse_input_fasta_file(args.query, args.kmer, fmi_index, args.output)


if __name__ == "__main__":
    main()
