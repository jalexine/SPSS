#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FM-index implementation for genome indexing and querying.
"""

import pickle
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../pysuffix3"))
from pysuffix3 import tools_karkkainen_sanders as tks

def load_fm_index(filename):
    """
    Load an FM-index object from a given file.

    This function deserializes and returns a precomputed FM-index from the specified file.

    Parameters:
        filename (str): The path to the file containing the serialized FM-index object.

    Returns:
        FmIndex: The deserialized FM-index instance loaded from the file.
    """
    with open(filename, 'rb') as f:
        fm_index = pickle.load(f)
    return fm_index

class FmIndex:
    """
    A class for constructing and querying an FM-index of a given sequence.

    The FM-index allows efficient substring queries (contains) on the indexed sequence.
    """

    def __init__(self, sequence):
        """
        Initialize the FM-index from the provided sequence.

        Parameters:
            sequence (str): The input sequence (e.g., a DNA sequence or any string).
        
        Returns:
            None
        """
        self.suffix_array = tks.simple_kark_sort(sequence)
        self.sequence = sequence
        self.bwt = self.set_bwt()
        self.n, self.rank = self.set_n_and_ranks()

    def save(self, filename):
        """
        Serialize the current FM-index instance and save it to a file.

        Parameters:
            filename (str): The path to the file where the FM-index should be saved.
        
        Returns:
            None
        """
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    def set_bwt(self):
        """
        Compute the Burrows–Wheeler Transform (BWT) of the indexed sequence.

        The BWT is derived using the suffix array, taking the character preceding each suffix.

        Parameters:
            None

        Returns:
            str: The BWT string of the indexed sequence.
        """
        bwt = "".join(self.sequence[i - 1] if i > 0 else self.sequence[-1] for i in self.suffix_array)
        return bwt

    def set_n_and_ranks(self):
        """
        Compute the character mapping and rank arrays needed for the FM-index.

        The 'n' dictionary maps each character to its starting position in the BWT's sorted array.
        The 'rank' dictionary provides cumulative counts of each character at each position in the BWT.

        Parameters:
            None
        
        Returns:
            tuple:
                dict: A dictionary `n` where n[char] gives the starting index of char in the sorted BWT.
                dict: A dictionary `rank` where rank[char] is a list of cumulative counts of char up to each position.
        """
        n, rank = {}, {}
        for char in set(self.bwt):
            # Position of char in the sorted BWT (number of chars before char in sorted order)
            n[char] = len(''.join(sorted(self.bwt)).split(char)[0])
            rank[char] = [0] * (len(self.bwt) + 1)

        # Build cumulative ranks
        for i, char in enumerate(self.bwt):
            for c in rank:
                rank[c][i + 1] = rank[c][i] + (1 if c == char else 0)
        return n, rank

    def lf(self, alpha, k):
        """
        Compute the LF mapping for a given character and occurrence number.

        LF(alpha, k) finds the position in the suffix array that corresponds to the k-th occurrence
        of character alpha in the BWT.

        Parameters:
            alpha (str): The character to map.
            k (int): The occurrence count (1-based) of alpha.
        
        Returns:
            int: The LF-mapped index in the suffix array for the given character occurrence.
        """
        lf_index = self.n[alpha] + k - 1
        return lf_index

    def find_next(self, alpha, l):
        """
        Find the next occurrence of a given character in the BWT starting at index l or after.

        Parameters:
            alpha (str): The character to search for in the BWT.
            l (int): The starting index for the search.
        
        Returns:
            int: The index of the next occurrence of alpha at or after l, or -1 if not found.
        """
        for i in range(l, len(self.bwt)):
            if self.bwt[i] == alpha:
                return i
        return -1

    def find_prev(self, alpha, l):
        """
        Find the previous occurrence of a given character in the BWT up to index l.

        Parameters:
            alpha (str): The character to search for in the BWT.
            l (int): The ending index for the backward search.
        
        Returns:
            int: The index of the occurrence of alpha at or before l, or -1 if not found.
        """
        for i in range(l, -1, -1):
            if self.bwt[i] == alpha:
                return i
        return -1
    
    def occ(self, n, i):
        """
        Count the occurrences of a specific character up to a given position in the BWT.

        Parameters:
            n (str): The character for which occurrences are counted.
            i (int): The end index (0-based) up to which occurrences are counted.
        
        Returns:
            int: The number of times character n appears in self.bwt up to index i.
        """
        if n not in self.rank:
            return 0
        return self.rank[n][i + 1]

    def contains(self, q):
        """
        Check if a given query substring q occurs in the indexed sequence.

        The query is processed in reverse using backward search on the BWT and cumulative counts.

        Parameters:
            q (str): The query substring to search for.
        
        Returns:
            bool: True if q is found in the indexed sequence, False otherwise.
        """
        l, r = 0, len(self.bwt) - 1
        for char in reversed(q):
            if char not in self.n:
                return False
            l = self.n[char] + (self.occ(char, l - 1) if l > 0 else 0)
            r = self.n[char] + self.occ(char, r) - 1
            if l > r:
                return False
        return l <= r