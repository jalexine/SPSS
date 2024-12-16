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
    # Open the file in binary read mode and load the FM-index object using pickle
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
        # Generate the suffix array using the Karkkainen-Sanders algorithm
        self.suffix_array = tks.simple_kark_sort(sequence)
        self.sequence = sequence
        # Compute the Burrows-Wheeler Transform (BWT) for the sequence
        self.bwt = self.set_bwt()
        # Compute the necessary mapping and rank arrays for the FM-index
        self.n, self.rank = self.set_n_and_ranks()

    def save(self, filename):
        """
        Serialize the current FM-index instance and save it to a file.

        Parameters:
            filename (str): The path to the file where the FM-index should be saved.
        
        Returns:
            None
        """
        # Open the file in binary write mode and serialize the current FM-index object
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    def set_bwt(self):
        """
        Compute the Burrows–Wheeler Transform (BWT) of the indexed sequence.

        The BWT is derived using the suffix array, 
        taking the character preceding each suffix.

        Parameters:
            None

        Returns:
            str: The BWT string of the indexed sequence.
        """
        # Build the BWT by iterating over the suffix array and getting the preceding character
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
        n = {}
        rank = {}
        cumulative_count = {}
        # Sort the BWT string and prepare for ranking

        sorted_bwt = sorted(self.bwt)
        # For each character in the BWT, initialize the rank and n values
        for char in set(self.bwt):
            n[char] = 0
            rank[char] = [0] * (len(self.bwt) + 1)  # rank will store cumulative counts
            cumulative_count[char] = 0
        # Update 'n' and 'rank' based on the BWT
        for i, char in enumerate(self.bwt):
            if cumulative_count[char] == 0:
                n[char] = sorted_bwt.index(char)
        # Increment the occurrence count of char
            cumulative_count[char] += 1
            for c in rank:
                rank[c][i + 1] = cumulative_count[c]
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
        # Find the LF-mapped index using the starting position 'n' and the occurrence 'k'
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
        # Iterate through the BWT starting from index 'l' to find the next occurrence of 'alpha'
        for i in range(l, len(self.bwt)):
            if self.bwt[i] == alpha:
                return i # Return the index of the next occurrence of 'alpha'
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
    # Search backwards in the BWT for the character 'alpha'
        for i in range(l, -1, -1):
            if self.bwt[i] == alpha:
                return i # Return index if 'alpha' is found
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
            return 0 # Return 0 if the character 'n' is not found in the rank
        return self.rank[n][i + 1] # Return the cumulative count for character 'n' 

    def contains(self, q):
        """
        Check if a given query substring q occurs in the indexed sequence.

        The query is processed in reverse using backward search on the BWT and cumulative counts.

        Parameters:
            q (str): The query substring to search for.
        
        Returns:
            bool: True if q is found in the indexed sequence, False otherwise.
        """
        # Initialize the search range for the BWT
        l, r = 0, len(self.bwt) - 1
        # Process the query in reverse order
        for char in reversed(q):
            # If the character is not in the FM-index, return False
            if char not in self.n:
                return False
            # Update the lower bound of the range based on the character's occurrence
            l = self.n[char] + (self.occ(char, l - 1) if l > 0 else 0)
            # Update the upper bound of the range based on the character's occurrence
            r = self.n[char] + self.occ(char, r) - 1
            # If the lower bound exceeds the upper bound = substring not found
            if l > r:
                return False
        return l <= r