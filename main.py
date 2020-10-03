# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 13:58:11 2020

@author: Kenneth
"""

"""
Naive exact matching algorithm with performance tracking.
"""
def naive_and_counting(p, t):
    occurrences = []
    N_alignment = 0
    N_comparison = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        N_alignment += 1
        match = True
        for j in range(len(p)):  # loop over characters
            N_comparison += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, N_alignment, N_comparison

"""
This module provides the BoyerMoore class, which encapsulates the preprocessing info used by the boyer_moore function.
It preprocesses the pattern P into the tables needed to execute the bad character and good suffix rules.
"""
from bm_preproc import *

"""
Boyer-Moore algorithm with performance tracking.
"""
def boyer_moore_and_count(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    N_alignment = 0
    N_comparison = 0    
    while i < len(t) - len(p) + 1:
        N_alignment += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            N_comparison += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, N_alignment, N_comparison

"""
This module builds a k-mer index.
"""
from kmer_index import *

"""
This function finds approximate matches by taking advantage of a k-mer index and implementing the pigeonhole principle.
"""
def approximate_match_index_pigeon(p, t, n, k):    
    index = Index(t, k)
    segment_length = round(len(p)/(n+1))
    matches = []
    index_hits = 0
    for i in range(n+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))        
        segment = p[start:end]
        hits = index.query(segment)
        if len(hits) > 0:
            index_hits += len(hits)
            for j in hits:
                if j < start or j-start+len(p) > len(t):
                    continue
                mismatch = 0
                for dummy in range(0, start):
                    if p[dummy] != t[j-start+dummy]:
                        mismatch += 1
                        if mismatch > n:
                            break
                for dummy in range(end, len(p)):
                    if p[dummy] != t[j-start+dummy]:
                        mismatch += 1
                        if mismatch > n:
                            break                
                if mismatch <= n:
                    matches.append(j - i*segment_length)
        matches = list(set(matches))
    return matches, index_hits

"""
This module builds an index that handles subsequences that take every lth character.
"""
from SubseqIndex import *

"""
This function finds approximate matches by taking advantage of a general k-mer index and implementing the pigeonhole principle.
It consider subsequences that take every lth character.
"""
def approximate_match_general_index_pigeon(p, t, n, k, l):    
    index = SubseqIndex(t, k, l)
    matches = []
    index_hits = 0
    for i in range(n+1):     
        segment = p[i:]
        hits = index.query(segment)
        if len(hits) > 0:
            index_hits += len(hits)
            for j in hits:
                if j-i < 0 or j-i+len(p) > len(t):
                    continue
                mismatch = 0
                for dummy in range(0, len(p)):
                    if p[dummy] != t[j-i+dummy]:
                        mismatch += 1
                        if mismatch > n:
                            break              
                if mismatch <= n:
                    matches.append(j-i)
        matches = list(set(matches))
    return matches, index_hits

"""
Read the provided FASTA file.
"""
file = open('chr1.GRCh38.excerpt.fasta', 'r')
seq = ''
for line in file:
    if line[0] != '>':
        dummy = line.rstrip()
        seq += dummy
        
"""
Question 1.
1. Align p with seq using the naive exact matching algorithm.
2. Identify all the matches.
3. Measure the performance by calculating the alignment and comparison counts.
"""
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
O, A, C = naive_and_counting(p, seq)
print(A)

"""
Question 2.
"""
print(C)

"""
Question 3.
1. Align p with seq using the Boyer-Moore algorithm.
2. Identify all the matches.
3. Measure the performance by calculating the alignment and comparison counts.
"""
uppercase_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ '
p_bm = BoyerMoore(p, uppercase_alphabet)
O, A, C = boyer_moore_and_count(p, p_bm, seq)
print(A)

"""
Question 4.
1. Build a k-mer index.
2. Apply the pigeonhole principle by dividing p into n+1 segments.
3. Find all the index hits for each segment.
4. For each index hit, verify that it is a match.
5. Identify the unique set of matches and the total number of index hits for all segments.
"""
n = 2
k = 8
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
matches, index_hits = approximate_match_index_pigeon(p, seq, n, k)
print(len(matches))

"""
Question 5.
"""
print(index_hits)

"""
Question 6.
Same as question 4, but use a generalised index which considers subsequences that take every lth character.
Also, instead of dividing p into n+1 segments, consider n+1 reading frames.
"""
n = 2
k = 8
l = 3
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
matches, index_hits = approximate_match_general_index_pigeon(p, seq, n, k, l)
print(index_hits)