# Boyer_Moore_Index_Pigeonhole
Context: By using and modifying Python programs provided in the Coursera course entitled 'Algorithms for DNA Sequencing', I aligned DNA sequencing reads with an excerpt of human chromosome 1. I used the naive exact matching algorithm, the Boyer-Moore algorithm, two k-mer indices (substrings and subsequences that consider every lth character), and the pigeonhole principle.

About:
1. naive_and_counting() aligns the pattern 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG' with the excerpt using the naive exact matching algorithm, identifies every match, and measures the algorithm's performance by counting instances of alignment and comparison.
2. boyer_moore_and_count() does the same with the Boyer-Moore algorithm.
3. approximate_match_index_pigeon() builds a k-mer index of substrings in the excerpt, applies the pigeonhole principle by dividing the pattern 'GGCGCGGTGGCTCACGCCTGTAAT' into n+1 segments, finds all the index hits for each segment, verifies each hit as to whether it is a match, and identifies the unique set of matches and the total number of index hits for all segments.
4. approximate_match_general_index_pigeon() does the same by building and using a k-mer index of subsequences that consider every lth character. Instead of dividing p into n+1 segments, it considers n+1 reading frames.

Files:
1. main.py, SubseqIndex.py, kmer_index.py, bm_preproc.py, and chr1.GRCh38.excerpt.fa must be in the same directory.
2. main.py should be implemented in Python 3.7.
3. chr1.GRCh38.excerpt.fa is a FASTA file containing an excerpt of human chromosome 1.

Modules:
1. bisect.
2. unittest.
