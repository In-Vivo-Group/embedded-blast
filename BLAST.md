# Understanding Current BLAST: Literature Review and Stage Analysis

## 1. BLAST Literature Review

### Key Papers:

1. Altschul, S.F., et al. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403-410.
   - This seminal paper introduced the BLAST algorithm.

2. Altschul, S.F., et al. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Research, 25(17), 3389-3402.
   - Introduced improvements including gapped alignments and position-specific scoring matrices.

3. Camacho, C., et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10, 421.
   - Describes the reimplementation of BLAST with a focus on modularity and performance.

### Key Concepts from Literature:

1. Sequence Similarity: BLAST is fundamentally about finding similar sequences in a database.
2. Local Alignment: Unlike global alignment algorithms, BLAST focuses on finding local regions of similarity.
3. Heuristic Approach: BLAST sacrifices some sensitivity for speed, making it suitable for large-scale searches.
4. Statistical Significance: BLAST provides E-values to indicate the statistical significance of matches.

## 2. BLAST Algorithm Stages

### Stage 1: Seeding

- Purpose: Quickly identify potential matches.
- Process:
  1. Break query sequence into words (typically 3 amino acids or 11 nucleotides).
  2. Generate list of similar words based on a scoring matrix.
  3. Search database for exact matches to these words.

### Stage 2: Extension

- Purpose: Extend seed matches to create local alignments.
- Process:
  1. Extend matches in both directions without allowing gaps.
  2. Continue extension until the score falls below a threshold.
  3. For nucleotide BLAST, this creates High-scoring Segment Pairs (HSPs).

### Stage 3: Evaluation

- Purpose: Assess statistical significance of alignments.
- Process:
  1. Calculate bit score for each HSP.
  2. Compute E-value (expected number of HSPs with this score by chance).
  3. Filter alignments based on E-value threshold.

### Stage 4: Gapped Alignment (for protein BLAST and newer versions)

- Purpose: Refine alignments by allowing gaps.
- Process:
  1. Use dynamic programming to introduce gaps in promising alignments.
  2. Recalculate scores and E-values for gapped alignments.

### Stage 5: Trace-back and Display

- Purpose: Generate human-readable alignments.
- Process:
  1. Perform trace-back to reconstruct full alignments.
  2. Format alignments for display, showing matches, mismatches, and gaps.

## 3. Key Algorithmic Components

1. Scoring Matrices:
   - BLOSUM62 (default for protein BLAST)
   - PAM matrices (alternative for protein sequences)
   - Simple match/mismatch scores for nucleotide BLAST

2. Word Size:
   - Protein BLAST typically uses 3-letter words
   - Nucleotide BLAST typically uses 11-letter words

3. E-value Calculation:
   - Based on extreme value distribution
   - Considers database size and query length

4. Filters:
   - Low-complexity filters (e.g., SEG for proteins, DUST for nucleotides)
   - Mask repetitive elements to avoid spurious matches

## 4. Variants and Improvements

1. PSI-BLAST (Position-Specific Iterated BLAST):
   - Uses position-specific scoring matrices
   - Iterative process to detect remote homologs

2. PHI-BLAST (Pattern-Hit Initiated BLAST):
   - Combines pattern matching with local alignment

3. DELTA-BLAST (Domain Enhanced Lookup Time Accelerated BLAST):
   - Uses conserved domain information to improve sensitivity

4. MegaBLAST:
   - Optimized for highly similar sequences
   - Uses larger word size (28 by default)

## 5. Computational Considerations

1. Time Complexity:
   - Roughly linear with database size
   - Significant speedup compared to dynamic programming approaches

2. Space Complexity:
   - Main memory usage is for index structures
   - Can be optimized for large databases using disk-based approaches

3. Parallelization:
   - Embarrassingly parallel across different queries
   - Can be parallelized for single queries using techniques like GPU acceleration

## Conclusion

Understanding these aspects of the current BLAST algorithm provides a solid foundation for developing a new method using embedding space. Key areas to focus on for improvement include:

1. Seeding: How can embedding space be used to identify potential matches more efficiently?
2. Extension: Can embeddings provide a more nuanced way to extend matches beyond simple linear extension?
3. Scoring: How can similarity in embedding space be translated into biologically meaningful scores?
4. Statistical Significance: What new statistical models might be needed for embedding-based alignments?

By addressing these questions, your new method has the potential to significantly advance the field of sequence similarity search.