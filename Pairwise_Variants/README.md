
# Reproducing neXtProt Pairwise Variant Analysis

The following commands will reproduce all figures and additional files pertaining to the neXtProt pairwise variants within high-scoring human PrLDs (Fig 2A-D, and Additional file 1).

### Dependencies
Scripts used in these analyses rely on the following external dependencies:\
Biopython\
Matplotlib\
Seaborn\
Numpy\
Pandas

### Instructions
1. Download mPAPA.py
2. Download all files in CELL_REPORTS/Pairwise_Variants
3. Extract "nextprot_all.peff" from compressed file
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python GENERATE_NEXTPROT_FASTA.py nextprot_all.peff

>\>python mPAPA.py --threshold 0.0 -o NextProt_proteome_mPAPA_results.tsv NextProt_protein_sequences.FASTA

>\>python GENERATE_NEXTPROT_PAIRWISE_VARIANT_SEQS.py NextProt_proteome_mPAPA_results.tsv nextprot_all.peff

>\>python mPAPA.py --threshold 0.0 -o NextProt_PairwiseVariant_ProteinSeqs_mPAPA_results.tsv NextProt_PairwiseVariant_ProteinSeqs.FASTA

>\>python GET_NEXTPROT_PAIRWISE_VARIANT_SCORES.PY NextProt_proteome_mPAPA_results.tsv NextProt_PairwiseVariant_ProteinSeqs_mPAPA_results.tsv

NOTE: Generating pairwise variants and running mPAPA for *all* pairwise variants may require parallel computing, but representative distributions of PAPA scores and PAPA score ranges depicted in Fig 2 can be obtained using only a subset of variants.
