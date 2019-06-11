
# Reproducing Isoform Comparison Analysis

The following commands will reproduce all figures and additional files pertaining to the comparison of "high-confidence" human protein isoforms containing PrLDs from ActiveDriverDB (Fig 3A, Fig 3B, and Additional file 2).

### Dependencies
Scripts used in these analyses rely on the following external dependencies:\
Biopython\
Matplotlib\
Seaborn

### Instructions
1. Download mPAPA.py
2. Download all files in CELL_REPORTS/Isoform_Comparison
3. Navigate via command line to appropriate folder containing downloaded files.
4. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python mPAPA.py --threshold 0.0 -o Human_Protein_Seqs_ActiveDriverDB_PAPA_results.tsv Human_Protein_Seqs_ActiveDriverDB.FASTA

>\>python HUMAN_ISOFORM_PAPA_ANALYSIS_PIPELINE.py Human_Protein_Seqs_ActiveDriverDB_PAPA_results.tsv Human_Protein_Seqs_ActiveDriverDB_PLAAC_results.csv Human_Protein_Seqs_ActiveDriverDB.FASTA
