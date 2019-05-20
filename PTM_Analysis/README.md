
# Reproducing PTM Analysis

The following commands will reproduce all figures and additional files pertaining to the analysis of PTMs within human PrLDs (Fig 5A, Fig 5B, Additional file 4, and Additional file 5).

### Dependencies
Scripts used in these analyses rely on the following external dependencies:\
Biopython\
Numpy\
Scipy\
Statsmodels\
Matplotlib\
Seaborn

### Instructions
1. Download mPAPA.py
2. Download all files in GENOME_BIOLOGY/PTM_Analysis
3. Navigate to appropriate folder via command line.
4. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python mPAPA.py --threshold 0.0 -o Human_Protein_Seqs_ActiveDriverDB_PAPA_results.tsv Human_Protein_Seqs_ActiveDriverDB.FASTA

>\>python HUMAN_PTM_ANALYSIS_PIPELINE.py Human_Protein_Seqs_ActiveDriverDB.FASTA site_table.tsv Human_Protein_Seqs_ActiveDriverDB_PAPA_results.tsv Human_Protein_Seqs_ActiveDriverDB_PLAAC_results.csv
