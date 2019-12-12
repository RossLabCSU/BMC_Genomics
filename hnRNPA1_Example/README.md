
# Reproducing hnRNPA1 Plotting

The following commands will reproduce the PTM mapping and plotting for hnRNPA1 and FUS in Fig 6A and Fig 6B.

### Dependencies and Required Files
Scripts used in these analyses rely on the following external dependencies:\
Biopython\
Matplotlib\
Seaborn\
Numpy\
Additional file 4 (from PTM_Analysis pipeline)

### Instructions
1. Download mPAPA.py
2. Download all files in BMC_Genomics/hnRNPA1_and_FUS_Examples
3. Copy "Additional file 4.tsv" (from PTM_Analysis pipeline) into the folder containing the hnRNPA1 scripts.
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

>\>python mPAPA.py --threshold 0.0 --verbose --ignore_fold_index -o HNRNPA1_All_Variants_PAPA_results_verbose.tsv HNRNPA1_All_variants.FASTA

>\>python mPAPA.py --threshold 0.0 --verbose --ignore_fold_index -o FUS_All_Variants_PAPA_results_verbose.tsv FUS_All_variants.FASTA

>\>python hnRNPA1_AND_FUS_EXAMPLE_ANALYSIS_PIPELINE.py "Additional file 4.tsv"
