
#README for Reproducing ClinVar Disease-Associated Mutant Analysis

##Analysis pipeline for ClinVar SNVs
The following commands will reproduce all figures and additional files pertaining to the ClinVar disease-associated mutations within human PrLDs (Fig 4A, Fig 4B, and Additional file 3).

###Dependencies
Scripts used in these analyses rely on the following external dependencies:
Biopython
Matplotlib
Seaborn

###Instructions
1. Download mPAPA.py
2. Download all files in GENOME_BIOLOGY/Disease_Associated_Mutants
3. Extract 'variant_summary.txt' from compressed file
4. Navigate to appropriate folder via command line.
4. Run the following commands in sequence:

>\>python GENERATE_ClinVar_SNV_FASTA.py Human_Protein_Seqs_ActiveDriverDB.FASTA variant_summary.txt

>\>python mPAPA.py --threshold 0.0 -o Human_proteome_with_SNVs_PAPA_results.tsv Human_proteome_with_SNVs.FASTA

>\>python HUMAN_ClinVar_SNV_ANALYSIS_PIPELINE.py Human_proteome_with_SNVs_PAPA_results.tsv Human_Protein_Seqs_ActiveDriverDB_PLAAC_results.csv variant_summary.txt
