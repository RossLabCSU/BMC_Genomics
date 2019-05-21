    
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

aa_conversion = {'Ala':'A',
                'Cys':'C',
                'Asp':'D',
                'Glu':'E',
                'Phe':'F',
                'Gly':'G',
                'His':'H',
                'Ile':'I',
                'Lys':'K',
                'Leu':'L',
                'Met':'M',
                'Asn':'N',
                'Pro':'P',
                'Gln':'Q',
                'Arg':'R',
                'Ser':'S',
                'Thr':'T',
                'Val':'V',
                'Trp':'W',
                'Tyr':'Y',
                'Ter':'*'}
                
def get_SNV_info(snv_info_file):

    h = open(snv_info_file)
    header = h.readline().rstrip().split('\t')

    #snv_info contains each accession number as keys and a set as the value, 
    #where each item in the set contains: 
    #1) Common Protein Name, 2) Amino acid BEFORE mutation, 3) Mutation position, 4) Amino acid AFTER mutation, 5) Phenotype/disease 
    snv_info = {}
    
    for line in h:
        items = line.rstrip().split('\t')
        mut_type = items[1]
        clin_sig_simple = items[7]  #ClinSigSimple is a simplified representation of ClinicalSignificance. Value=1 if at least one reference classified the variant as "Pathogenic" or "Likely pathogenic". Value=0 otherwise.
        if mut_type == 'single nucleotide variant' and clin_sig_simple == '1':
            name_info = items[2].split('.')
            if ')' not in name_info[-1]:
                continue
            accession = name_info[0]
            mut = name_info[-1][:-1]
            common_name = items[4]
            phenotype = items[13]
            if mut[0:3] not in aa_conversion or mut[-3:] not in aa_conversion:
                continue
            aa_before_mut = aa_conversion[mut[0:3]]
            position = mut[3:-3]
            aa_after_mut = aa_conversion[mut[-3:]]
            
            snv_info[accession] = snv_info.get(accession, set())
            snv_info[accession].add( (common_name, aa_before_mut, position, aa_after_mut, phenotype) )

    h.close()

    return snv_info
    
    
def get_papa_prld_boundaries(papa_results_file):

    papa_prlds = {}
    h = open(papa_results_file)
    header = h.readline()
    for line in h:  
        items = line.rstrip().split('\t')
        accession = items[0]
        
        #Get PrLD boundaries for wild-type proteins only
        if '_' in accession.replace('NM_', ''):
            continue
            
        if len(items) == 5:
            prlds = items[4].replace('"', '').split('_;_')
            for prld in prlds:
                papa_start, papa_end = prld.split(', ')
                papa_start = int(papa_start[1:])
                papa_end = int(papa_end[:-1])
                
                papa_prlds[accession] = papa_prlds.get(accession, [] )
                papa_prlds[accession].append( (papa_start, papa_end) )
        
    return papa_prlds
    
    
def get_plaac_prld_boundaries(plaac_results_file, acc2common):
    """gets PLAAC PrLD boundaries from the PLAAC output file
    Only stores proteins that have a PrLD (i.e. those with a COREscore != NaN in the PLAAC output)
    Each key = protein name
    Each value = (PrLD_start, PrLD_end)
    """
    
    h = open(plaac_results_file)
    
    prld_boundaries = {}
    plaac_positives = []
    
    for line in h:
        if line.startswith('#') or line.startswith('SEQid'):
            continue
            
        items = line.rstrip().split(',')
        if len(items) != 38:
            continue
            
        SEQid, MW, MWstart, MWend, MWlen, LLR, LLRstart, LLRend, LLRlen, NLLR, VITmaxrun, COREscore, COREstart, COREend, CORElen, PRDscore, PRDstart, PRDend, PRDlen, PROTlen, HMMall, HMMvit, COREaa, STARTaa, ENDaa, PRDaa, FInumaa, FImeanhydro, FImeancharge, FImeancombo, FImaxrun, PAPAcombo, PAPAprop, PAPAfi, PAPAllr, PAPAllr2, PAPAcen, PAPAaa = items
        if int(PRDstart) != 0:
            plaac_positives.append(SEQid)
            prld_boundaries[SEQid] = prld_boundaries.get( SEQid, (int(PRDstart), int(PRDend)) )
        
    h.close()

    return prld_boundaries, plaac_positives
    
    
def check_prld_overlap(papa_boundaries, plaac_start, plaac_end):
    toggle = False
    for prld in papa_boundaries:
        papa_start, papa_end = prld
        papa_start = int(papa_start)
        papa_end = int(papa_end)
        if (plaac_start >= papa_start and plaac_start <= papa_end) or (plaac_end >= papa_start and plaac_end <= papa_end) or (plaac_start <= papa_start and plaac_end >= papa_end):
            toggle = True
            
    return toggle
  
def add_annotations(snv_papa_results_file, clinvar_info, acc2common, plaac_prlds):

    output = open('Additional file 3.tsv', 'w')
    output.write('Accession\tCommon Protein Name\tMutation\tHighest PAPA score\tPosition of Highest Scoring Window\tOriginal PAPA score\tMutant PAPAscore - Original PAPAscore\tPhenotype/Disease\tOverlapping PLAAC PrLD?\tPLAAC PrLD Indices\n')

    fig4a_plotting_df = { 'Accessions':[],
                    'Common names':[],
                    'Orig scores':[],
                    'Mutant scores':[] }
                    
    fig4b_plotting_df = { 'Accessions':[],
                    'Common names':[],
                    'Orig scores':[],
                    'Mutant scores':[] }
                    
    h = open(snv_papa_results_file)
    header = h.readline()
    scores = {}
    for line in h:
        if 'protein length below window size' in line:
            continue

        items = line.rstrip().split('\t')
        accession = items[0]
        papa_score = items[1]
        
        if len(items) == 5:
            position = items[2]
            prld_boundaries = []
            prlds = items[4].split('_;_')
            for prld in prlds:
                prld_start, prld_end = prld[1:-1].split(', ')
                prld_boundaries.append( (prld_start, prld_end) )
        else:
            position = '0'
        
        if '_' not in accession.replace('NM_', ''):
            scores[accession] = scores.get(accession, [ (papa_score,position) ])

        else:
            temp = accession.split('_')
            acc = temp[0] + '_' + temp[1]
            mut = temp[2]
            
            #SKIP TRUNCATION MUTATIONS
            if '*' in mut:
                continue

            #delta_papa_score represents mutant PAPA score minus original PAPA score (positive values indicate an increase in prion propensity with the given mutation)
            delta_PAPA_score = float(papa_score) - float(scores[acc][0][0])
            scores[acc].append( (accession, papa_score, position, str(delta_PAPA_score) ) )
            phenotype = ''
            for mutant in clinvar_info[acc]:
                if mut == (mutant[1] + mutant[2] + mutant[3]):
                    phenotype = mutant[4]
                    
            orig_PAPAscore = float(papa_score) - float(delta_PAPA_score)
            
            if acc in acc2common:
                common_name = acc2common[acc]
            else:
                common_name = acc
                
            if orig_PAPAscore < 0.05 and float(papa_score) > 0.05:
                fig4a_plotting_df['Accessions'].append(acc)
                fig4a_plotting_df['Common names'].append(common_name)
                fig4a_plotting_df['Orig scores'].append(orig_PAPAscore)
                fig4a_plotting_df['Mutant scores'].append( float(papa_score) )
                
            if orig_PAPAscore > 0.05 and float(papa_score) > orig_PAPAscore:
                fig4b_plotting_df['Accessions'].append(acc)
                fig4b_plotting_df['Common names'].append(common_name)
                fig4b_plotting_df['Orig scores'].append(orig_PAPAscore)
                fig4b_plotting_df['Mutant scores'].append( float(papa_score) )

                    
            output.write(accession + '\t')
            output.write(common_name + '\t' + mut + '\t' + papa_score + '\t' + position + '\t' + str(orig_PAPAscore) + '\t' + str(delta_PAPA_score) + '\t' + phenotype + '\t')
            
            if acc in plaac_prlds:
                plaac_start, plaac_end = plaac_prlds[acc]
                is_overlapping = check_prld_overlap(prld_boundaries, plaac_start, plaac_end)
                if is_overlapping:
                    output.write('1\t' + str(plaac_prlds[acc]) + '\n')
                else:
                    output.write('0\t' + str(plaac_prlds[acc]) + '\n')
            else:
                output.write('0\t-\n')
                
    output.close()
    
    return fig4a_plotting_df, fig4b_plotting_df
    
    
def plot_fig4a(small_df, plaac_positives):
            
    #get co-sorted lists
    current_iso_PAPAscore, max_iso_PAPAscore, common_name, accession = zip(*sorted(zip(small_df['Orig scores'], small_df['Mutant scores'], small_df['Common names'], small_df['Accessions'])))
    
    #plot barchart
    palette = sns.color_palette()
    for i in range(len(small_df['Common names'])):
        plt.bar(i, current_iso_PAPAscore[i], color=palette[0])
        if current_iso_PAPAscore[i] < 0:
            b = 0
            plt.bar(i, max_iso_PAPAscore[i], bottom=b, color=palette[1])
        else:
            b = current_iso_PAPAscore[i]
            plt.bar(i, (max_iso_PAPAscore[i] - current_iso_PAPAscore[i]), bottom=b, color=palette[1])
        if accession[i] in plaac_positives:
            plt.text(i, (max_iso_PAPAscore[i]), '*', fontsize=18, ha='center', va='bottom')

        ax = plt.gca()
        ax.set_xticklabels(common_name)
        plt.xticks(range(len(common_name)), rotation=90, fontsize=12, fontname='Arial')
    plt.legend(['Wild-type', 'Mutation'], loc='lower right', prop={'family':'Arial', 'size':20})
    ax = plt.gca()
    ax.set_xticklabels(common_name, fontname='Arial')
    plt.yticks(fontsize=12, fontname='Arial')
    plt.xlabel('Proteins', fontname='Arial', fontsize=14)
    plt.ylabel('PAPA Scores', fontname='Arial', fontsize=14)
    plt.tight_layout()
    plt.plot([0.05 for i in range(len(common_name))], color='0.2', linestyle='--')
    fig = plt.gcf()
    fig.set_size_inches(15, 7.5)
    plt.savefig('Fig 4A.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    

def plot_fig4b(small_df, plaac_positives):
    
    #get co-sorted lists
    current_iso_PAPAscore, max_iso_PAPAscore, common_name, accession = zip(*sorted(zip(small_df['Orig scores'], small_df['Mutant scores'], small_df['Common names'], small_df['Accessions'])))
    
    #plot barchart
    palette = sns.color_palette()
    for i in range(len(small_df['Common names'])):
        plt.bar(i, current_iso_PAPAscore[i], color=palette[0])
        if current_iso_PAPAscore[i] < 0:
            b = 0
            plt.bar(i, max_iso_PAPAscore[i], bottom=b, color=palette[1])
        else:
            b = current_iso_PAPAscore[i]
            plt.bar(i, (max_iso_PAPAscore[i] - current_iso_PAPAscore[i]), bottom=b, color=palette[1])
        if accession[i] in plaac_positives:
            plt.text(i, (max_iso_PAPAscore[i]), '*', ha='center', va='bottom')
        ax = plt.gca()
        ax.set_xticklabels(common_name)
        plt.xticks(range(len(common_name)), rotation=90, fontsize=10, fontname='Arial')
    plt.legend(['Wild-type', 'Mutation'], loc='upper left', prop={'family':'Arial', 'size':12})
    ax = plt.gca()
    ax.set_xticklabels(common_name, fontname='Arial')
    plt.yticks(fontsize=10, fontname='Arial')
    plt.xticks(fontsize=6, fontname='Arial')
    plt.xlabel('Proteins', fontname='Arial', fontsize=12)
    plt.ylabel('PAPA Scores', fontname='Arial', fontsize=12)
    plt.tight_layout()
    plt.plot([0.05 for i in range(len(common_name))], color='0.2', linestyle='--')
    fig = plt.gcf()
    fig.set_size_inches(10, 5)
    plt.savefig('Fig 4B.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def main(snv_papa_results_file, snv_plaac_results_file, snv_info_file):
    
    #IMPORT ACCESSION TO COMMON NAME CONVERSION DICTIONARY
    acc2common = pickle.load(open('Human_prots_GenBankAccession_to_Common_GeneName.dat', 'rb'))
    
    #GET SNV INFO
    snv_info = get_SNV_info(snv_info_file)
    
    #GET PLAAC PrLDs
    plaac_boundaries, plaac_positives = get_plaac_prld_boundaries(snv_plaac_results_file, acc2common)
    
    #ADD ANNOTATIONS AND OUTPUT ADDITIONAL FILE 3
    fig4a_plotting_df, fig4b_plotting_df = add_annotations(snv_papa_results_file, snv_info, acc2common, plaac_boundaries)
    
    #PLOTTING FOR FIGURE 3
    plot_fig4a(fig4a_plotting_df, plaac_positives)
    plot_fig4b(fig4b_plotting_df, plaac_positives)
    
        
if __name__ == '__main__':
    import sys
    #snv_papa_results_file --> 'Human_proteome_with_SNVs_PAPA_results.tsv' (the result of running mPAPA.py on Human_proteome_with_SNVs.FASTA with threshold=0.0)
    #snv_plaac_results_file --> 'Human_Protein_Seqs_ActiveDriverDB_PLAAC_results.csv' (the result of running plaac.jar on 'Human_Protein_Seqs_ActiveDriverDB.FASTA' with alpha=0 and the FASTA file used to calculate background amino acid frequencies)
    #snv_info_file --> 'variant_summary.txt' (uncompressed after download)
    snv_papa_results_file, snv_plaac_results_file, snv_info_file = sys.argv[1:]
    main(snv_papa_results_file, snv_plaac_results_file, snv_info_file)
    
