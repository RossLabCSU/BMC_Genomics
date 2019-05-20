
import pickle
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def compare_isoforms(papa_scores, papa_prlds, plaac_prlds, high_score_positions, seqs):

    output = open('Additional file 2.tsv', 'w')
    output.write('Accession\tCommon Protein Name\tLength of CURRENT Isoform\tPAPA score of CURRENT isoform\tPosition of Highest Scoring Window\tLongest Isoform Length\tHighest-Scoring Isoform PAPA score\tHighest-scoring Isoform PAPAscore - CURRENT Isoform PAPAscore\tOverlapping PLAAC PrLD?\tPLAAC PrLD Indices\n')

    fig3a_plotting_df = {'Common_name':[],
        'Min PAPAscore':[],
        'Max PAPAscore':[]}
        
    fig3b_plotting_df = {'Common_name':[],
        'Min PAPAscore':[],
        'Max PAPAscore':[]}
        
    for protein in papa_scores:
        isoform_accessions, isoform_lengths, isoform_max_scores, isoform_max_score_positions, delta_papascores, highest_iso_score = calculate_deltaPAPA(protein, papa_scores, seqs, high_score_positions)
        if max(isoform_max_scores) > 0.05 and min(isoform_max_scores) < 0.05:
            fig3a_plotting_df['Common_name'].append(protein)
            fig3a_plotting_df['Min PAPAscore'].append( min(isoform_max_scores) )
            fig3a_plotting_df['Max PAPAscore'].append( max(isoform_max_scores) )
            
        high_scores = set()
        for score in isoform_max_scores:
            if score > 0.05:
                high_scores.add(score)
                
        if len(high_scores) > 1:
            fig3b_plotting_df['Common_name'].append(protein)
            fig3b_plotting_df['Min PAPAscore'].append( min(high_scores) )
            fig3b_plotting_df['Max PAPAscore'].append( max(high_scores) )
        
        for i in range( len(isoform_accessions) ):
            accession = isoform_accessions[i]
            output.write(accession + '\t' + protein + '\t' + str( len( seqs[isoform_accessions[i]] )) + '\t' + str(isoform_max_scores[i]) + '\t' + str(isoform_max_score_positions[i]) + '\t' + str( max(isoform_lengths) ) + '\t' + str(highest_iso_score) + '\t' + str( delta_papascores[i] ) + '\t')

            if accession in plaac_prlds:
                plaac_start, plaac_end = plaac_prlds[accession]
                if accession in papa_prlds:
                    is_overlapping = check_prld_overlap(papa_prlds[accession], plaac_start, plaac_end)
                else:
                    is_overlapping = False
                if is_overlapping:
                    output.write('1\t' + str(plaac_prlds[accession]) + '\n')
                else:
                    output.write('0\t' + str(plaac_prlds[accession]) + '\n')
            else:
                output.write('0\t-\n')
            
    output.close()
    # print(plotting_df)
    return fig3a_plotting_df, fig3b_plotting_df

        
def calculate_deltaPAPA(protein_name, papa_scores, seqs, high_score_positions):

    isoform_accessions = []
    isoform_max_scores = []
    isoform_lengths = []
    isoform_max_score_positions = []
        
    for accession in papa_scores[protein_name]:
        isoform_accessions.append(accession)
        isoform_max_scores.append(papa_scores[protein_name][accession])
        isoform_lengths.append(len(seqs[accession]))
        isoform_max_score_positions.append(high_score_positions[protein_name][accession])

            
    longest_iso_index = isoform_lengths.index( max(isoform_lengths) )
    highest_iso_score_index = isoform_max_scores.index( max(isoform_max_scores) )
    highest_iso_score = isoform_max_scores[ highest_iso_score_index ]
    
    delta_papascores = []
    for score in isoform_max_scores:
        d_papa = highest_iso_score - score
        delta_papascores.append(d_papa)
    
    return isoform_accessions, isoform_lengths, isoform_max_scores, isoform_max_score_positions, delta_papascores, highest_iso_score
    
    
def get_len_longest_isoform(prot_dict):
    max_len = 0
    for prot in prot_dict:
        length = len(prot_dict[prot])
        if length > max_len:
            max_len = length
            
    return max_len
    
def get_protein_seqs(fasta_file, acc2common):
    
    h = open(fasta_file)
    seqs = {}
    for seq_record in SeqIO.parse(h, 'fasta'):
        accession = str(seq_record.id)
        if accession in acc2common:
            common_name = acc2common[accession]
        else:
            common_name = accession
        seq = str(seq_record.seq)
        
        seqs[accession] = seqs.get(accession, seq)
        
    return seqs

    
def get_papa_scores(papa_results_file, acc2common):

    papa_scores = {}
    papa_prlds = {}
    high_score_positions = {}
    h = open(papa_results_file)
    header = h.readline()
    for line in h:  
        items = line.rstrip().split('\t')
        accession = items[0]
        if accession in acc2common:
            common_name = acc2common[accession]
        else:
            common_name = accession

        if items[1] == 'protein length below window size':
            papa_score = -1
        else:
            papa_score = float(items[1])
        
        papa_scores[common_name] = papa_scores.get(common_name, {})
        papa_scores[common_name][accession] = papa_scores[common_name].get(accession, papa_score)
        
        if len(items) < 5:
            position = '-'
        else:
            position = int(items[2])
        high_score_positions[common_name] = high_score_positions.get(common_name, {})
        high_score_positions[common_name][accession] = high_score_positions[common_name].get(accession, position)
        
        if len(items) == 5:
            prlds = items[4].replace('"', '').split('_;_')
            for prld in prlds:
                papa_start, papa_end = prld.split(', ')
                papa_start = int(papa_start[1:])
                papa_end = int(papa_end[:-1])
                
                papa_prlds[accession] = papa_prlds.get(accession, [] )
                papa_prlds[accession].append( (papa_start, papa_end) )
        
    return papa_scores, papa_prlds, high_score_positions
    
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
            plaac_positives.append(acc2common[SEQid])
            prld_boundaries[SEQid] = prld_boundaries.get( SEQid, (int(PRDstart), int(PRDend)) )
        
    h.close()

    return prld_boundaries, plaac_positives
    
    
def check_prld_overlap(papa_prlds, plaac_start, plaac_end):
    toggle = False
    for prld in papa_prlds:
        papa_start, papa_end = prld
        papa_start = int(papa_start)
        papa_end = int(papa_end)
        if (plaac_start >= papa_start and plaac_start <= papa_end) or (plaac_end >= papa_start and plaac_end <= papa_end) or (plaac_start <= papa_start and plaac_end >= papa_end):
            toggle = True
            
    return toggle

def plot_Fig3A(new_small_df, plaac_positives):
         
    #get co-sorted lists
    current_iso_PAPAscore, max_iso_PAPAscore, common_name = zip(*sorted(zip(new_small_df['Min PAPAscore'], new_small_df['Max PAPAscore'], new_small_df['Common_name'])))
    common_name = list(common_name)
    
    #plot barchart
    palette = sns.color_palette()
    for i in range(len(new_small_df['Common_name'])):
        plt.bar(i, current_iso_PAPAscore[i], color=palette[0])
        if current_iso_PAPAscore[i] < 0:
            b = 0
            plt.bar(i, max_iso_PAPAscore[i], bottom=b, color=palette[1])
        else:
            b = current_iso_PAPAscore[i]
            plt.bar(i, (max_iso_PAPAscore[i] - current_iso_PAPAscore[i]), bottom=b, color=palette[1])
        if common_name[i] in plaac_positives:
            plt.text(i, (max_iso_PAPAscore[i]), '*', ha='center', va='bottom', fontsize=12)

        ax = plt.gca()
    plt.xticks(range(len(common_name)), rotation=90, fontsize=8, fontname='Arial')
    names = []
    for x in common_name:
        if x == '7-Mar':
            names.append('MARCH7')
        else:
            names.append(x)
    common_name = [x.upper() if x not in plaac_positives else (x + ' (*)') for x in names]
    minor_labels = []
    major_labels = []
    for i in range(len(common_name)):
        if i%2 == 0:
            major_labels.append(common_name[i])
        else:
            minor_labels.append(common_name[i])
    plt.legend(['Low-scoring Isoform', 'High-scoring Isoform'], loc='right', prop={'family':'Arial', 'size':20})
    ax = plt.gca()
    ax.set_xticklabels(common_name, fontsize=11, family='monospace')
    counter = 1
    for tick in ax.xaxis.get_major_ticks()[1::2]:
        #NOTE: these values were the result of a trial-and-error testing to get the labels on the inside of the axis to offset by similar amounts
        #Basically, it takes the length of the gene name, raises it to the 0.85 power, and multiplies it be -10.5
        #The power portion is a scaling factor based on the length of the gene name
        #The *-10.5 does the main offsetting to get it on the upper side of the axis
        #This scaling is also heavily dependent on the font size, and ABSOLUTELY REQUIRES a monospaced font
        tick.set_pad((len(common_name[counter]) ** 0.85) * -10.5)
        counter += 2

    plt.ylim(-1.2,0.3)
    plt.yticks(fontsize=12, fontname='Arial')
    plt.xlabel('Proteins', fontname='Arial', fontsize=14)
    plt.ylabel('PAPA Scores', fontname='Arial', fontsize=14)
    plt.tight_layout()
    plt.plot([0.05 for i in range(len(common_name))], color='0.2', linestyle='--')
    fig = plt.gcf()
    fig.set_size_inches(20, 10)
    plt.savefig('Fig 3A.tiff', bbox_inches='tight', dpi=300)
    plt.close()
    
def plot_Fig3B(small_df, plaac_positives):

    #get co-sorted lists
    current_iso_PAPAscore, max_iso_PAPAscore, common_name = zip(*sorted(zip(small_df['Min PAPAscore'], small_df['Max PAPAscore'], small_df['Common_name'])))
    
    #plot barchart
    palette = sns.color_palette()
    for i in range(len(small_df['Common_name'])):
        plt.bar(i, current_iso_PAPAscore[i], color=palette[0])
        if current_iso_PAPAscore[i] < 0:
            b = 0
            plt.bar(i, max_iso_PAPAscore[i], bottom=b, color=palette[1])
        else:
            b = current_iso_PAPAscore[i]
            plt.bar(i, (max_iso_PAPAscore[i] - current_iso_PAPAscore[i]), bottom=b, color=palette[1])
        if common_name[i] in plaac_positives:
            plt.text(i, (max_iso_PAPAscore[i]), '*', ha='center', va='bottom')
        ax = plt.gca()
        ax.set_xticklabels(common_name)
        plt.xticks(range(len(common_name)), rotation=90, fontsize=8, fontname='Arial')
    plt.legend(['Low-scoring Isoform', 'High-scoring Isoform'], loc='upper left', prop={'family':'Arial', 'size':10})
    ax = plt.gca()
    ax.set_xticklabels(common_name, fontname='Arial')
    plt.yticks(fontsize=12, fontname='Arial')
    plt.xticks(fontsize=10, fontname='Arial')
    plt.xlabel('Proteins', fontname='Arial', fontsize=14)
    plt.ylabel('PAPA Scores', fontname='Arial', fontsize=14)
    plt.tight_layout()
    plt.plot([0.05 for i in range(len(common_name))], color='0.2', linestyle='--')
    fig = plt.gcf()
    fig.set_size_inches(8, 4)
    plt.savefig('Fig 3B.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def main(papa_results_file, plaac_results_file, fasta_file):
    
    #Load accession to common name conversion dictionary
    acc2common = pickle.load(open('Human_prots_GenBankAccession_to_Common_GeneName.dat', 'rb'))
    
    #Get protein sequences
    seqs = get_protein_seqs(fasta_file, acc2common)
    
    #Get PAPA scores
    papa_scores, papa_prlds, high_score_positions = get_papa_scores(papa_results_file, acc2common)
    
    #Get PLAAC scores
    plaac_prlds, plaac_positives = get_plaac_prld_boundaries(plaac_results_file, acc2common)

    #Compare isoform scores and output Additional file 2
    fig3a_plotting_df, fig3b_plotting_df = compare_isoforms(papa_scores, papa_prlds, plaac_prlds, high_score_positions, seqs)
            
    #Plot high-scoring isoforms
    plot_Fig3A(fig3a_plotting_df, plaac_positives)
    plot_Fig3B(fig3b_plotting_df, plaac_positives)
    

if __name__ == '__main__':
    import sys
    #papa_results_file --> 'Human_Protein_Seqs_ActiveDriverDB_PAPA_results.tsv' (the result of running mPAPA.py on 'Human_Protein_Seqs_ActiveDriverDB.FASTA' with threshold=0.0)
    #plaac_results_file --> 'Human_Protein_Seqs_ActiveDriverDB_PLAAC_results.csv' (the result of running plaac.jar on 'Human_Protein_Seqs_ActiveDriverDB.FASTA' with alpha=0 and the FASTA file used to calculate background amino acid frequencies)
    #fasta_file --> 'Human_Protein_Seqs_ActiveDriverDB.FASTA' (original file from ActiveDriverDB)
    papa_results_file, plaac_results_file, fasta_file = sys.argv[1:]
    main(papa_results_file, plaac_results_file, fasta_file)
