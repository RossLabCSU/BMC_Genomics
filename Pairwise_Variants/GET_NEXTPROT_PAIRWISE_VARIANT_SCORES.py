
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import pandas as pd
import numpy as np
import statistics

    
def plot_num_unique_scores_histo(df):

    unique_scores_per_gene = []
    for gene in df:
        scores_set = set(df[gene])
        unique_scores_per_gene.append(len(scores_set))
        
    scores_counts = [unique_scores_per_gene.count(score) for score in unique_scores_per_gene ]
    arr=plt.hist(unique_scores_per_gene, bins=50, align='mid')
    plt.ylabel('Frequency', fontname='Arial', fontsize=14)
    plt.xlabel('# of Unique PAPA Scores', fontname='Arial', fontsize=14)
    ax = plt.gca()
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.tight_layout()
    plt.savefig('Fig 2A.tiff', bbox_inches='tight', dpi=600)
    plt.close()
  

def plot_minmax_orig_distributions(minmax_scores, orig_scores):

    mins = []
    maxes = []
    origs = []
    for prot in orig_scores:
        if prot in minmax_scores:
            min, max = minmax_scores[prot]
            #FILTER OUT PROTEINS FOR WHICH THE MINIMUM VALUE IS -1.0 (CORRESPONDING TO PROTS WITH VARIANTS THAT ELIMINATE DISORDERED REGION)
            if min == -1.0:
                continue
                
            origs.append(orig_scores[prot])
            mins.append(min)
            maxes.append(max)

    sns.distplot(origs, label='Original Scores', kde=False)
    sns.distplot(mins, label='Estimated Minimum Scores', bins=75, kde=False)
    sns.distplot(maxes, label='Estimated Maximum Scores', kde=False)
    
    plt.legend(prop={'family':'Arial', 'size':12})
    plt.xlabel('PAPA Score', fontname='Arial', fontsize=14)
    plt.ylabel('Frequency', fontname='Arial', fontsize=14)
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.tight_layout()
    plt.savefig('Fig 2C.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
def plot_range_distribution(minmax_scores):

    diffs = []
    for prot in minmax_scores:
        min, max = minmax_scores[prot]
        diffs.append( float(max) - float(min) )
            
    #FILTER OUT RANGES GREATER THAN 1.0 (CORRESPONDING TO PROTS WITH VARIANTS THAT ELIMINATE DISORDERED REGION)
    diffs = [x for x in diffs if x<1.0]

    plt.hist(diffs, bins=40)
    
    plt.xlabel('PAPA Score Range', fontname='Arial', fontsize=14)
    plt.ylabel('Frequency', fontname='Arial', fontsize=14)
    plt.xticks(fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.tight_layout()
    plt.savefig('Fig 2B.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    #Calculate mean, median, and standard deviation for differences in PAPA score
    ave = statistics.mean(diffs)
    med = statistics.median(diffs)
    std = statistics.stdev(diffs)
    
    print('\nAverage: ' + str(ave))
    print('Median: ' + str(med))
    print('Standard Deviation: ' + str(std) + '\n')
    

def plot_minmax_prot_set(var_scores, minmax_scores, orig_scores):

    proteins = ['EWSR1', 'TAF15', 'TARDBP', 'FUS', 'HNRNPA1', 'HNRNPA2B1', 'HNRNPDL', 'TIA1']
    genename_to_nextprot = pickle.load(open('CommonGeneName_to_NextProtID.dat', 'rb'))
    
    index = 0
    x_labels = []
    for prot in proteins:
        nextprot_ids = []
        genenames = []
        minmaxes = []
        for g in genename_to_nextprot:
            if prot in g and genename_to_nextprot[g] in var_scores:
                nextprot_ids.append( genename_to_nextprot[g] )
                genenames.append( g )
                minmaxes.append( (min(var_scores[genename_to_nextprot[g]]), max(var_scores[genename_to_nextprot[g]])) )

        x_labels, index = plot_individ_prot(nextprot_ids, genenames, index, x_labels, minmax_scores, orig_scores)
        
    
    plt.xticks([x for x in range(len(x_labels))], x_labels, rotation=90, fontname='Arial', fontsize=12)
    plt.yticks(fontname='Arial', fontsize=12)
    plt.xlabel('Protein', fontname='Arial', fontsize=14)
    plt.ylabel('PAPA Score', fontname='Arial', fontsize=14)
    plt.ylim(-0.025, 0.199)
    plt.legend(['Estimated Minimum Score', 'Estimated Maximum Score'], prop={'family':'Arial', 'size':12})
    fig = plt.gcf()
    fig.set_size_inches(10, 7)
    plt.tight_layout()
    plt.savefig('Fig 2D.tiff', bbox_inches='tight', dpi=600)
    plt.close()
                
        
def plot_individ_prot(nextprot_ids, genenames, index, x_labels, minmax_scores, orig_scores):

    colors = sns.color_palette(palette='pastel')

    for i in range(len(nextprot_ids)):
        genename = genenames[i]
        nextprot_id = nextprot_ids[i]
        if nextprot_id not in minmax_scores or 'HNRNPA1L' in genename:
            continue
        
        orig_score = orig_scores[nextprot_id]
        min, max = minmax_scores[nextprot_id]

        plt.bar(index, orig_score - min, bottom=min, color=colors[1])
        plt.bar(index, max - orig_score, bottom=orig_score, color=colors[2])
        
        x_labels.append(genename.replace('_', ' '))
        index += 1
    
    return x_labels, index
    
def get_orig_scores(orig_papa_results_file):

    df = {}
    h = open(orig_papa_results_file)
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        id = items[0]
        score = items[1]
        df[id] = df.get(id, float(score))
        
    h.close()

    return df


def get_scores(var_papa_results_file):

    orig_scores = {}
    var_scores = {}
    # for split in range(1,3):
        # h = open('NextProt_protein_sequences_PairwiseComboVariantSpace_ZeroThreshold_PrLDs_Split' + str(split) + '_PAPA_Scores.tsv')
        #SWITCH TO THIS LINE WHEN FINALIZED:
    h = open(var_papa_results_file)

    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        id = items[0]
        score = items[1]

        if 'Variant' in id:
            #NOTE: the '_' delimiter also splits the accession number used as gene identifiers (e.g. 'NX_Q01844')
            acc_part1, acc_part2, var_num, muts = id.split('_')
            gene = acc_part1 + '_' + acc_part2

            var_scores[gene] = var_scores.get(gene, [])
            var_scores[gene].append(float(score))
            
        else:
            orig_scores[id] = orig_scores.get(id, float(score) )
        
    h.close()

    minmax_scores = {}
    for prot in var_scores:
        minmax_scores[prot] = minmax_scores.get(prot, ( min(var_scores[prot]), max(var_scores[prot]) ) )
    
    return var_scores, minmax_scores, orig_scores
        
def output_nextprot_summary(minmax_scores, orig_scores, var_scores):

    nextprotID_to_common_name = pickle.load(open('NextProtID_to_CommonGeneName.dat', 'rb'))

    prld_candidates = []
    for prot in orig_scores:
        if float(orig_scores[prot]) > 0.0:
            prld_candidates.append(prot)

    num_variants_dict = {}
    for prot in prld_candidates:
        if prot in var_scores:
            num_variants_dict[prot] = num_variants_dict.get(prot, len(var_scores[prot]))
        else:
            num_variants_dict[prot] = num_variants_dict.get(prot, 0)
    
    #OUTPUT
    output = open('Additional file 1.tsv', 'w')
    output.write('ProteinID\tOriginal PAPA Score\tSeq Variants Minimum PAPA Score\tSeq Variants Maximum PAPA Score\tNumber of Seq Variants Analyzed\n')
    for prot in prld_candidates:
        output.write(prot + '\t')
        if prot in nextprotID_to_common_name:
            output.write(nextprotID_to_common_name[prot] + '\t')
        else:
            output.write(prot + '\t')
            
        if num_variants_dict[prot] > 0:
            output.write(str(orig_scores[prot]) +  '\t' + str(minmax_scores[prot][0]) + '\t' + str(minmax_scores[prot][1]) + '\t' + str(num_variants_dict[prot]) + '\n')
        else:
            output.write(str(orig_scores[prot]) + '\tNaN\tNaN\t0\n')
    output.close()
    

def main(orig_papa_results_file, var_papa_results_file):

    # orig_scores = get_orig_scores(orig_papa_results_file)
    var_scores, minmax_scores, orig_scores = get_scores(var_papa_results_file)
    
    #Plot Fig 2A
    plot_num_unique_scores_histo(var_scores)
    
    #Plot Fig 2B
    plot_range_distribution(minmax_scores)
    
    #Plot Fig 2C
    plot_minmax_orig_distributions(minmax_scores, orig_scores)
    
    #Plot Fig 2D
    plot_minmax_prot_set(var_scores, minmax_scores, orig_scores)
    
    #Output Additional file 1
    output_nextprot_summary(minmax_scores, orig_scores, var_scores)
    
if __name__ == '__main__':
    import sys
    #orig_papa_results_file --> "NextProt_proteome_PAPA_results.tsv"
    #var_papa_results_file --> "NextProt_PairwiseVariant_ProteinSeqs_PAPA_results.tsv"
    orig_papa_results_file, var_papa_results_file = sys.argv[1:]
    main(orig_papa_results_file, var_papa_results_file)
