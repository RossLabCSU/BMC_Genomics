
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import seaborn as sns
import numpy as np


def get_len_longest_isoform(prot_dict):
    max_len = 0
    for prot in prot_dict:
        length = len(prot_dict[prot])
        if length > max_len:
            max_len = length
            
    return max_len
  
        
def get_variant_PAPAscores(mutant_papa_results_file):
    
    papa_scores = {}
    fold_index = {}

    h = open(mutant_papa_results_file)
    header = h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        
        if len(items) < 5:
            continue

        accession = items[0]
        scores = items[3][1:-1].split(' ')
        papa_scores[accession] = papa_scores.get(accession, scores)
        
        foldindex_scores = items[4][1:-1].split(' ')
        fold_index[accession] = fold_index.get(accession, foldindex_scores)
            
    h.close()
    
    return papa_scores, fold_index
    
def plotting_v3(protein_name, prot_dict, foldindex_dict, ptm_positions, full_ptm_info):

    colors_index = 0
    legend_labels = []
    legend_colors = []
    isoform_counter = 0
    
    ptm_linestyles = {'Phosphorylation':'solid',
                    'Methylation':'dotted',
                    'Ubiquitination':'dashed',
                    'Acetylation':'dashdot'}
    
    len_longest_isoform = get_len_longest_isoform(prot_dict)
    plt.plot([i for i in range( len_longest_isoform )], [0.05 for i in range( len_longest_isoform )], color='0.50', linestyle='--', linewidth=1.25)
    for i in range(len(ptm_positions)):
        if full_ptm_info[i][-2] == 350 and protein_name == 'HNRNPA1':
            plt.plot([ptm_positions[i]]*10, [j for j in np.linspace(-0.75, 0.75, 10)], color='red', linewidth=1.0, linestyle=ptm_linestyles[ 'Acetylation' ])
        else:
            if full_ptm_info[i][-1] == 'Phosphorylation':
                line = plt.plot([ptm_positions[i]]*10, [j for j in np.linspace(-0.75, 0.75, 10)], color='red', linewidth=0.6, linestyle=ptm_linestyles[ full_ptm_info[i][-1]])
            else:
                line = plt.plot([ptm_positions[i]]*10, [j for j in np.linspace(-0.75, 0.75, 10)], color='red', linewidth=1.0, linestyle=ptm_linestyles[ full_ptm_info[i][-1]])
            
    legend_ghostplots = []
    index = 0
    
    #creates a gradient color palette using seaborns default color palette with variations on the "light" (l) parameter
    #first creates a list that covers the linear space from 0.3-0.9 in increments of (number_of_lines_to_plot / 10)
    #then cycles through the seaborn husl palette and changes the lightness each time the 10 colors are cycled through
    #NOTE: this should properly scale the colors regardless of how many lines you are plotting, although the difference between lines will become less apparent as the # of lines increases
    color_palette = []
    if len(prot_dict)%10 == 0:
        vals = np.linspace(0.3,0.8, len(prot_dict)/10)
    else:
        vals = np.linspace(0.3,0.8, int(len(prot_dict)/10)+1)
    for i in range(len(vals)):
        pal = sns.husl_palette(10, l=vals[i])
        for color in pal:
            color_palette.append(color)
            
    #reverse color_palette list so that dark colors are plotted last
    color_palette = color_palette[::-1]
            
    
    #looping through variants in prot_dict
    for variant in prot_dict:
        #only variants that indicate mutants should have 2 underscores in their name
        if variant.count('_') == 2:
            nm, number, mut = variant.split('_')
            legend_labels.append(mut)
        if variant == 'NM_031157' or variant == 'NM_004960':
            legend_labels.append('Long Isoform')
        if variant == 'NM_001170937':      #FUS only
            legend_labels.append('Medium Isoform')
        if variant == 'NM_002136' or variant == 'NM_001170634':
            legend_labels.append('Short Isoform')
        papa_scores = prot_dict[variant]
        foldindex_scores = foldindex_dict[variant]

        #get coordinates for papascores vs protein position for predicted ordered and disordered regions
        xvals_lol, yvals_lol, xvals_ordered_regions_lol, yvals_ordered_regions_lol = split_list_v2(papa_scores, foldindex_scores)
        if len(xvals_lol) == 0:
            continue

        #main plotting of scores for ordered and disordered regions
        for i in range(len(xvals_ordered_regions_lol)):
            plt.plot(xvals_ordered_regions_lol[i], yvals_ordered_regions_lol[i], color='0.8', linewidth=0.4)
        for i in range(len(xvals_lol)):
            plt.plot(xvals_lol[i], yvals_lol[i], color=color_palette[colors_index], label=legend_labels[index], linewidth=1.25)

        #ghost plots for legend formatting
        legend_ghostplots.append(Line2D([0], [0], color=color_palette[colors_index]))
            
        colors_index += 1
        index += 1
    
    #basic formatting
    ptm_lines = []
    if protein_name == 'HNRNPA1':
        ptms = ['Phosphorylation', 'Methylation', 'Ubiquitination', 'Ubiquitination +\nAcetylation']
        i = 0
        for ls in ['solid', 'dotted', 'dashed', 'dashdot']:
            line = mlines.Line2D([], [], color='red', markersize=15, linestyle=ls, label=ptms[i])
            ptm_lines.append(line)
            i += 1
            
    else:
        ptms = ['Phosphorylation', 'Methylation']
        i = 0
        for ls in ['solid', 'dotted']:
            line = mlines.Line2D([], [], color='red', markersize=15, linestyle=ls, label=ptms[i])
            ptm_lines.append(line)
            i += 1

    second_legend = plt.legend(handles=ptm_lines, bbox_to_anchor=(1.01, 0), loc=3, borderaxespad=0.1, fontsize=9, title='Modifications')

    plt.legend(legend_ghostplots, legend_labels, bbox_to_anchor=(1.01, 1), ncol=2, loc=2, borderaxespad=0.1, fontsize=7, title='Mutations')
    plt.gca().add_artist( second_legend )
    plt.title(protein_name, fontname='Arial', fontsize=16)

    plt.yticks([-0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.10, 0.15], fontsize=12, fontname='Arial')
    plt.title(protein_name, fontname='Arial', fontsize=16)
    plt.xlabel('Amino Acid Position', fontname='Arial', fontsize=14)
    plt.ylabel('PAPA Score', fontname='Arial', fontsize=14)
    
    if protein_name == 'HNRNPA1':
        plt.ylim(-0.35, 0.19)
        plt.xticks([x for x in range(0, 400, 50)], fontsize=12, fontname='Arial')
    else:
        plt.ylim(-0.4, 0.13)
        plt.xticks([x for x in range(0, 600, 50)], fontsize=12, fontname='Arial')

    fig = plt.gcf()
    fig.set_size_inches(6.5, 4)
    
    if protein_name == 'HNRNPA1':
        plt.savefig('Fig 6A.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig 6B.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def split_list_v2(papa_scores, foldindex_scores):
    """Function that splits a list of verbose PAPA scores at 'None' values into separate lists
    for plotting. This allows for plotting with line breaks where the protein was not scored
    by PAPA due to a positive FoldIndex.
    Returns:
        1) a list of lists for x-values, where each value is a position within the protein.
        2) a list of lists for y-values, where each value is the PAPA score at the corresponding position."""
        
    yvals_lol = []
    xvals_lol = []
    yvals_ordered_regions_lol = []
    xvals_ordered_regions_lol = []
    i = 1
    
    while len(papa_scores) > 0:
        yval_nones = []
        xval_nones = []
        while len(papa_scores) > 0 and float(foldindex_scores[0]) > 0:
            yval_nones.append( float(papa_scores.pop(0)) )
            xval_nones.append(i)
            foldindex_scores.pop(0)
            i += 1
        if len(xval_nones) == 1:
            xval_nones.append(i)
            yval_nones.append(float(papa_scores[0]))
        if len(xval_nones) > 0:
            yvals_ordered_regions_lol.append( yval_nones )
            xvals_ordered_regions_lol.append( xval_nones )
            
        yvals = []
        xvals = []
        while len(papa_scores) > 0 and float(foldindex_scores[0]) <= 0:
            yvals.append( float(papa_scores.pop(0)) )
            xvals.append(i)
            foldindex_scores.pop(0)
            i += 1
        if len(xvals) == 1:
            xvals.append(i)
            yvals.append(float(papa_scores[0]))
        if len(xvals) > 0:
            yvals_lol.append(yvals)
            xvals_lol.append(xvals)

    return xvals_lol, yvals_lol, xvals_ordered_regions_lol, yvals_ordered_regions_lol

        
def get_PTM_sites(gene, ptm_results_file):
    
    h = open(ptm_results_file)
    header = h.readline()
    
    if gene == 'HNRNPA1':
        #ONLY GET PTMs FOR THE LONG ISOFORMS
        accession = 'NM_031157'
    else:
        accession = 'NM_004960'
    full_ptm_info = []
    ptm_positions_only = []
    for line in h:
        items = line.rstrip().split('\t')
        if len(items) < 5:
            continue
        acc = items[0]
        if acc == accession:
            common_name = items[1]
            max_papa_score = items[2]
            papa_zerothreshold_prld_seq = items[3]
            prld_boundaries = items[4]
            ptm_sites = items[5]

            mods = ptm_sites.split( '); (' )
            mods = [x.replace('(', '') for x in mods]
            mods = [x.replace(')', '') for x in mods]
            for mod in mods:
                res, location, type = mod.split('_')
                location = int(location)
                tup = (res, location, type)
                full_ptm_info.append(tup)
                ptm_positions_only.append(int(location))
                
    h.close()
    
    return ptm_positions_only, full_ptm_info
            
    
def main(ptm_results_file):

    for gene in ['HNRNPA1', 'FUS']:
    
        #GET PTM INFO FOR HNRNPA1
        ptm_positions_only, full_ptm_info = get_PTM_sites(gene, ptm_results_file)
        
        #GET PAPA SCORES AND CORRESPONDING FOLDINDEX SCORES FOR HNRNPA1
        scores, foldindex_scores = get_variant_PAPAscores(gene + '_All_Variants_PAPA_results_verbose.tsv')
        
        #PLOT HNRNPA1 SCORES AND PTMs
        plotting_v3(gene, scores, foldindex_scores, ptm_positions_only, full_ptm_info)
    

if __name__ == '__main__':
    import sys
    #ptm_results_file --> 'Additional file 4.tsv'
    ptm_results_file = sys.argv[1]
    main(ptm_results_file)
