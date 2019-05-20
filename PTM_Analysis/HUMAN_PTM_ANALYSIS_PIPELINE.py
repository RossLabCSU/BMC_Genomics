
#Global imports
import pickle
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy import stats
import math
from statsmodels.stats import multitest


def get_human_mods(ptm_file):
    """Accepts the ptm_file downloaded from the ActiveDriverDB, parses the PTM info,
    and stores it in a dictionary (mod_dict), where each key is a modification
    type in ['Acetylation', 'Phosphorylation', 'Methylation', 'Ubiquitination'] and 
    each value is a sub-dictionary.
    
    Within the sub-dictionary, each key is a gene/protein accession number, and each value
    is a set of non-redundant modifications for the given protein.
    
    Returns the mod_dict"""
    
    mod_dict = {'Acetylation':{},
                'Phosphorylation':{},
                'Methylation':{},
                'Ubiquitination':{}}
    
    for ptm in mod_dict:
        h = open(ptm_file)
        header = h.readline()

        for line in h:
            accession, position, residue, enzyme, pubmedID, mod_type = line.rstrip().split('\t')
            if ',' in mod_type:
                items = mod_type.split(',')
                for i in items:
                    i = i.title()
                    if i == ptm:
                        mod_dict[ptm][accession] = mod_dict[ptm].get(accession, set())
                        mod_dict[ptm][accession].add( (residue, position) )

            else:
                mod_type = mod_type.title()
                if mod_type == ptm:
                    mod_dict[ptm][accession] = mod_dict[ptm].get(accession, set())
                    mod_dict[ptm][accession].add( (residue, position) )

        h.close()

    return mod_dict


def mod_enrichment_stats(aa_counts_proteome, perc_modified_proteome, aa_counts_prlds, perc_modified_prlds):

    output = open('Additional file 5.csv' ,'w')
    mod_info = {'Modification type':[],
                'Amino acid':[],
                'Num mod proteome':[],
                'Percent mod proteome':[],
                'Residue total proteome':[],
                'Num mod prlds':[],
                'Percent mod prlds':[],
                'Residue total prlds':[],
                'Odds ratio':[],
                'ln(OR)':[],
                'P-value':[],
                'BH adjusted pval':[]
                }
                
    output.write('Modification Type,Amino Acid,# of Residues Modified (Whole Proteome),Percent Modified (Whole Proteome),Residue Total in Proteome,# of Residues Modified (PrLDs),Percent Modified (PrLDs),Residue Total in PrLDs,Odds Ratio,ln(OR),P-value,Benjamini-Hochberg Adjusted p-value\n')
    for mod_type in perc_modified_proteome:
        index = -1
        for res_tuple in perc_modified_proteome[mod_type]:
            index += 1

            aa = res_tuple[0]
            num_mod_PrLD = round( (perc_modified_prlds[mod_type][index][1]/100) * aa_counts_prlds[aa] )
            num_mod_nonPrLD = round( (res_tuple[1]/100) * aa_counts_proteome[aa] ) - num_mod_PrLD
            num_NOTmod_PrLD = aa_counts_prlds[aa] - num_mod_PrLD
            num_NOTmod_nonPrLD = (aa_counts_proteome[aa] - aa_counts_prlds[aa]) - num_mod_nonPrLD
            
            if round( (res_tuple[1]/100) * aa_counts_proteome[aa] ) < 100:
                continue
                
            #form contingency table for Fishers exact test
            table = []
            table.append([num_mod_PrLD, num_mod_nonPrLD])
            table.append([num_NOTmod_PrLD, num_NOTmod_nonPrLD])
            oddsratio, pvalue = stats.fisher_exact(table)
            
            if oddsratio != 0:
                lnOR = math.log(oddsratio)
            else:
                lnOR = 'NaN'
                
            mod_info['Modification type'].append(mod_type)
            mod_info['Amino acid'].append(aa)
            mod_info['Num mod proteome'].append(str(num_mod_nonPrLD))
            mod_info['Percent mod proteome'].append(str(res_tuple[1]))
            mod_info['Residue total proteome'].append(str(aa_counts_proteome[aa]))
            mod_info['Num mod prlds'].append(str(num_mod_PrLD))
            mod_info['Percent mod prlds'].append(str(perc_modified_prlds[mod_type][index][1]))
            mod_info['Residue total prlds'].append(str(aa_counts_prlds[aa]))
            mod_info['Odds ratio'].append(str(oddsratio))
            mod_info['ln(OR)'].append(str(lnOR))
            mod_info['P-value'].append(pvalue)
            

    reject_H0, adjusted_pvals, sidak_new_alpha, bonferroni_new_alpha = multitest.multipletests(mod_info['P-value'], method='fdr_bh')
    for pval in adjusted_pvals:
        mod_info['BH adjusted pval'].append(str(pval))
        
    mod_info['P-value'] = [str(x) for x in mod_info['P-value']]
    
    for i in range(len(mod_info['Modification type'])):
        output.write(mod_info['Modification type'][i])
        for key in list(mod_info.keys())[1:]:
            output.write(',' + mod_info[key][i])
        output.write('\n')
    output.close()
    
    return mod_info
    

def check_mods(gene, window_indices, mod_dict):
    """function that accepts indicies for the highest-scoring PAPA window and 
        returns a list of all modifications within that window"""

    #each hit will be a tuple, where the first value is the entry in the mod dictionary and the second value is the name of the dataset it came from
    master_hits = []
    all_mods = []
    for window in window_indices:
        hits = []
        for mod_name in mod_dict:
            d = mod_dict[mod_name]
            if gene in d:
                for mod in d[gene]:
                    res = mod[0]
                    position = mod[1]
                    mod_index = int(mod[-1]) - 1
                    if mod_index >= window[0] and mod_index <= window[1]:
                        master_hits.append((res, position, mod_name))
                    all_mods.append((res, position, mod_name))

    return master_hits, all_mods
           
           
def add_PTM_and_PLAAC_annotations(papa_results_file, plaac_results_file, mod_dict):

    plaac_prlds = get_plaac_prld_boundaries(plaac_results_file)
    papa_prlds = {}
    
    accession2common_name = pickle.load(open('Human_prots_GenBankAccession_to_Common_GeneName.dat', 'rb'))
    
    h = open(papa_results_file)
    output = open('Additional file 4.tsv', 'w')
    output.write('Accession Number\tGene\tMax PAPA Score\tHighest-scoring Window Sequence(s)\tHighest-scoring Window(s) Indices\tModifications within Highest-scoring or all PAPA-positive Window(s)\n')

    #tracks values for each gene
    master_dict = {}
    
    #run loop for each line in PAPA output file
    #each line corresponds to a separate gene
    h.readline()
    for line in h:
        items = line.rstrip().split('\t')
        
        #proteins that are shorter than 41aa are designated 'protein length below window size'
        #proteins without a disordered region as classified by FoldIndex are assigned a place-holder score of '-1.0'
        #output the gene and 'protein length below window size' or '-1.0' before continuing to the next iteration
        if items[1] == 'protein length below window size' or items[1] == '-1.0':
            if gene in accession2common_name:
                output.write(items[0] + '\t' + accession2common_name[gene] + '\t' + items[1] + '\t-\t-\t-')
                if gene in plaac_prlds:
                    output.write('\t1\t' + str(plaac_prlds[gene]) + '\n')
                else:
                    output.write('\t0\t' + '-' + '\n')
            else:
                output.write(items[0] + '\t' + items[0] + '\t' + items[1] + '\t-\t-\t-')
                if gene in plaac_prlds:
                    output.write('\t1\t' + str(plaac_prlds[gene]) + '\n')
                else:
                    output.write('\t0\t' + '-' + '\n')
            continue        

        gene, high_score, position, prld_seq, windows = items
        if float(high_score) < 0.0:
            prld_seq = '-'
        parsed_windows = windows.replace(', ', '-')
        items = [gene, high_score, prld_seq, windows]
        
        parsed_windows = parsed_windows.split('_;_')
        window_indices = []
        for window in parsed_windows:
            lower_bound, upper_bound = window[1:-1].split('-')
            window_indices.append( (int(lower_bound), int(upper_bound)) )
        
        #get posttranslational modification sites within the PAPA-positive windows
        mods_prld, all_mods = check_mods(gene, window_indices, mod_dict)
        if float(high_score) > 0.0:
            papa_prlds[gene] = papa_prlds.get(gene, [high_score, prld_seq, windows, str(mods_prld)[1:-1]])
        
        #START OUTPUT===========================
        if gene in accession2common_name:
            output.write(gene + '\t' + accession2common_name[gene])
        else:
            output.write(gene + '\t' + gene)
        
        for item in items[1:]:
            output.write('\t' + item)
            
        if len(mods_prld) >= 2:
            output.write('\t' + str(mods_prld[0]).replace(', ', '_').replace("'", '') )
            for mod in mods_prld[1:]:
                output.write('; ' + str(mod).replace(', ', '_').replace("'", '') )
        elif len(mods_prld) == 1:
            output.write('\t' + str(mod).replace(', ', '_').replace("'", '') )
        else:
            output.write('\t-')
            
        if gene in plaac_prlds:
            output.write('\t1\t' + str(plaac_prlds[gene]) + '\n')
        else:
            output.write('\t0\t' + '-' + '\n')
 
    output.close()
    h.close()

    return papa_prlds


def get_plaac_prld_boundaries(plaac_results_file):
    """gets PLAAC PrLD boundaries from the PLAAC output file
    Only stores proteins that have a PrLD (i.e. those with a COREscore != NaN in the PLAAC output)
    Each key = protein name
    Each value = (PrLD_start, PrLD_end)
    """
    
    h = open('Human_Protein_Seqs_Corresponding_to_PTM_dataset_PLAAC_Output.csv')
    
    prld_boundaries = {}
    
    for line in h:

        if line.startswith('#') or line.startswith('SEQid'):
            continue
            
        items = line.rstrip().split(',')
        if len(items) != 38:
            continue
            
        SEQid, MW, MWstart, MWend, MWlen, LLR, LLRstart, LLRend, LLRlen, NLLR, VITmaxrun, COREscore, COREstart, COREend, CORElen, PRDscore, PRDstart, PRDend, PRDlen, PROTlen, HMMall, HMMvit, COREaa, STARTaa, ENDaa, PRDaa, FInumaa, FImeanhydro, FImeancharge, FImeancombo, FImaxrun, PAPAcombo, PAPAprop, PAPAfi, PAPAllr, PAPAllr2, PAPAcen, PAPAaa = items
       
        if int(PRDstart) != 0:
            prld_boundaries[SEQid] = prld_boundaries.get( SEQid, (int(PRDstart), int(PRDend)) )
        
    h.close()

    return prld_boundaries
    

def barchart_lnOR_PTMs(mod_info):
              
    lnORs = []
    standard_errors = []
    confidence_ints = []
    mods = []

    for i in range(len(mod_info['Modification type'])):
        mod_type = mod_info['Modification type'][i]
        mod_res = mod_info['Amino acid'][i]
        #integers in the file are being read as floats with many trailing zero's, so I had to first convert them to floats, then convert them to integers (which rounds off the 0's)
        num_mod_proteome = int(mod_info['Num mod proteome'][i])
        modifiable_proteome = int(mod_info['Residue total proteome'][i])
        num_mod_prld = int(mod_info['Num mod prlds'][i])
        modifiable_prld = int(mod_info['Residue total prlds'][i])
        lnOR = float(mod_info['ln(OR)'][i])
        adj_pval = float(mod_info['BH adjusted pval'][i])
        standard_error = (( (1/num_mod_prld) + (1/(modifiable_prld-num_mod_prld)) + (1/num_mod_proteome) + (1/(modifiable_proteome-num_mod_proteome)) )**0.5)
        
        lnORs.append(lnOR)
        standard_errors.append(standard_error)
        mods.append(mod_res + ', ' + mod_type)
        confidence_ints.append(standard_error*1.96)
        
    master_list = []
    for i in range(len(lnORs)):
        master_list.append( (lnORs[i], standard_errors[i], confidence_ints[i], mods[i]) )
        
    master_list.sort()
    
    #PLOTTING
    labels = []
    for i in range(len(lnORs)):
        lnOR = master_list[i][0]
        standard_error = master_list[i][1]
        confidence_int = master_list[i][2]
        mod = master_list[i][3]
        labels.append(mod)
        if lnOR > 0:
            plt.bar(mod, lnOR, color='xkcd:sky blue', yerr=standard_error, ecolor='black')
        else:
            plt.bar(mod, lnOR, color='xkcd:scarlet', yerr=standard_error, ecolor='black')

    ax = plt.gca()
    ax.set_xticklabels(labels)
    plt.xlabel('Modification Type', fontname='Arial', fontsize=14)
    plt.xticks([x for x in range(len(lnORs))], fontname='Arial', rotation=90)
    plt.yticks(fontname='Arial')
    plt.ylabel('Degree of Enrichment in PrLDs', fontname='Arial', fontsize=14)
    plt.savefig('Fig 5B.tiff', bbox_inches = 'tight', dpi=600)
    

def get_modifiable(mod_types, mod_dict):
    """Returns a dictionary where each key is a modification type and the corresponding value is
    a set of modifiable residues for each modification type.
    e.g. dictionary['Phosphorylation'] might contain {'S', 'T', 'Y'}"""
    
    modifiable = {}
    for mod_type in mod_types:
        modified_residues = set()
        for gene in mod_dict[mod_type]:
            for mod in mod_dict[mod_type][gene]:
                modified_residues.add(mod[0])

        modifiable[mod_type] = modifiable.get(mod_type, modified_residues)
        
    return modifiable
    
        
def get_percentage_modified_proteome(proteome_file, modifiable, mod_types, mod_dict):
    """For each modification type, this gets the percentage of possible modifiable residues
    that are actually modified in the proteome.
    
    e.g. this might get the percentage of serine residues that are phosphorylated.
    
    Returns a dictionary with modification types as keys and a list of tuples as values,
    where each tuple contains the amino acid as the first element and the
    percentage modified (from 0 to 100) as the second element.
    Also returns a dictionary with each amino acid as the keys and the raw counts
    within the proteome as the values."""
        
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_counts_proteome = {x:0 for x in amino_acids}
    h = open(proteome_file)
    for seq_record in SeqIO.parse(h, 'fasta'):
        seq = str(seq_record.seq)
        for aa in amino_acids:
            aa_counts_proteome[aa] += seq.count(aa)
    h.close()
    
    perc_modified = {}
    for mod_type in mod_types:
        perc_modified[mod_type] = perc_modified.get(mod_type, [])
        for aa in modifiable[mod_type]:
            num_modified = 0
            for gene in mod_dict[mod_type]:
                for mod in mod_dict[mod_type][gene]:
                    if mod[0] == aa:
                        num_modified += 1

            perc = (num_modified / aa_counts_proteome[aa]) * 100
            perc_modified[mod_type].append( (aa, perc) )
    
    return aa_counts_proteome, perc_modified
    
def get_percentage_modified_prlds(modifiable, prlds, mod_types, mod_dict):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aa_counts_prlds = {x:0 for x in amino_acids}
    for gene in prlds:
        prld_seq = prlds[gene][1]
        for aa in amino_acids:
            aa_counts_prlds[aa] += prld_seq.count(aa)
    
    perc_modified = {}
    total_modified = {}
    for mod_type in mod_types:
        perc_modified[mod_type] = perc_modified.get(mod_type, [])
        total_modified[mod_type] = total_modified.get(mod_type, [])

        for aa in modifiable[mod_type]:
            num_modified = 0
            for gene in prlds:
                # skip PrLDs without mods
                if len(prlds[gene]) < 4:
                    continue
                    
                mods = prlds[gene][-1].split( '), (' )
                mods = [x.replace('(', '') for x in mods]
                mods = [x.replace("'", '') for x in mods]
                mods = [x.replace(')', '') for x in mods]

                for mod in mods:
                    if mod == '':
                        continue
                    res, location, type = mod.split(', ')

                    if res == aa and type == mod_type:
                        num_modified += 1

            perc = (num_modified / aa_counts_prlds[aa]) * 100
            perc_modified[mod_type].append( (aa, perc) )
            total_modified[mod_type].append( (aa, num_modified) )
    
    return aa_counts_prlds, perc_modified, total_modified
 

def crosscheck_mod_with_PrLDs(gene_mod, gene, prlds):

    #prld info with len < 4 are proteins for which there are no PTMs
    if len(prlds[gene]) < 4:
        return False
        
    else:
        prld_bounds = prlds[gene][2].split(';')

        for bound in prld_bounds:
            bound = bound.replace('(', '')
            bound = bound.replace(')', '')
            bound = bound.replace('_', '')
            lower_bound, upper_bound = bound.split(', ')

            if int(gene_mod[-1]) >= int(lower_bound) and int(gene_mod[-1]) <= int(upper_bound):
                return True
                
        return False

    
def plot_histos(prld_modcounts):

    for mod_type in prld_modcounts:
        nonzero_counts = [x for x in prld_modcounts[mod_type] if x != 0 ]
        if len(nonzero_counts) ==0:
            continue

        plt.hist(nonzero_counts, bins=np.arange(1,np.max(nonzero_counts)+2)-0.5, rwidth=0.8)

        plt.title(mod_type, fontsize=18, fontname='Arial')
        
        plt.xlim(0, max(nonzero_counts)+1)
        plt.yticks(fontname='Arial', fontsize=14)
        
        #if the x-axis is small (<=20), force discrete integer ticks by setting the range and step size (otherwise, some ticks would be in decimal increments).
        #if the x-axis is large, allow auto-scaling of the axis
        if max(nonzero_counts) <= 20:
            plt.xticks(range(0, max(nonzero_counts)+1, 1),fontname='Arial', fontsize=14)
        else:
            plt.xticks(fontname='Arial')
        plt.xlabel('# of Modifications within PrLD', fontname='Arial', fontsize=16)
        plt.ylabel('Frequency', fontname='Arial', fontsize=16)

        
        plt.savefig(mod_type + '_Fig5A.tiff', bbox_inches ='tight', pad_inches=0, dpi=600)
        plt.close()


def main_PTM_desc_stats(mod_dict, proteome_file, prld_info):

    mod_types = list(mod_dict.keys())

    modifiable = get_modifiable(mod_types, mod_dict)
    aa_counts_proteome, perc_modified_proteome = get_percentage_modified_proteome(proteome_file, modifiable, mod_types, mod_dict)
    aa_counts_prlds, perc_modified_prlds, total_modified_prlds = get_percentage_modified_prlds(modifiable, prld_info, mod_types, mod_dict)
    
    total_mods = 0
    mods_in_prlds = 0
    mods_outside_prlds = 0
    prld_modcounts = {}
    for mod in mod_types:
        prld_modcounts[mod] = prld_modcounts.get(mod, [])
        prld_total = 0
        prld_density = 0
        non_prld_total = 0
        non_prld_density = 0
        
        h = open(proteome_file)
        for seq_record in SeqIO.parse(h, 'fasta'):
            gene = str(seq_record.id)
            if gene in mod_dict[mod]:
                mods = mod_dict[mod][gene]
                for gene_mod in mods:
                    total_mods += 1
                if gene in prld_info:
                    mod_count = 0
                    for gene_mod in mods:
                        mod_in_prld = crosscheck_mod_with_PrLDs(gene_mod, gene, prld_info)
                        if mod_in_prld is True:
                            mod_count += 1
                            mods_in_prlds += 1
                        else:
                            mods_outside_prlds +=1
                    prld_modcounts[mod].append(mod_count)
        h.close()

    plot_histos(prld_modcounts)
    
    return aa_counts_proteome, perc_modified_proteome, aa_counts_prlds, perc_modified_prlds


def main(proteome_file, ptm_file, papa_results_file, plaac_results_file):
    
    #need to specify file name downloaded from ActiveDriver with PTM info (link titled "Protein PTM sites*" on https://www.activedriverdb.org/download/
    #parse PTM sites from ActiveDriverDB downloaded file
    mod_dict = get_human_mods(ptm_file)
    
    #Add PTM and PLAAC info to the original mPAPA results file, and gather prld info
    prld_info = add_PTM_and_PLAAC_annotations(papa_results_file,plaac_results_file, mod_dict)
    
    #Run descriptive stats and plot PTM distributions within mPAPA=0.0 threshold PrLDs
    aa_counts_proteome, perc_modified_proteome, aa_counts_prlds, perc_modified_prlds = main_PTM_desc_stats(mod_dict, proteome_file, prld_info)
    
    #Run PrLD enrichment/depletion stats within mPAPA=0.0 threshold PrLDs, with Benjamini-Hochberg correction for multiple hypothesis testing
    mod_info = mod_enrichment_stats(aa_counts_proteome, perc_modified_proteome, aa_counts_prlds, perc_modified_prlds)
    
    #Plot PTM enrichment/depletion barchart
    barchart_lnOR_PTMs(mod_info)
    
if __name__ == '__main__':
    import sys
    #proteome_file --> "Human_Protein_Seqs_ActiveDriverDB.FASTA"
    #ptm_file --> "site_table.tsv" (downloaded from ActiveDriverDB)
    #papa_results_file --> 'Human_Protein_Seqs_ActiveDriverDB_PAPA_results.tsv' (the result of running mPAPA.py on 'Human_Protein_Seqs_ActiveDriverDB.FASTA' with threshold=0.0)
    #plaac_results_file --> "Human_Protein_Seqs_ActiveDriverDB_PLAAC_Output.csv" (the result of running plaac.jar on 'Human_Protein_Seqs_ActiveDriverDB.FASTA' with alpha=0 and the FASTA file used to calculate background amino acid frequencies)
    proteome_file, ptm_file, papa_results_file, plaac_results_file = sys.argv[1:]
    main(proteome_file, ptm_file, papa_results_file, plaac_results_file)