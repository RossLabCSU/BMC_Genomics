
import random
from itertools import combinations
from Bio import SeqIO

def main(papa_results_file, nextprot_fasta_file):

    papa_scores = get_mPAPA_scores(papa_results_file)
    
    output_file = open('NextProt_PairwiseVariant_ProteinSeqs.FASTA', 'w')
    
    h = open(nextprot_fasta_file)
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        nxp, id = id.split(':')
        seq = str(seq_record.seq)
        if papa_scores[id] == 'protein length below window size' or papa_scores[id] < 0.0:
            continue
            
        desc = str(seq_record.description)
        items = desc.split('\\')
        
        for item in items:
            if 'VariantSimple' in item:
                title, vars = item.split('=')
                vars = vars[1:-1]
                vars = vars.split(')(')
                vars = [x for x in vars if '*' not in x]
                
                
        positions_set = set()
        for var in vars:
            pos, aa = var.split('|')
            positions_set.add(pos)

        #Get all pairwise combinations of variants
        final_variants = get_all_pairwise_combinations(vars)
            
        #Create and output sequence variants
        output_variant_seqs(output_file, final_variants, id, seq)
                  
    
def get_mPAPA_scores(papa_results_file):
    h = open(papa_results_file)
    header = h.readline()
    papa_scores = {}
    for line in h:
        items = line.rstrip().split('\t')
        id = items[0]
        score = float(items[1])
        papa_scores[id] = papa_scores.get(id, score)
        
    return papa_scores


def get_all_pairwise_combinations(vars):
    
    pairwise_var_combs = []
    for var in vars:
        pos, aa = var.split('|')
        for var2 in vars:
            pos2, aa2 = var2.split('|')
            if pos == pos2:
                continue
            pairwise_var_combs.append( (var, var2) )

    return pairwise_var_combs
             

def output_variant_seqs(output_file, final_variants, gene, seq):

    orig_seq = seq[:]
    output_file.write('>' + gene + '\n')
    output_file.write(orig_seq + '\n')
    
    index = 1
    for var_set in final_variants:
        variant_seq = orig_seq[:]
        for var in var_set:
            pos, mut_aa = var.split('|')
            pos = int(pos) - 1
            variant_seq = variant_seq[:pos] + mut_aa + variant_seq[pos+1:]
        
        var_string = ','.join('('+var+')' for var in var_set)
        output_file.write('>' + gene + '_Variant' + str(index) + '_' + var_string + '\n')
        output_file.write(variant_seq + '\n')
        
        index += 1

if __name__ == '__main__':
    import sys
    #papa_results_file --> NextProt_Human_Proteome_PAPA_ZeroThreshold_MergedWindows.tsv, but will be specified by user at mPAPA runtime
    #nextprot_file --> nextprot_all.peff
    papa_results_file, nextprot_fasta_file = sys.argv[1:]
    main(papa_results_file, nextprot_fasta_file)
