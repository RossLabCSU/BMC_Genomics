
from Bio import SeqIO

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

    
def generate_snv_fasta(fasta_file, snv_info_file):
    """Generates a FASTA file with a separate sequence for each disease-associated single nucleotide variant
    NOTE: this is also derived from the FASTA file containing all "high-confidence" protein isoforms."""

    h = open(fasta_file)
    snvs = get_SNV_info(snv_info_file)
    
    output = open('Human_proteome_with_SNVs.fasta', 'w')

    for seq_record in SeqIO.parse(h, 'fasta'):
        prot = str(seq_record.id)
        seq = str(seq_record.seq)
        output.write('>' + prot + '\n')
        output.write(seq + '\n')
        if prot in snvs:
            snv_seqs = []
            for snv in snvs[prot]:
                common_name, aa_before_mut, position, aa_after_mut, phenotype = snv

                #Excludes truncation mutations
                if aa_after_mut == '*':
                    continue
                
                #Excludes mutations that are not single-amino acid substitutions
                try:
                    int(position)
                except:
                    continue
                    
                #necessary because the position of some SNVs was greater than the length of the sequence (e.g. NM_152336 -> length=1067, but annotated with an R1074* mutation)
                if int(position) > len(seq):
                    continue

                mut_ind = int(position) - 1
                if seq[mut_ind] == aa_before_mut:
                    snv_seq = seq[ : mut_ind] + aa_after_mut + seq[mut_ind + 1 : ]
                    snv_seqs.append(snv_seq)
                    output.write('>' + prot + '_' + aa_before_mut + position + aa_after_mut + '\n')
                    output.write(snv_seq + '\n')
            
    output.close()

    
if __name__ == '__main__':
    import sys
    #snv_info_file --> 'variant_summary.txt'
    #fasta_file --> 'Human_Protein_Seqs_ActiveDriverDB.FASTA'
    fasta_file, snv_info_file = sys.argv[1:]
    generate_snv_fasta(fasta_file, snv_info_file)