
from Bio import SeqIO
import pickle
        
def main(proteome_file):

    acc2common = pickle.load(open('Human_prots_GenBankAccession_to_Common_GeneName.dat', 'rb'))
    accessions = []
    for acc in acc2common:
        if acc2common[acc] == 'HNRNPA1':
            accessions.append(acc)
    
    h = open(proteome_file)
    output = open('HNRNPA1_All_variants.FASTA', 'w')
    for seq_record in SeqIO.parse(h, 'fasta'):
        seqid = str(seq_record.id)
        seq = str(seq_record.seq)
        for acc in accessions:
            if acc in seqid:
                output.write('>' + seqid + '\n')
                output.write(seq + '\n')
                
    output.close()


if __name__ == '__main__':
    import sys
    #proteome_file --> "Human_proteome_with_SNVs.FASTA"
    proteome_file = sys.argv[1]
    main(proteome_file)
