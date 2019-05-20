
from Bio import SeqIO

def get_nextprot_seqs(nextprot_file):

    h = open(nextprot_file)
    output = open('NextProt_protein_sequences.FASTA', 'w')
    
    #Read through header lines
    for i in range(11):
        h.readline()

    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        temp = id.split(':')
        id = temp[1]
        seq = str(seq_record.seq)
        
        output.write('>' + id + '\n')
        output.write(seq + '\n')
        
    output.close()

if __name__ == '__main__':
    import sys
    #nextprot_file --> nextprot_all.peff
    nextprot_file = sys.argv[1]
    get_nextprot_seqs(nextprot_file)
    
