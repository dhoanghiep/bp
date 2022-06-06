from Bio import SeqIO
import pickle
import os

def extractChromosomeFromGenome(genome_path):
    hg19 = SeqIO.parse(open(genome_path),'fasta')
    for chromosome in hg19:
        name = chromosome.id
        sequence = str(chromosome.seq)
        with open(name+'.fa', 'w') as f:
            f.write('>'+name+'\n'+sequence)


def hg19ToPythonDict(genome_directory = '../hg19/'):
    chromosome_sequences = {}
    for chr_fa in os.listdir(genome_directory):
        if chr_fa[-2:] == 'fa':
            chromosome = chr_fa[:-3]
            print(chromosome)
            with open(genome_directory+chr_fa, 'r') as f:
                chromosome_sequences[chromosome] = f.read().split('\n')[1]
    with open(genome_directory+'hg19.pickle', 'wb') as f:
        pickle.dump(chromosome_sequences, f, protocol=pickle.HIGHEST_PROTOCOL)
