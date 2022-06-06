import pandas as pd
import numpy as np
import pickle
from Bio import SeqIO


def reverseComplement(string):
    complement_string = str()
    for nu in string:
        nu = nu.upper()
        if nu == 'A':
            complement_string += 'T'
        elif nu == 'T':
            complement_string += 'A'
        elif nu == 'G':
            complement_string += 'C'
        elif nu == 'C':
            complement_string += 'G'
    return complement_string[::-1]


def openHg19Pickle(pickle_path = '../hg19/hg19.pickle'):
    global chromosome_sequences
    with open(pickle_path, 'rb') as f:
        chromosome_sequences = pickle.load(f)


def getChromosomeSequence(chromosome, genome_dir = 'hg19/',memory ='RAM'):
    if memory == 'RAM':
        global chromosome_sequences
        return chromosome_sequences[chromosome]
    elif memory == 'disk':
        choromosome_path = genome_dir + chromosome + '.fa' 
        parsing_file = SeqIO.parse(open(choromosome_path),'fasta')
        for parsing_seq in parsing_file:
            return str(parsing_seq.seq)


class BedFile():
    def __init__(self,file_path):
        self.path = file_path
        self.importBEDFile(self.path)
        
    def importBEDFile(self,file_path):
        file_extension = file_path[-3:].lower()
        if file_extension == 'bed':
            data = pd.read_csv(file_path, sep = '\t', header = None)
        elif file_extension == 'csv':
            data = pd.read_csv(file_path, sep= ',',header = None)
            data = data.drop(0, axis = 0)
        id_column = 3   
        data = data.drop([id_column],axis = 1)
        data.columns = ['chr','start','end','other','strand']
        self.data = data


    def getSequence(self,padding = 0, offset = 0, side = 'all_span', span = 0):
        global chromosome_sequences
        if 'chromosome_sequences' not in globals():
          openHg19Pickle()
        
        sequences = []
        sample_size = self.data.shape[0]
        for i in range(sample_size):
            chromosome = self.data['chr'][i]
            start = self.data['start'][i]
            end = self.data['end'][i]
            strand = self.data['strand'][i]

            if strand == '+':
                if side == 'all_span':
                    begin = start + offset - padding
                    stop = end + offset + padding
                if side == '5_prime':
                    begin = start + offset
                    stop = start + offset + span
                if side == '3_prime':
                    begin = end + offset - span
                    stop = end + offset
                extracted_sequence = getChromosomeSequence(chromosome)[begin:stop].upper()
                sequences.append(extracted_sequence)

            if strand == '-':
                if side == 'all_span':
                    begin = start + offset - padding
                    stop = end + offset + padding
                if side == '3_prime':
                    begin = start + offset
                    stop = start + offset + span
                if side == '5_prime':
                    begin = end + offset - span
                    stop = end + offset
                extracted_sequence = getChromosomeSequence(chromosome)[begin:stop].upper()
                sequences.append(reverseComplement(extracted_sequence))
            
        return sequences
