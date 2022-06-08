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
        elif nu == 'N':
            complement_string += 'N'
    return complement_string[::-1]


def openHg19Pickle(pickle_path = 'hg19/hg19.pickle'):
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
            data = data.drop(3, axis = 1)
        elif file_extension == 'csv':
            data = pd.read_csv(file_path, sep= ',',header = None)
        data.rename(columns = {0:'chr', 1:'start' ,2:'end'},inplace = True)
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


class ExperimentDataset(BedFile):
    def __init__(self,dataset_path):
        super().__init__(dataset_path)
        self.renameData()
    
    def renameData(self,new_column = ['chr','start','end','bp_position','strand','score']):
        self.data.columns = new_column

    def makeAnnotatedData(self,window_size):
        data = pd.DataFrame()
        data['chromosome'] = self.data['chr']
        data['sequence'] = self.getSequence(side = '3_prime',span = window_size)
        data['annotation'] = self.getAnnotation(window_size)
        data['score'] = self.data['score']
        return data

    def getAnnotation(self,window_size):
        annotations = []
        for i in range(self.data.shape[0]):
            annotation = [0]*window_size
            branch_point_distance_to_3_prime = self.getDistanceTo3Prime(index = i)
                      
            if branch_point_distance_to_3_prime <= window_size:
                annotation[-branch_point_distance_to_3_prime] = 1
            else:
                annotation = [branch_point_distance_to_3_prime]
            annotations.append(np.array(annotation))
        return annotations
       
    def getDistanceTo3Prime(self, index):
        strand = self.data['strand'][index]
        start = self.data['start'][index]
        end = self.data['end'][index]
        bp = self.data['bp_position'][index]
        
        if strand == '+':
            distance_to_3_prime = end - bp 
        elif strand == '-':
            distance_to_3_prime = bp - start + 1        
        return distance_to_3_prime

    def makeDeduplicatedData(self,window_size):
        if 'annotated_data' not in self.__dir__():
            self.annotated_data = self.makeAnnotatedData(window_size)
        data = self.annotated_data

        deduplicated_sequences = ['pseudo-sequence'] # initialize with pseudo-sequence
        merged_annotations = []
        chromosome = []
        merged_scores = []
        for i in range(data.shape[0]):
            main_sequence = deduplicated_sequences[-1]
            current_sequence = data['sequence'][i]
            current_annotation = data['annotation'][i]
            current_chromosome = data['chromosome'][i]
            current_score = current_annotation*data['score'][i]

            sequence_is_annotated = (len(current_annotation) != 1)
            if sequence_is_annotated:
                if current_sequence != main_sequence:
                    deduplicated_sequences.append(current_sequence)
                    merged_annotations.append(current_annotation)
                    chromosome.append(current_chromosome)
                    merged_scores.append(current_score)
                elif current_sequence == main_sequence:
                    merged_annotations[-1] = np.logical_or(merged_annotations[-1],current_annotation).astype('f')
                    merged_scores[-1] = np.maximum(merged_scores[-1],current_score)
            else:
                continue
        deduplicated_sequences = deduplicated_sequences[1:] # remove pseudo-sequence

        df = pd.DataFrame()
        df['chromosome'] = chromosome
        df['sequence'] = deduplicated_sequences
        df['annotation'] = merged_annotations
        df['score'] = merged_scores
        return df
