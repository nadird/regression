import sys
import os


import os
import numpy as np
from Bio import SeqIO
import glob
from pymsa import MSA, Entropy, PercentageOfNonGaps, PercentageOfTotallyConservedColumns, Star, SumOfPairs
from pymsa import PAM250, Blosum62, FileMatrix
from pymsa.util.fasta import print_alignment


def run_all_scores(sequences: list) :
    my_list = []
    aligned_sequences = list(pair[1] for pair in sequences)
    sequences_id = list(pair[0] for pair in sequences)

    msa = MSA(aligned_sequences, sequences_id)
    #print_alignment(msa)

     # Entropy
    value = Entropy(msa).compute()
    my_list.append(value)    
    #print("Entropy score: {0}".format(value))
    return my_list



# Read sequences from the provided file
fileFasta =sys.argv[1]

with open(fileFasta, "r") as input_handle:
    # Utiliser une list comprehension pour accélérer la lecture des séquences
    sequences = [[seq_record.id, str(seq_record.seq)] for seq_record in SeqIO.parse(input_handle, "fasta")]

# Calcul des scores
my_list = run_all_scores(sequences)

print(str(my_list[0]))





