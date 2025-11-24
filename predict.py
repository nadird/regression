import sys
import os
from tensorflow import keras
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import os
import numpy as np
from Bio import SeqIO

from itertools import combinations
import csv
import statistics
from Bio import Align
import glob


def sequence_similarity(seq1, seq2):
    # Initialiser l'objet PairwiseAligner
    aligner = Align.PairwiseAligner()
    # Effectuer l'alignement entre les deux séquences
    alignments = aligner.align(seq1, seq2)
    # Sélectionner le meilleur alignement
    best_alignment = alignments[0]
    # Calculer le degré de similarité entre les séquences
    similarity = best_alignment.score / max(len(seq1), len(seq2))
    return similarity



# Ensure a file path is provided
if len(sys.argv) < 1:
   # print("Usage: python msascorecalculate.py <sequence_file>")
    sys.exit(1)

# Read sequences from the provided file
sequence_file = sys.argv[1]
try:
    with open(sequence_file, "r") as f:
        sequences = [line.strip() for line in f.readlines() if line.strip()]  # Remove empty lines
except Exception as e:
    #print(f"Error reading file: {e}")
    sys.exit(1)

# Ensure there are sequences before processing
if not sequences:
    #print("No sequences found in file.")
    sys.exit(1)



lengths=[]
N=len(sequences)
for record in sequences:
    lengths.append(len(record))
MLen=np.mean(lengths)
StdLen=np.std(lengths)

# Calculate distances between sequences
pairs = combinations(sequences, 2)
listx=[]
for pair in pairs:
    seq1, seq2 = pair
    listx.append([sequence_similarity(seq1, seq2)])

MPw=np.mean(listx)
StdPw=np.std(listx)

#X_enc=str(N)+','+str(MLen)+','+str(StdLen)+','+str(MPw)+','+str(StdPw)
X_enc = np.array([N, MLen, StdLen, MPw, StdPw])
#print(X_enc)

#print (sequences)
model1 = keras.models.load_model('Model_ANN-Stat.keras')
#print("Model loaded")
prediction1=model1.predict(X_enc.reshape(1,5))


#print("- Entropy Score                         : ",prediction1)
#print("- % of Totally Conserved Colummns Score : ",prediction2)


print("CodeBegin",prediction1)


