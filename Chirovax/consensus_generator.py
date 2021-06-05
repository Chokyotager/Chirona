from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO

import numpy as np
import json

weight_distribution = json.load(open("weight_distribution.json"))
blosum62 = json.load(open("BLOSUM62.json"))

alignment = AlignIO.read("../Mutations/pangolin_output.aln.fasta", "fasta")
nucleotide_alignment = AlignIO.read("../Mutations/nucleotide_pangolin_spike.aln.fasta", "fasta")

original_sequence = SeqIO.read("../Mutations/spike_protein_COVID.fasta", "fasta")
alignment_size = len(alignment[0])

probabilities = list()
highlighted = list()

for i in range(alignment_size):

    residue_probabilities = dict()
    highlighted_probabilities = dict()

    for sequence in alignment:

        residue = sequence[i]
        variant = sequence.description

        if residue == "-":
            # BLOSUM differences
            residue = "*"

        highlighted_score = 0

        for j in range(len(weight_distribution["variants"])):

            if weight_distribution["variants"][j]["pango"] == variant:

                highlighted_score = weight_distribution["variants"][j]["score"]
                break

        if residue not in residue_probabilities.keys():

            residue_probabilities[residue] = int()
            highlighted_probabilities[residue] = int()

        residue_probabilities[residue] += 1
        highlighted_probabilities[residue] += highlighted_score

    highlighted.append(highlighted_probabilities)
    probabilities.append(residue_probabilities)

possible_residues = list(blosum62.keys())
datafile = [["#"] + possible_residues]

for i in range(len(probabilities)):

    residues = list(probabilities[i].keys())
    index_probabilities = [0] * len(possible_residues)

    residue_count = sum(probabilities[i].values())
    highlighted_count = sum(highlighted[i].values())

    for j in range(len(residues)):

        if residues[j] not in possible_residues:
            continue

        index = possible_residues.index(residues[j])

        # Include highlighted possibilities
        index_probabilities[index] += ((probabilities[i][residues[j]] / residue_count) * (1 - weight_distribution["highlighted-weightage"]) + (highlighted[i][residues[j]] / highlighted_count) * weight_distribution["highlighted-weightage"]) * 100

    datafile.append([i + 1] + index_probabilities)

open("Sequence matrix adjusted probability.tsv", "w+").write("\n".join(["\t".join([str(y) for y in x]) for x in datafile]))
