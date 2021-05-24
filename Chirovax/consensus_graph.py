from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO

import numpy as np
import json
from matplotlib import pyplot as plt

weight_distribution = json.load(open("weight_distribution.json"))
blosum62 = json.load(open("BLOSUM62.json"))

alignment = AlignIO.read("../Mutations/pangolin_output.aln.fasta", "fasta")
nucleotide_alignment = AlignIO.read("../Mutations/nucleotide_pangolin_spike.aln.fasta", "fasta")

original_sequence = SeqIO.read("../Mutations/spike_protein_COVID.fasta", "fasta")

summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus(threshold=0)
pssm = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore=None)

nucleotide_summary_align = AlignInfo.SummaryInfo(nucleotide_alignment)
nucleotide_consensus = nucleotide_summary_align.dumb_consensus(threshold=0)
nucleotide_pssm = nucleotide_summary_align.pos_specific_score_matrix(nucleotide_consensus, chars_to_ignore=None)

graph = list()

for amino_acid in pssm:

    graph.append(max(amino_acid.values()) / len(alignment))

figure = plt.figure()
plt.plot(graph)
plt.ylim(0.95, 1)
figure.savefig("Frequencies unweighted.png")

graph = list()

for nucleotide in nucleotide_pssm:

    graph.append(max(nucleotide.values()) / len(alignment))

figure = plt.figure()
plt.plot(graph)
plt.ylim(0.95, 1)
figure.savefig("Frequencies unweighted nucleotide.png")
