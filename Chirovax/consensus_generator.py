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

summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus(threshold=0)

# Multiply by the weights
for sequence in alignment:
    raw_sequence = str(sequence.seq)

    for i in range(len(raw_sequence)):

        raw_sequence[i]
        consensus[i]
