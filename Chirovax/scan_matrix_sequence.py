import json

blosum62 = json.load(open("BLOSUM62.json"))

data = [x.split("\t") for x in open("Sequence matrix adjusted probability.tsv").read().split("\n") if len(x.split("\t")) > 1]

index = data[0][1:]
incidence_cutoff = 2

mutations = list()

# Required for residue number changes
skips = int()

for i in range(1, len(data)):

    residue_number = int(data[i][0])

    residue_probabilities = [float(x) for x in data[i][1:]]
    highest_probability = max(residue_probabilities)

    consensus_resiude = index[residue_probabilities.index(highest_probability)]

    if consensus_resiude == "*":
        skips += 1

    for j in range(len(residue_probabilities)):

        probability = residue_probabilities[j]

        if probability < highest_probability and probability >= incidence_cutoff:

            mutated_residue = index[j]
            mutation = consensus_resiude + str(residue_number - skips) + mutated_residue
            blosum62_value = blosum62[consensus_resiude][mutated_residue]

            blosum62_original = blosum62[consensus_resiude][consensus_resiude]
            blosum62_delta = blosum62_original - blosum62_value

            mutations.append({"mutation": mutation, "probability": probability, "score": probability * blosum62_delta, "blosum62": blosum62_delta, "position": residue_number})

            #print(consensus_resiude + residue_number + index[j] + "     [" + str(probability) + "%]")

mutations = sorted(mutations, key=lambda x: -x["score"])
writable = ["Mutation\tAlignment position\tScore\tAdjusted incidence\tBLOSUM62 delta"]

for i in range(len(mutations)):

    print("{}    [Score: {}] [Adjusted incidence: {}%] [BLOSUM62 delta: {}]".format(mutations[i]["mutation"], mutations[i]["score"], mutations[i]["probability"], mutations[i]["blosum62"]))
    writable.append("\t".join([mutations[i]["mutation"], str(mutations[i]["position"]), str(mutations[i]["score"]), str(mutations[i]["probability"]), str(mutations[i]["blosum62"])]))

open("mutations_of_interest.tsv", "w+").write("\n".join(writable))
