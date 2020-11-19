import numpy as np
import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.io.genbank as gb
import biotite.sequence.graphics as graphics
import biotite.application.clustalo as clustalo
import biotite.database.entrez as entrez

"""def _get_entropy(sequences):
    alphabet = sequences[0].get_alphabet()
    freq = np.zeros((len(sequences[0].code), len(alphabet)))

    for i in range(len(sequences[0].code)):
        for j in range(len(sequences)):
            code = sequences[j].code[i]
            freq[i, code] += 1
    freq = freq / np.sum(freq, axis=1)[:, np.newaxis]
    # 0 * log2(0) = 0 -> Convert NaN to 0
    no_zeros = freq != 0
    pre_entropies = np.zeros(freq.shape)
    pre_entropies[no_zeros] \
        = freq[no_zeros] * np.log2(freq[no_zeros])
    entropies = -np.sum(pre_entropies, axis=1)
    max_entropy = np.log2(len(alphabet))
    return freq, entropies, max_entropy"""

alignment = fasta.FastaFile.read(r"C:\Users\Rickman\Labrotation\Sequences.txt")
alignment = fasta.get_alignment(alignment)


sequences = []
for seq_str in alignment.get_gapped_sequences():
    sequences.append(seq.ProteinSequence(seq_str))
    print(seq_str)


#print(_get_entropy(sequences))

fig = plt.figure(figsize=(8.0, 3.0))
ax = fig.add_subplot(111)
graphics.plot_sequence_logo(ax, alignment)
ax.set_xticks([5,10,15,20])
ax.set_xlabel("Residue position")
ax.set_ylabel("Bits")
# Only show left and bottom spine
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
fig.tight_layout()
# sphinx_gallery_thumbnail_number = 2

plt.show()