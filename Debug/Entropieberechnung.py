import numpy as np
import timeit
import csv
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
from biotite.sequence.alphabet import LetterAlphabet
"""import matplotlib.pyplot as plt
from biotite.visualize import set_font_size_in_coord
from biotite.sequence.graphics.colorschemes import get_color_scheme"""


def _get_entropy(alignment):

    alphabet = alignment.sequences[0].get_alphabet()
    trace = alignment.trace
    sequences = alignment.sequences
    freq = np.zeros((len(trace), len(alphabet)))
    for i in range(trace.shape[0]):
        for j in range(trace.shape[1]):
            index = trace[i,j]
            if index != -1:
                code = sequences[j].code[index]
                freq[i, code] += 1
    freq = freq / np.sum(freq, axis=1)[:, np.newaxis]
    no_zeros = freq != 0
    pre_entropies = np.zeros(freq.shape)
    pre_entropies[no_zeros] \
        = freq[no_zeros] * np.log2(freq[no_zeros])
    entropies = -np.sum(pre_entropies, axis=1)
    max_entropy = np.log2(len(alphabet))

    with open('Entropies.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['Residue', 'Entropy'])
        count = 1
        for entropy in entropies:
            writer.writerow([count, entropy])
            count += 1

    end = timeit.default_timer()
    print(end - start)

    return freq, entropies, max_entropy

i = 0
while i <= 0:

    start = timeit.default_timer()

    alignment = fasta.FastaFile.read(r"C:\Users\Rickman\Documents\GitHub\Labrotation\Testsequenzen\alan_sequences_aligned_squashed.txt")
    alignment = fasta.get_alignment(alignment)

    frequencies, entropies, max_entropy = _get_entropy(alignment)

    """fig = plt.figure(figsize=(8.0, 3.0))
    ax = fig.add_subplot(111)
    plot_sequence_logo(ax, alignment)"""

    i += 1


"""ax.set_xticks([5,10,15,20])
ax.set_xlabel("Residue position")
ax.set_ylabel("Bits")
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
fig.tight_layout()

plt.show()

def plot_sequence_logo(axes, alignment, scheme=None, **kwargs):

    from matplotlib.text import Text

    sequences = alignment.sequences
    alphabet = sequences[0].get_alphabet()
    for seq in sequences:
        if seq.get_alphabet() != alphabet:
            raise ValueError("Alphabets of the sequences in the alignment "
                             "are not equal")
    if not isinstance(alphabet, LetterAlphabet):
        raise TypeError("The sequences' alphabet must be a letter alphabet")

    if scheme is None:
        colors = get_color_scheme("buried", alphabet)
    elif isinstance(scheme, str):
        colors = get_color_scheme(scheme, alphabet)
    else:
        colors = scheme
    
    kwargs.pop("color", None)
    kwargs.pop("size",  None)
    
    frequencies, entropies, max_entropy = _get_entropy(alignment)
    
    stack_heights = (max_entropy - entropies)
    symbols_heights = stack_heights[:, np.newaxis] * frequencies
    index_order = np.argsort(symbols_heights, axis=1)
    for i in range(symbols_heights.shape[0]):
        index_order = np.argsort(symbols_heights)
        start_height = 0
        for j in index_order[i]:
            height = symbols_heights[i,j]
            if height > 0:
                symbol = alphabet.decode(j)
                text = axes.text(
                    i+0.5, start_height, symbol,
                    ha="left", va="bottom", color=colors[j],
                    size=1,
                    **kwargs
                )
                text.set_clip_on(True)
                set_font_size_in_coord(text, width=1, height=height)
                start_height += height
    
    axes.set_xlim(0.5, len(alignment)+0.5)
    axes.set_ylim(0, max_entropy)"""