import numpy as np
import biotite.sequence as seq
import biotite.application.mafft as mafft
import biotite.sequence.align as align

alphabet = seq.Alphabet(("foo", "bar", 42))
sequences = [seq.GeneralSequence(alphabet, sequence) for sequence in [
    ["foo", "bar", 42, "foo", "foo", 42, 42],
    ["foo", 42, "foo", "bar", "foo", 42, 42],
]]
matrix = align.SubstitutionMatrix(
    alphabet, alphabet, np.array([
        [ 100, -100, -100],
        [-100,  100, -100],
        [-100, -100,  100]
    ])
)
alignment = mafft.MafftApp.align(sequences, bin_path='C:/Users/Rickman/Downloads/mafft-7.471-win64-signed/mafft-win/mafft.bat', matrix=matrix)
# As the alphabet do not has characters as symbols
# the alignment cannot be directly printed
# However, we can print the trace
print(alignment.trace)