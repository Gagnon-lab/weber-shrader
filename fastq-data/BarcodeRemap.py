import pandas as pd
import re

def get_sequence(bc: str, ref: str) -> str:
    """
    Get the original sequence from HMID barcode. Applies edits to the reference sequence to obtain a 'clean' version
    of the original sequence that produced the barcode when run through GESTALT.

    :param bc: The HMID barcode
    :param ref: The reference sequence
    :return: The original sequence

    *Example:*\n
    >>>get_sequence('NONE_3D+10_NONE', 'AAACCCGGGTTTAAA')\n
    AAACCCGGGTAA
    """
    delimiters = r'[_&]'
    split = re.split(delimiters, bc)
    split_set = set(split)
    mutation_set = set()

    for edit in split_set:
        if edit == 'NONE':
            continue
        elif 'D' in edit:
            del_length, del_index = map(int, edit.split('D+'))
            mutation_set.add(Mutator(del_index, -del_length, lambda seq, idx, dl=del_length: seq[0:idx] + seq[idx + dl:]))
        elif 'I' in edit:
            ins_length, ins_index, ins_seq = edit.replace('I', '').split('+')
            mutation_set.add(Mutator(int(ins_index), int(ins_length), lambda seq, idx, iseq=ins_seq: seq[0:idx] + iseq + seq[idx:]))
        elif 'S' in edit:
            swp_length, swp_index, swp_seq = edit.replace('S', '').split('+')
            mutation_set.add(Mutator(int(swp_index), 0, lambda seq, idx, sl=int(swp_length), sseq=swp_seq: seq[0:idx] + sseq + seq[idx + sl:]))

    sorted_mutation_set = sorted(mutation_set, key=lambda m: m.start_index)
    index_offset = 0
    for mut in sorted_mutation_set:
        ref = mut.apply(ref, index_offset)
        index_offset += mut.offset_index
    return ref

class Mutator:
    """
    Helper class for applying mutations to a sequence.
    """
    def __init__(self, start_index, offset_index, function):
        self.start_index = start_index
        self.offset_index = offset_index
        self.function = function

    def apply(self, sequence, idx_offset):
        return self.function(sequence, self.start_index + idx_offset)

class BarcodeRemap:
    # The default reference sequence is the zebrafish reference sequence used in the GESTALT pipeline
    zebrafish_reference_sequence = 'AGATCTCGAGCTCaagcttcgGACAGCAGTATCATCGACTATGGagtcGAGAGCGCGCTCGTCGACTATGGagtcGTCAGCAGTACTACTGACGATGGagtcGACAGCAGTGTGTGAGTCTATGGagtcGAGAGCATAGACATCGAGTATGGagtcGACTACAGTCGCTACGACTATGGagtcGACAGAGATATCATGCAGTATGGagtcGACAGCAGTATCTGCTGTCATGGagtcGACTGCACGACAGTCGACTATGGagtcGCGAGCGCTATGAGCGACTATGGGCATGCaattaagtctcagccaccgttaactgaacatgtcaaaacc'

    def __init__(self, reference_sequence=zebrafish_reference_sequence):
        self.reference_sequence = reference_sequence

    def get_all_sequences(self, file_path):
        """
        Finds the 'original' sequences for a set of HMID barcodes derived from the output to the GESTALT pipeline (.allReadCounts or .topReadCounts files).
        Takes the path to a .allReadCounts or .topReadCounts file and returns a pandas dataframe representing the contents of the file, plus a new column
        named 'sequence'. The 'sequence' column contains the original DNA sequences for each barcode.

        :param file_path: The path to the .allReadCounts file
        :return: A pandas dataframe representing the .allReadCounts file with the 'sequence' column added.
        """
        read_count_df = pd.read_csv(file_path, sep='\t')
        read_count_df['sequence'] = read_count_df['event'].apply(lambda x: get_sequence(x, ref=self.reference_sequence))
        return read_count_df