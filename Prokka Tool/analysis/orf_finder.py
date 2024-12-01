def find_orfs(sequence, min_length=300):
    """
    Find open reading frames (ORFs) in a DNA sequence with a minimum length.
    """
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    sequence_len = len(sequence)

    for i in range(sequence_len - 3):
        if sequence[i:i+3] == start_codon:  # Check for start codon
            for j in range(i + 3, sequence_len - 3, 3):
                if sequence[j:j+3] in stop_codons:  # Check for stop codon
                    orf = sequence[i:j+3]
                    if len(orf) >= min_length:
                        orfs.append((i, j+3, orf))
                    break
    return orfs
