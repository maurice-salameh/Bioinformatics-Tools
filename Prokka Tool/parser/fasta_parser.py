from Bio import SeqIO

def parse_fasta(fasta_file, min_length=200):
    """
    Parse the FASTA file and return sequences longer than min_length.
    """
    sequences = [rec for rec in SeqIO.parse(fasta_file, "fasta") if len(rec) >= min_length]
    print(f"Loaded {len(sequences)} sequences longer than {min_length} bp.")
    return sequences
