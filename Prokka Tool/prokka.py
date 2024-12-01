from parser.fasta_parser import parse_fasta
from analysis.orf_finder import find_orfs
from analysis.orf_filter import filter_overlapping_orfs
from output.gff_writer import write_gff

def main(fasta_file, gff_output):
    sequences = parse_fasta(fasta_file)
    all_orfs = []
    for record in sequences:
        print(f"Processing {record.id}...")
        orfs = find_orfs(str(record.seq))
        filtered_orfs = filter_overlapping_orfs(orfs)
        all_orfs.append(filtered_orfs)
    write_gff(gff_output, sequences, all_orfs)

if __name__ == "__main__":
    input_fasta = "data/SRR11528307.fasta"  # Replace with your FASTA file
    output_gff = "data/annotations.gff3"
    main(input_fasta, output_gff)
