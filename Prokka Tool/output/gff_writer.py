def write_gff(output_file, sequences, orfs):
    """
    Write the GFF3 file with ORF annotations.
    """
    with open(output_file, "w") as gff:
        gff.write("##gff-version 3\n")
        for record, record_orfs in zip(sequences, orfs):
            for i, (start, end, orf) in enumerate(record_orfs):
                product = f"hypothetical protein {i+1}"  # Example product
                inference = "ab initio prediction"  # Example inference
                gff.write(
                    f"{record.id}\tProkka-like\tgene\t{start+1}\t{end}\t.\t+\t.\t"
                    f"ID=gene{i+1};Name=ORF{i+1};Product={product};Inference={inference}\n"
                )
        gff.write("##FASTA\n")
        for record in sequences:
            gff.write(f">{record.id}\n{record.seq}\n")
    print(f"GFF3 written to {output_file}")
