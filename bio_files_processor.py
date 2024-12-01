def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str
) -> dict[str, str]:
    """
    Action: Reading the provided FASTA file, where the sequence (DNA/RNA/protein, etc.)
    may be split across multiple lines, and saving the data to a new FASTA file
    where each sequence is contained in a single line.

    Parameters:
    path to the input FASTA file containing multiline sequences (input_fasta argument);
    path to the output FASTA file where single-line sequences will be written (output_fasta argument);

    Result:
    a file containing sequences written in a single line.
    """

    sequences_oneline: dict[str, list[str]] = {}
    with open(input_fasta) as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                name = line
                sequences_oneline[name] = []
            else:
                seq = line
                sequences_oneline[name].append(seq)

    for name, seq in sequences_oneline.items():
        sequences_oneline[name] = "".join(seq)

    with open(output_fasta, mode="w") as file:
        for name, seq in sequences_oneline.items():
            file.write(f">{name}\n")
            file.write(f"{seq}\n")

    return sequences_oneline


input_fasta = "c:/Users/kozlo/Downloads/example_multiline_fasta.fasta"
output_fasta = "c:/Users/kozlo/Downloads/multiline_fasta_to_oneline/multiline_fasta_to_oneline.fastq"
result = convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
print(result)


def parse_blast_output(input_file: str, output_file: str) -> list[str]:
    """
    Action:
    extracting protein names from a text file with BLAST results that correspond to the best match
    and saving the names in alphabetical order to a new file.

    Parameters:

    path to the input text file with BLAST results (input_file argument);
    path to the output text file where the best matches will be written
    (output_file argument).

    Result:
    a file containing a list of proteins corresponding to the best matches with the database,
    sorted in alphabetical order.
    """
    proteins: list[str] = []
    with open(input_file) as file:
        for line in file:
            line = line.strip()
            if line.startswith("Sequences producing"):
                next(file)
                next(file)
                protein = next(file)
                protein = protein.split("  ")[0]
                proteins.append(protein)
                proteins = sorted(proteins, key=str.lower)
    with open(output_file, mode="w") as file:
        for protein in proteins:
            file.write(f">{protein}\n")
    return proteins


input_file = "c:/Users/kozlo/Downloads/example_blast_results.txt"
output_file = "c:/Users/kozlo/Downloads/parse_blast_output/parse_blast_output.txt"

result = parse_blast_output(input_file, output_file)
print(result)
