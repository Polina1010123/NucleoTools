import sys

sys.path.append("")

import dop2
import dopskript


def filter_fastq(
    input_fastq: str = "c:/Users/kozlo/Downloads/example_fastq.fastq",
    quality_threshold: float = 20,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[float, float] = (0, 100),
    output_fastq: str = "c:/Users/kozlo/Downloads/filtered/filtered_sequences.fastq",
) -> dict[str, tuple[str, str]]:
    """
    Action:
    filtering DNA read sequences (input_fastq argument) based on specified parameters.

    Parameters:

    path to the input file with unfiltered sequences (input_fastq argument);
    average quality of the read (quality_threshold argument);
    GC content of the read (gc_bounds argument);
    length of the read (length_bounds argument);
    path to the output file with filtered sequences (output_fastq argument).

    Result:
    a file containing filtered sequences that meet the quality, GC content, and length criteria.

    """

    filtered_sequences = dop2.filter_fastq(
        input_fastq, quality_threshold, gc_bounds, length_bounds, output_fastq
    )
    dop2.save_filtered_sequences(filtered_sequences, output_fastq)
    return filtered_sequences


def run_dna_rna_tools(seq: str, action: str) -> str:
    """
    Action:
    transformation of a DNA or RNA sequence (seq argument).

    Types of Transformation:
    converting a coding DNA strand to mRNA (transcribe argument);
    writing the sequence in reverse order (reverse argument);
    creating a complementary sequence (complement argument);
    writing the complementary sequence in reverse order (reverse_complement argument).

    Returns:
    The result of the specified action performed on the sequence.

    """

    result: str = dopskript.run_dna_rna_tools(seq, action)
    return result


if __name__ == "__main__":

    filtered_results = filter_fastq()
    print(filtered_results)

    result = run_dna_rna_tools("ATGAAAA", "reverse")
    print(result)
