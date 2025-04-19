from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import argparse
import logging
import os

logging.basicConfig(
    filename="fastq_filter.log",
    filemode="w",
    level=logging.DEBUG,
)


def filter_fastq(
    input_fastq="example_fastq.fastq",
    quality_threshold=20,
    gc_bounds=(0, 100),
    length_bounds=(0, 100),
    output_fastq="filtered_sequences.fastq",
):

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    if gc_bounds[0] > gc_bounds[1] or length_bounds[0] > length_bounds[1]:
        logging.error(
            f"Incorrect bounds: gc_bounds={gc_bounds}, length_bounds={length_bounds}"
        )
        raise ValueError(
            f"Incorrect bounds: gc_bounds={gc_bounds}, length_bounds={length_bounds}"
        )

    result = {}

    good_reads = (
        rec
        for rec in SeqIO.parse(input_fastq, "fastq")
        if (
            sum(rec.letter_annotations["phred_quality"])
            / len(rec.letter_annotations["phred_quality"])
        )
        >= quality_threshold
        and gc_bounds[0] <= gc_fraction(rec.seq) * 100 <= gc_bounds[1]
        and length_bounds[0] <= len(rec.seq) <= length_bounds[1]
    )

    filtered_records = list(good_reads)

    SeqIO.write(filtered_records, output_fastq, "fastq")

    for record in filtered_records:
        seq = record.seq
        quality_str = "".join(
            chr(q + 33) for q in record.letter_annotations["phred_quality"]
        )
        result[record.id] = (str(seq), quality_str)

    logging.info(
        f"The sequences that meet the specified length, quality, and GC content criteria, in the amount of {len(filtered_records)} can be found in {os.path.basename(output_fastq)}"
    )

    return result


parser = argparse.ArgumentParser(
    prog="Fastq filter",
    description="This tool is needed to filter DNA sequences based on several parameters: average read quality, GC content of the read, Length of the read",
    epilog='Usage example: python fastq_filtrator.py --input_fastq "C:/Users/kozlo/OneDrive/Рабочий стол/ИБ/Питон/HW18_инструменты_разработчика/example_fastq.fastq" --quality 25 --gc_bounds 40 60 --length_bounds 50 150 --output_fastq "C:/Users/kozlo/OneDrive/Рабочий стол/ИБ/Питон/HW18_инструменты_разработчика/filtered_sequences.fastq"',
)


parser.add_argument(
    "--input_fastq",
    default="example_fastq.fastq",
    help="Path to the input FASTQ file",
)
parser.add_argument(
    "--quality", type=float, default=20, help="Minimum average quality threshold"
)
parser.add_argument(
    "--gc_bounds",
    nargs=2,
    type=float,
    default=[0, 100],
    help="GC content bounds, e.g. --gc_bounds 40 60",
)
parser.add_argument(
    "--length_bounds",
    nargs=2,
    type=int,
    default=[0, 100],
    help="Length bounds, e.g. --length_bounds 50 150",
)
parser.add_argument(
    "--output_fastq",
    default="filtered_sequences.fastq",
    help="Output FASTQ file path",
)


if __name__ == "__main__":
    args = parser.parse_args()
    print(args)


filtered_results = filter_fastq()
print(filtered_results)
