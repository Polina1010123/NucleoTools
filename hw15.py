from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


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

    return result


filtered_results = filter_fastq()
print(filtered_results)
