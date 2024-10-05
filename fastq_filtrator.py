import sys

sys.path.append("/home/polina/проверка_дз4/modules")
import dop
import dopskript


def filter_fastq(
    seqs, quality_threshold=0, gc_bounds=(0, 100), length_bounds=(0, 2**32)
):
    filtered_seqs_gc_len = dop.filter_fastq(
        seqs, quality_threshold, gc_bounds, length_bounds
    )
    return filtered_seqs_gc_len


def run_dna_rna_tools(seq, action):
    result = dopskript.run_dna_rna_tools(seq, action)
    return result


if __name__ == "__main__":
    seqs = {
        "@SRX079801": ("AAAAAA", "@@@@@@@@@"),
        "@SRX079802": ("GCGCGC", "FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEF"),
        "@SRX079803": ("TTTTTTT", "0000000"),
    }
    filtered_results = filter_fastq(seqs)
    print(filtered_results)
    result = run_dna_rna_tools("ATG", "transcribe")
    print(result)
