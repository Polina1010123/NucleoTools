import sys

sys.path.append("/home/polina/проверка_дз4/проверка2/modules")

import dop
import dopskript


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    quality_threshold: float = 0,
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: float | tuple[float, float] = (0, 2**32),
) -> dict[str, tuple[str, str]]:
    """
    Действие:
    фильтрование последовательности ридов ДНК (аргумент seqs)
    по заданным параметрам.

    Параметры:
    среднее качества рида (аргумент quality_threshold)
    GC-состав рида (аргумент gc_bounds)
    длина рида (аргумент length_bounds)

    Возвращает:
    отфильтрованные последовательности, которые удовлетворяют критериям
    качества, GC-содержания и длины.

    """

    filtered_seqs_gc_len: dict[str, tuple[str, str]] = dop.filter_fastq(
        seqs, quality_threshold, gc_bounds, length_bounds
    )

    return filtered_seqs_gc_len


def run_dna_rna_tools(seq: str, action: str) -> str:
    """
    Действие:
    преобразование последовательности ДНК или РНК (аргумент seq).

    Виды преобразования:
    превращение кодирующей цепи ДНК в мРНК (аргумент transcribe)
    записывание последовательности в обратном порядке (аргумент reverse)
    создание комплементарной последовательности (аргумент complement)
    записывание комплеменентарной последовательности в обратном порядке
    (аргумент reverse_complement)

    Возвращает:
    результат выполнения указанного действия над последовательностью.

    """

    result: str = dopskript.run_dna_rna_tools(seq, action)
    return result


if __name__ == "__main__":
    seqs = {
        "@SRX079801": ("AAAAAA", "@@@@@@@@@"),
        "@SRX079802": ("GCGCGC", "FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEF"),
        "@SRX079803": ("TTTTTTT", "0000000"),
    }

    filtered_results = filter_fastq(seqs)
    print(filtered_results)

    result = run_dna_rna_tools("ATGAAAA", "reverse")
    print(result)
