def filter_fastq(
    input_fastq: str,
    quality_threshold: float,
    gc_bounds: tuple[float, float] | float,
    length_bounds: tuple[int, int] | int,
    output_fastq: str,
) -> dict[str, tuple[str, str]]:
    non_filtered_sequences: dict[str, tuple[str, str]] = {}

    with open(input_fastq) as file:
        for i, line in enumerate(file):
            line = line.strip()

            if i % 4 == 0:
                seq_name = line[1:]
            elif i % 4 == 1:
                sequence = line
            elif i % 4 == 3:
                quality = line
                non_filtered_sequences[seq_name] = [sequence, quality]

    filtered_by_quality = filter_quality(non_filtered_sequences, quality_threshold)
    filtered_by_gc = filter_gc_content(filtered_by_quality, gc_bounds)
    filtered_results = filter_length(filtered_by_gc, length_bounds)

    return filtered_results


def filter_quality(
    non_filtered_sequences: dict[str, tuple[str, str]], quality_threshold: float
) -> dict[str, tuple[str, str]]:
    filtered_seqs: dict[str, tuple[str, str]] = {}

    for name, (sequence, quality) in non_filtered_sequences.items():
        quality_num = [ord(i) - 33 for i in quality]
        mean_quality = sum(quality_num) / len(quality_num)

        if mean_quality >= quality_threshold:
            filtered_seqs[name] = (sequence, quality)

    return filtered_seqs


def filter_gc_content(
    filtered_by_quality: dict[str, tuple[str, str]],
    gc_bounds: tuple[float, float] | float,
) -> dict[str, tuple[str, str]]:
    filtered_seqs_gc: dict[str, tuple[str, str]] = {}

    for name, (sequence, quality) in filtered_by_quality.items():
        gc_count = (
            sequence.count("G")
            + sequence.count("C")
            + sequence.count("g")
            + sequence.count("c")
        )
        gc_content = (gc_count / len(sequence)) * 100

        if isinstance(gc_bounds, tuple):
            if gc_bounds[0] <= gc_content <= gc_bounds[1]:
                filtered_seqs_gc[name] = (sequence, quality)
        else:
            if gc_content <= gc_bounds:
                filtered_seqs_gc[name] = (sequence, quality)

    return filtered_seqs_gc


def filter_length(
    filtered_by_gc: dict[str, tuple[str, str]], length_bounds: tuple[int, int] | int
) -> dict[str, tuple[str, str]]:
    filtered_seqs_gc_len: dict[str, tuple[str, str]] = {}

    for name, (sequence, quality) in filtered_by_gc.items():
        length_read = len(sequence)
        if isinstance(length_bounds, tuple):
            if length_bounds[0] <= length_read <= length_bounds[1]:
                filtered_seqs_gc_len[name] = (sequence, quality)
        else:
            if length_read <= length_bounds:
                filtered_seqs_gc_len[name] = (sequence, quality)

    return filtered_seqs_gc_len


def save_filtered_sequences(
    filtered_sequences: dict[str, tuple[str, str]], output_fastq: str
) -> None:
    with open(output_fastq, mode="w") as file:
        for name, (sequence, quality) in filtered_sequences.items():
            file.write(f"@{name}\n")
            file.write(f"{sequence}\n")
            file.write("+\n")
            file.write(f"{quality}\n")
