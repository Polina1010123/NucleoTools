def filter_fastq(seqs, quality_threshold, gc_bounds, length_bounds):
    filtered_seqs = {}

    # Первый фильтр
    for name, (sequence, quality) in seqs.items():
        quality_num = [ord(i) - 33 for i in quality]
        mean_quality = sum(quality_num) / len(quality_num)

        if mean_quality >= quality_threshold:
            filtered_seqs[name] = (sequence, quality)

    # Второй фильтр
    filtered_seqs_gc = {}
    for name, (sequence, quality) in filtered_seqs.items():
        gc_count = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
        gc_content = (gc_count / len(sequence)) * 100

        if isinstance(gc_bounds, tuple):
            if gc_bounds[0] <= gc_content <= gc_bounds[1]:
                filtered_seqs_gc[name] = (sequence, quality)
        else:
            if gc_content <= gc_bounds:
                filtered_seqs_gc[name] = (sequence, quality)

    # Третий фильтр
    filtered_seqs_gc_len = {}
    for name, (sequence, quality) in filtered_seqs_gc.items():
        length_read = len(sequence)
        if isinstance(length_bounds, tuple):
            if length_bounds[0] <= length_read <= length_bounds[1]:
                filtered_seqs_gc_len[name] = (sequence, quality)
        else:
            if length_read <= length_bounds:
                filtered_seqs_gc_len[name] = (sequence, quality)

    return filtered_seqs_gc_len
