def run_dna_rna_tools(*args):
    *seqs, action = args
    if action == "transcribe":
        result = transcribe(seqs)
        return "".join(result) if len(seqs) == 1 else result
    elif action == "reverse":
        result = reverse(seqs)
        return "".join(result) if len(seqs) == 1 else result
    elif action == "complement":
        result = complement(seqs)
        return "".join(result) if len(seqs) == 1 else result
    elif action == "reverse_complement":
        result = reverse_complement(seqs)
        return "".join(result) if len(seqs) == 1 else result

def transcribe(seqs):
    result = []
    for i in seqs:
        if "U" in i or "u" in i:
            result.append("Введите ДНК")
        else:
            result.append(i.replace("T", "U").replace("t", "u"))
    return result

def reverse(seqs):
    return [i[::-1] for i in seqs]

def complement(seqs):
    results = []
    for i in seqs:
        table = str.maketrans("ATGCatgc", "TACGtacg")
        results.append(i.translate(table))
    return results

def reverse_complement(seqs):
    results = []
    for i in seqs:
        table = str.maketrans("ATGCatgc", "TACGtacg")
        results.append(i.translate(table)[::-1])
    return results
