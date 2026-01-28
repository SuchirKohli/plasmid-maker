from typing import List, Tuple, Dict, Set
from collections import defaultdict

def compute_gc_skew(sequence: str, window_size: int = 1000) -> List[float]:
    """
    Compute GC skew over sliding windows.

    GC skew = (G - C)

    Parameters
    ----------
    sequence : str
        DNA sequence
    window_size : int
        Window size for skew calculation

    Returns
    -------
    List[float]
        GC skew values for each window
    """
    skew_values = []

    for i in range(0, len(sequence) - window_size + 1, window_size):
        window = sequence[i : i + window_size]
        g = window.count("G")
        c = window.count("C")

        if g + c == 0:
            skew = 0.0
        else:
            skew = (g - c) 

        skew_values.append(skew)

    return skew_values

def cumulative_gc_skew(skew_values: List[float]) -> List[float]:
    """
    Compute cumulative GC skew.

    Parameters
    ----------
    skew_values : List[float]

    Returns
    -------
    List[float]
        Cumulative GC skew
    """
    cumulative = []
    total = 0.0

    for value in skew_values:
        total += value
        cumulative.append(total)

    return cumulative

def find_ori_position(sequence: str, window_size: int = 1000) -> Tuple[int, int]:
    """
    Find ORI position based on minimum cumulative GC skew.

    Returns
    -------
    Tuple[int, int]
        (ori_start, ori_end)
    """
    skew = compute_gc_skew(sequence, window_size)
    cumulative = cumulative_gc_skew(skew)

    min_index = cumulative.index(min(cumulative))

    ori_start = min_index * window_size
    ori_end = ori_start + window_size

    return ori_start, ori_end

def extract_ori_region(sequence: str, window_size: int = 1000, flank: int = 5000,) -> Tuple[str, int, int]:
    """
    Extract a candidate ORI region based on GC skew.

    Returns
    -------
    Tuple[str, int, int]
        (ori_sequence, start, end)
    """
    ori_start, ori_end = find_ori_position(sequence, window_size)
    ori_center = ori_start + window_size // 2

    start = max(0, ori_center - flank)
    end = min(len(sequence), ori_center + flank)

    return sequence[start:end], start, end

def find_clumps(sequence: str, k: int = 9, window_length: int = 1000, min_occurrence: int = 3,) -> Dict[Tuple[int, int], Set[str]]:
    """
    Find frequent k-mer clumps in a sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence
    k : int
        k-mer size
    window_length : int
        Window size for clump finding
    min_occurrence : int
        Minimum frequency to qualify as a clump

    Returns
    -------
    Dict[(start, end), Set[str]]
        Window positions mapped to clumped k-mers
    """
    clumps = {}

    for i in range(0, len(sequence) - window_length + 1):
        window = sequence[i : i + window_length]
        kmer_counts = defaultdict(int)

        for j in range(len(window) - k + 1):
            kmer = window[j : j + k]
            kmer_counts[kmer] += 1

        frequent_kmers = {
            kmer for kmer, count in kmer_counts.items()
            if count >= min_occurrence
        }

        if frequent_kmers:
            clumps[(i, i + window_length)] = frequent_kmers

    return clumps

def find_ori_with_clumps(
    sequence: str,
    skew_window: int = 1000,
    flank: int = 5000,
    k: int = 9,
    clump_window: int = 1000,
    min_occurrence: int = 3,
):
    """
    Identify ORI using GC skew and refine using clump finding.

    Returns
    -------
    dict
        ORI information including region and clumps
    """
    ori_region, start, end = extract_ori_region(
        sequence,
        window_size=skew_window,
        flank=flank,
    )

    clumps = find_clumps(
        ori_region,
        k=k,
        window_length=clump_window,
        min_occurrence=min_occurrence,
    )

    return {
        "ori_start": start,
        "ori_end": end,
        "ori_sequence": ori_region,
        "clumps": clumps,
    }
