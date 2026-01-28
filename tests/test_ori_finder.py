import sys
from pathlib import Path

# Ensure project root is on PYTHONPATH
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.fasta_utils import read_fasta
from src.ori_finder import (
    compute_gc_skew,
    cumulative_gc_skew,
    find_ori_position,
    extract_ori_region,
    find_clumps,
    find_ori_with_clumps,
)


FASTA_FILE = ROOT / "data" / "input" / "pUC19.fa"


def test_compute_gc_skew_runs():
    seq = read_fasta(FASTA_FILE)
    skew = compute_gc_skew(seq, window_size=500)

    assert isinstance(skew, list)
    assert len(skew) > 0
    assert all(isinstance(x, (int, float)) for x in skew)


def test_cumulative_gc_skew_monotonic_length():
    skew = [1, -2, 3, -1]
    cumulative = cumulative_gc_skew(skew)

    assert len(cumulative) == len(skew)
    assert cumulative[-1] == sum(skew)


def test_find_ori_position_returns_valid_indices():
    seq = read_fasta(FASTA_FILE)
    start, end = find_ori_position(seq, window_size=500)

    assert isinstance(start, int)
    assert isinstance(end, int)
    assert 0 <= start < len(seq)
    assert start < end <= len(seq)


def test_extract_ori_region_returns_sequence_and_coords():
    seq = read_fasta(FASTA_FILE)
    ori_seq, start, end = extract_ori_region(seq, window_size=500, flank=2000)

    assert isinstance(ori_seq, str)
    assert len(ori_seq) > 0
    assert 0 <= start < end <= len(seq)
    assert len(ori_seq) == end - start


def test_find_clumps_runs_on_small_sequence():
    test_seq = "ATG" * 400  # repetitive sequence
    clumps = find_clumps(
        test_seq,
        k=3,
        window_length=100,
        min_occurrence=5,
    )

    assert isinstance(clumps, dict)

    # If clumps exist, ensure structure is correct
    for (start, end), kmers in clumps.items():
        assert isinstance(start, int)
        assert isinstance(end, int)
        assert isinstance(kmers, set)
        assert all(len(k) == 3 for k in kmers)


def test_find_ori_with_clumps_integration():
    seq = read_fasta(FASTA_FILE)
    result = find_ori_with_clumps(
        seq,
        skew_window=500,
        flank=2000,
        k=9,
        clump_window=500,
        min_occurrence=3,
    )

    assert isinstance(result, dict)

    assert "ori_start" in result
    assert "ori_end" in result
    assert "ori_sequence" in result
    assert "clumps" in result

    assert isinstance(result["ori_sequence"], str)
    assert len(result["ori_sequence"]) > 0
    assert isinstance(result["clumps"], dict)
