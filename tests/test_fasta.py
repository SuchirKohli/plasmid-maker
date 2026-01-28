import sys
from pathlib import Path

# Add project root to PYTHONPATH
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.fasta_utils import read_fasta


def test_read_fasta_returns_sequence():
    seq = read_fasta("data/input/pUC19.fa")
    assert isinstance(seq, str)
    assert len(seq) > 0


def test_read_fasta_contains_only_dna():
    seq = read_fasta("data/input/pUC19.fa")
    valid_bases = set("ACGTN")
    assert set(seq).issubset(valid_bases)
