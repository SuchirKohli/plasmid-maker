from pathlib import Path


def read_fasta(fasta_path: str) -> str:
    """
    Reads a FASTA file and returns the DNA sequence as a single uppercase string.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    str
        DNA sequence with no headers or whitespace.
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    sequence_parts = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            sequence_parts.append(line)

    sequence = "".join(sequence_parts).upper()

    if not sequence:
        raise ValueError("No sequence found in FASTA file")

    return sequence
