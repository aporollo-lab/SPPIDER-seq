import argparse


DEFAULT_FASTA = """
>Example1
MTMTLHTKASGMALLHQIQGNELEPLNRPQLKIPLERPLGEVYLDSSKPAVYNYPEGAAY
>Example2
MVEIFDMLLATSSRFRMMNLQGEEFVCLKSIIL
"""


def parse_fasta(text):
    sequences = []
    header = None
    seq_lines = []
    for line in text.strip().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header and seq_lines:
                sequences.append((header, "".join(seq_lines)))
            header = line
            seq_lines = []
        else:
            seq_lines.append(line)
    if header and seq_lines:
        sequences.append((header, "".join(seq_lines)))
    return sequences


def generate_alanine_scans(fasta_text, window_size, stride):
    if window_size < 2:
        raise ValueError("Window size must be at least 2.")
    if stride < 1:
        raise ValueError("Stride must be at least 1.")

    parsed = parse_fasta(fasta_text)
    result = []
    for header, sequence in parsed:
        seq_len = len(sequence)
        for i in range(0, seq_len, stride):
            remaining = seq_len - i
            if remaining <= 0:
                break
            actual_window = min(window_size, remaining)
            modified = sequence[:i] + "A" * actual_window + sequence[i + actual_window:]
            result.append(f"{header}_AlaScan{i+1}\n{modified}")
            if actual_window < window_size:
                break  # Stop after first trimming event
    return "\n".join(result)


def read_input(path):
    if path:
        with open(path) as f:
            return f.read()
    return DEFAULT_FASTA


def main():
    parser = argparse.ArgumentParser(
        description="Generate sliding-window alanine-scanning FASTA variants."
    )
    parser.add_argument("--input", help="Input FASTA file.")
    parser.add_argument(
        "--output",
        default="alanine_scanned_sequences.fasta",
        help="Output FASTA file. Default: alanine_scanned_sequences.fasta",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=10,
        help="Alanine replacement window size. Default: 10",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=3,
        help="Sliding-window stride. Default: 3",
    )
    args = parser.parse_args()

    scanned = generate_alanine_scans(read_input(args.input), args.window_size, args.stride)
    with open(args.output, "w") as f:
        f.write(scanned)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
