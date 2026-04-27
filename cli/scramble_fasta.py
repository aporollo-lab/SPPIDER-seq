import argparse
import random


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


def scramble_sequence(seq):
    scrambled = list(seq)
    random.shuffle(scrambled)
    return "".join(scrambled)


def generate_scrambled_sequences(fasta_text, num_scrambles):
    parsed = parse_fasta(fasta_text)
    result = []
    for header, sequence in parsed:
        for i in range(num_scrambles):
            scrambled = scramble_sequence(sequence)
            result.append(f"{header}_scramble{i+1}\n{scrambled}")
    return "\n".join(result)


def read_input(path):
    if path:
        with open(path) as f:
            return f.read()
    return DEFAULT_FASTA


def main():
    parser = argparse.ArgumentParser(
        description="Scramble FASTA sequences while preserving amino acid composition."
    )
    parser.add_argument("--input", help="Input FASTA file.")
    parser.add_argument(
        "--output",
        default="scrambled_sequences.fasta",
        help="Output FASTA file. Default: scrambled_sequences.fasta",
    )
    parser.add_argument(
        "--num-scrambles",
        type=int,
        default=3,
        help="Scrambled variants per input sequence. Default: 3",
    )
    parser.add_argument("--seed", type=int, help="Optional random seed.")
    args = parser.parse_args()

    if args.num_scrambles <= 0:
        raise ValueError("Number of scrambles must be positive.")
    if args.seed is not None:
        random.seed(args.seed)

    scrambled = generate_scrambled_sequences(read_input(args.input), args.num_scrambles)
    with open(args.output, "w") as f:
        f.write(scrambled)

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
