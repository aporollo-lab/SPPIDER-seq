#!/usr/bin/env python3

import argparse
import csv
import os
import re
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from io import TextIOWrapper


CUTOFFS = (0.7, 0.5, 0.3, 0.1)

MODEL_SUFFIXES = (
    "query_as_receptor",
    "query_as_peptide",
    "receptor_centric",
    "peptide_centric",
)


@dataclass
class Prediction:
    source_name: str
    query: str
    partner: str
    model: str
    positions: list[int]
    amino_acids: list[str]
    probabilities: list[float]


def sanitize_filename(text):
    return re.sub(r"[^a-zA-Z0-9_.-]", "_", text)


def normalize_model(model):
    model = model.strip()

    if model.endswith("-centric"):
        model = model[: -len("-centric")]

    model = model.replace("-", "_")

    model_aliases = {
        "query_as_receptor": "receptor_centric",
        "query_as_peptide": "peptide_centric",
        "receptor": "receptor_centric",
        "peptide": "peptide_centric",
    }

    return model_aliases.get(model, model)


def parse_header(line):
    if not line.startswith("#"):
        return None

    match = re.search(
        r"Query:(.*?)\s+Partner:(.*?)\s+Model:(.*)$",
        line.strip(),
    )

    if not match:
        return None

    return {
        "query": match.group(1).strip(),
        "partner": match.group(2).strip(),
        "model": normalize_model(match.group(3).strip()),
    }


def parse_filename(name):
    base = os.path.basename(name)

    if not base.endswith(".txt"):
        return None

    stem = base[:-4]

    for model in MODEL_SUFFIXES:
        suffix = f"__{model}"

        if stem.endswith(suffix):
            pair_stem = stem[: -len(suffix)]

            if "__" not in pair_stem:
                return None

            query, partner = pair_stem.rsplit("__", 1)

            return {
                "query": query,
                "partner": partner,
                "model": model,
            }

    return None


# ----------------------------------------------------------------------
# UniProt -> gene mapping
# ----------------------------------------------------------------------

def detect_mapping_delimiter(mapping_file):
    """
    Determine whether a mapping file is CSV or TSV.

    The filename extension is used first:
        .csv -> comma
        .tsv -> tab

    For other extensions, csv.Sniffer is used to inspect the file.
    """

    ext = os.path.splitext(mapping_file)[1].lower()

    if ext == ".csv":
        return ","

    if ext == ".tsv":
        return "\t"

    # Fall back to delimiter detection for files with another extension.
    with open(mapping_file, "r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter

    except csv.Error:
        raise ValueError(
            f"Could not determine whether mapping file "
            f"{mapping_file!r} is CSV or TSV."
        )


def load_uniprot_to_gene_map(
    mapping_file,
    source_column="From",
    target_column="To",
):
    """
    Load a UniProt -> gene-symbol mapping from a CSV or TSV file.

    Parameters
    ----------
    mapping_file
        Path to the CSV or TSV mapping file.

    source_column
        Column containing UniProt identifiers.
        Default: From

    target_column
        Column containing gene symbols.
        Default: To

    Example
    -------
    From    To
    Q8BUN5  Smad3
    P50220  Nkx2-1
    """

    delimiter = detect_mapping_delimiter(mapping_file)

    id_to_gene = {}

    with open(
        mapping_file,
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(handle, delimiter=delimiter)

        if reader.fieldnames is None:
            raise ValueError(
                f"Mapping file has no header: {mapping_file}"
            )

        # Strip accidental whitespace around header names.
        reader.fieldnames = [
            field.strip() if field else field
            for field in reader.fieldnames
        ]

        missing = [
            column
            for column in (source_column, target_column)
            if column not in reader.fieldnames
        ]

        if missing:
            raise ValueError(
                f"Mapping file {mapping_file!r} is missing required "
                f"column(s): {', '.join(missing)}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )

        skipped = 0

        for row in reader:
            source_value = (row.get(source_column) or "").strip()
            target_value = (row.get(target_column) or "").strip()

            if not source_value or not target_value:
                skipped += 1
                continue

            # Normalize the source identifier here as well. This permits
            # mapping files containing either bare accessions or complete
            # UniProt identifiers.
            uniprot_id = extract_uniprot_id(source_value)

            id_to_gene[uniprot_id] = target_value

    delimiter_name = "TSV" if delimiter == "\t" else "CSV"

    print(
        f"Loaded {len(id_to_gene)} UniProt-to-gene mappings "
        f"from {mapping_file} ({delimiter_name}; "
        f"{source_column} -> {target_column})"
    )

    if skipped:
        print(
            f"Warning: skipped {skipped} mapping row(s) "
            f"with an empty {source_column!r} or {target_column!r} value."
        )

    return id_to_gene


def extract_uniprot_id(protein_id):
    """
    Normalize common representations of a UniProt identifier to
    the UniProt accession.

    Examples
    --------
    Q8BUN5
        -> Q8BUN5

    sp|Q8BUN5|SMAD3_MOUSE
        -> Q8BUN5

    tr|A0A123|ABC_MOUSE
        -> A0A123

    >sp|Q8BUN5|SMAD3_MOUSE
        -> Q8BUN5

    sp|Q8BUN5|SMAD3_MOUSE Mothers against decapentaplegic homolog 3
        -> Q8BUN5

    >sp|Q8BUN5|SMAD3_MOUSE Mothers against decapentaplegic homolog 3
        -> Q8BUN5
    """

    if protein_id is None:
        return ""

    protein_id = str(protein_id).strip()

    # Remove FASTA header marker if present.
    if protein_id.startswith(">"):
        protein_id = protein_id[1:].strip()

    # UniProt identifier is expected in the first whitespace-delimited
    # token even if a complete FASTA description is supplied.
    first_token = protein_id.split()[0] if protein_id else ""

    parts = first_token.split("|")

    # Canonical UniProt FASTA forms:
    # sp|ACCESSION|ENTRY_NAME
    # tr|ACCESSION|ENTRY_NAME
    if len(parts) >= 3 and parts[0].lower() in ("sp", "tr"):
        return parts[1].strip()

    # Handle an embedded sp|... or tr|... expression defensively.
    match = re.search(
        r"(?:^|[>\s])(?:sp|tr)\|([^|\s]+)\|",
        protein_id,
        flags=re.IGNORECASE,
    )

    if match:
        return match.group(1).strip()

    # Otherwise assume that the supplied value itself is the accession.
    return first_token


def map_to_gene_symbol(
    protein_id,
    id_to_gene,
    warn_missing=True,
):
    """
    Convert a protein identifier to a gene symbol.

    If no corresponding mapping exists, preserve the original protein
    identifier and optionally issue a warning.
    """

    if not id_to_gene:
        return protein_id

    uniprot_id = extract_uniprot_id(protein_id)

    gene_symbol = id_to_gene.get(uniprot_id)

    if gene_symbol:
        return gene_symbol

    if warn_missing:
        print(
            f"Warning: UniProt accession {uniprot_id!r} "
            f"extracted from {protein_id!r} was not found "
            f"in the mapping file."
        )

    # Preserve the original identifier so a failed mapping is obvious.
    return protein_id


def map_and_deduplicate(protein_ids, id_to_gene):
    seen = set()
    mapped = []

    for protein_id in protein_ids:
        gene_symbol = map_to_gene_symbol(
            protein_id,
            id_to_gene,
        )

        if gene_symbol in seen:
            continue

        mapped.append(gene_symbol)
        seen.add(gene_symbol)

    return mapped


# ----------------------------------------------------------------------
# SPPIDER-seq prediction parsing
# ----------------------------------------------------------------------

def parse_prediction_file(source_name, handle):
    first_line = handle.readline()

    metadata = (
        parse_header(first_line)
        or parse_filename(source_name)
    )

    if metadata is None:
        return None

    positions = []
    amino_acids = []
    probabilities = []

    for line in handle:
        line = line.strip()

        if (
            not line
            or line.startswith("#")
            or line.startswith("Position")
        ):
            continue

        parts = line.split("\t")

        if len(parts) < 3:
            continue

        try:
            position = int(parts[0])
            probability = float(parts[2])

        except ValueError:
            continue

        positions.append(position)
        amino_acids.append(parts[1])
        probabilities.append(probability)

    if not positions:
        return None

    return Prediction(
        source_name=source_name,
        query=metadata["query"],
        partner=metadata["partner"],
        model=metadata["model"],
        positions=positions,
        amino_acids=amino_acids,
        probabilities=probabilities,
    )


def iter_zip_prediction_files(zip_path, source_prefix=None):
    with zipfile.ZipFile(zip_path) as archive:

        for info in archive.infolist():

            if info.is_dir() or not info.filename.endswith(".txt"):
                continue

            source_name = info.filename

            if source_prefix:
                source_name = f"{source_prefix}/{info.filename}"

            with archive.open(info) as raw:

                text_handle = TextIOWrapper(
                    raw,
                    encoding="utf-8",
                    newline="",
                )

                yield source_name, text_handle


def iter_prediction_files(input_path):
    if os.path.isdir(input_path):

        for root, _, files in os.walk(input_path):

            for filename in files:
                path = os.path.join(root, filename)

                if filename.endswith(".txt"):

                    with open(path, newline="") as handle:
                        yield path, handle

                elif zipfile.is_zipfile(path):

                    yield from iter_zip_prediction_files(
                        path,
                        path,
                    )

        return

    if zipfile.is_zipfile(input_path):

        yield from iter_zip_prediction_files(
            input_path,
        )

        return

    raise ValueError(
        f"Input must be a directory or zip file: {input_path}"
    )


def read_predictions(input_path):
    predictions = []

    for source_name, handle in iter_prediction_files(input_path):

        with handle:
            prediction = parse_prediction_file(
                source_name,
                handle,
            )

        if prediction is not None:
            predictions.append(prediction)

    return predictions


def validate_prediction_group(query, model, predictions):
    expected_positions = predictions[0].positions
    expected_amino_acids = predictions[0].amino_acids

    for prediction in predictions[1:]:

        if prediction.positions != expected_positions:

            raise ValueError(
                f"Position mismatch for query={query!r}, "
                f"model={model!r}: "
                f"{prediction.source_name}"
            )

        if prediction.amino_acids != expected_amino_acids:

            raise ValueError(
                f"Amino-acid mismatch for query={query!r}, "
                f"model={model!r}: "
                f"{prediction.source_name}"
            )


def percentile(sorted_values, percent):
    if not sorted_values:
        return 0.0

    if len(sorted_values) == 1:
        return sorted_values[0]

    rank = (len(sorted_values) - 1) * percent

    lower = int(rank)
    upper = min(
        lower + 1,
        len(sorted_values) - 1,
    )

    fraction = rank - lower

    return (
        sorted_values[lower] * (1.0 - fraction)
        + sorted_values[upper] * fraction
    )


def aggregate_positions(predictions, id_to_gene=None):
    validate_prediction_group(
        predictions[0].query,
        predictions[0].model,
        predictions,
    )

    rows = []
    distributions = []

    for index, (position, amino_acid) in enumerate(
        zip(
            predictions[0].positions,
            predictions[0].amino_acids,
        )
    ):
        probabilities_by_partner = [
            (
                prediction.partner,
                prediction.probabilities[index],
            )
            for prediction in predictions
        ]

        row = {
            "position": position,
            "amino_acid": amino_acid,
        }

        for cutoff in CUTOFFS:

            partners = map_and_deduplicate(
                sorted(
                    partner
                    for partner, probability
                    in probabilities_by_partner
                    if probability >= cutoff
                ),
                id_to_gene,
            )

            row[f"pPPI_ge_{cutoff:g}"] = ";".join(partners)

        probabilities = sorted(
            probability
            for _, probability
            in probabilities_by_partner
        )

        distributions.append(
            {
                "position": position,
                "mean": sum(probabilities) / len(probabilities),
                "q25": percentile(probabilities, 0.25),
                "q75": percentile(probabilities, 0.75),
                "max": probabilities[-1],
            }
        )

        rows.append(row)

    return rows, distributions


# ----------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------

def write_csv(path, rows):
    fieldnames = [
        "position",
        "amino_acid",
        "pPPI_ge_0.7",
        "pPPI_ge_0.5",
        "pPPI_ge_0.3",
        "pPPI_ge_0.1",
    ]

    with open(path, "w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(rows)


def write_plot(path, query, model, distributions):
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    positions = [
        item["position"]
        for item in distributions
    ]

    means = [
        item["mean"]
        for item in distributions
    ]

    q25 = [
        item["q25"]
        for item in distributions
    ]

    q75 = [
        item["q75"]
        for item in distributions
    ]

    max_values = [
        item["max"]
        for item in distributions
    ]

    fig, ax = plt.subplots(
        figsize=(12, 4),
        dpi=300,
    )

    ax.fill_between(
        positions,
        q25,
        q75,
        color="#bdbdbd",
        alpha=0.45,
        label="Interquartile range",
    )

    ax.plot(
        positions,
        means,
        color="#525252",
        linewidth=1.8,
        label="Mean",
    )

    ax.plot(
        positions,
        max_values,
        color="#de2d26",
        linewidth=1.1,
        label="Max",
    )

    ax.set_xlabel(
        f"Residue position (L={len(positions)})"
    )

    ax.set_ylabel(
        "PPI site probability"
    )

    ax.set_ylim(0, 1.05)

    ax.set_xlim(
        min(positions),
        max(positions),
    )

    ax.set_title(
        f"Query: {query}\nModel: {model}"
    )

    ax.grid(
        True,
        linestyle=":",
        linewidth=0.7,
    )

    ax.legend(
        loc="upper right"
    )

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)



def format_track_name(prediction, id_to_gene=None):
    """Build a track identifier from query, partner, and model."""

    query = map_to_gene_symbol(
        prediction.query,
        id_to_gene,
        warn_missing=False,
    )
    partner = map_to_gene_symbol(
        prediction.partner,
        id_to_gene,
        warn_missing=False,
    )

    return f"{query}||{partner}||{prediction.model}"


def write_tracks(
    path,
    predictions,
    id_to_gene=None,
    binary=False,
    threshold=0.5,
):
    """
    Write predictions in three-line track format:

        >query||partner||model
        AMINOACIDSEQUENCE
        probability vector

    With binary=True, the third line is a 0/1 PPI vector using
    probability >= threshold.
    """

    with open(path, "w", encoding="utf-8") as handle:

        for prediction in sorted(
            predictions,
            key=lambda item: (
                item.query,
                item.partner,
                item.model,
                item.source_name,
            ),
        ):
            track_name = format_track_name(
                prediction,
                id_to_gene,
            )

            sequence = "".join(
                prediction.amino_acids
            )

            if binary:
                values = [
                    "1" if probability >= threshold else "0"
                    for probability in prediction.probabilities
                ]
            else:
                values = [
                    f"{probability:.3f}"
                    for probability in prediction.probabilities
                ]

            handle.write(f">{track_name}\n")
            handle.write(f"{sequence}\n")
            if binary:
                handle.write("".join(values) + "\n")
            else:
                handle.write(" ".join(values) + "\n")


def write_outputs(
    predictions,
    output_dir,
    id_to_gene=None,
    write_png=True,
    write_csv_output=False,
    write_track_output=False,
    binary_tracks=False,
    binary_threshold=0.5,
):
    os.makedirs(
        output_dir,
        exist_ok=True,
    )

    grouped = defaultdict(list)

    for prediction in predictions:
        grouped[
            (
                prediction.query,
                prediction.model,
            )
        ].append(prediction)

    written = []

    for (query, model), group in sorted(grouped.items()):

        rows, distributions = aggregate_positions(
            group,
            id_to_gene,
        )

        display_query = map_to_gene_symbol(
            query,
            id_to_gene,
        )

        base = (
            f"{sanitize_filename(display_query)}"
            f"__{sanitize_filename(model)}"
        )

        if write_png:
            png_path = os.path.join(
                output_dir,
                f"{base}.png",
            )

            write_plot(
                png_path,
                display_query,
                model,
                distributions,
            )

            written.append(png_path)

        if write_csv_output:
            csv_path = os.path.join(
                output_dir,
                f"{base}.csv",
            )

            write_csv(
                csv_path,
                rows,
            )

            written.append(csv_path)

        if write_track_output:
            suffix = (
                "binary_tracks"
                if binary_tracks
                else "tracks"
            )

            track_path = os.path.join(
                output_dir,
                f"{base}.{suffix}.txt",
            )

            write_tracks(
                track_path,
                group,
                id_to_gene=id_to_gene,
                binary=binary_tracks,
                threshold=binary_threshold,
            )

            written.append(track_path)

    return written

def print_verbose_predictions(predictions, id_to_gene=None):
    """
    Print detailed information for every parsed prediction to stdout.
    """

    print("\nParsed predictions:")
    print("Query\tPartner\tModel\tSource")

    for prediction in sorted(
        predictions,
        key=lambda item: (
            item.query,
            item.partner,
            item.model,
            item.source_name,
        ),
    ):
        query = map_to_gene_symbol(
            prediction.query,
            id_to_gene,
            warn_missing=False,
        )

        partner = map_to_gene_symbol(
            prediction.partner,
            id_to_gene,
            warn_missing=False,
        )

        print(
            f"{query}\t"
            f"{partner}\t"
            f"{prediction.model}\t"
            f"{prediction.source_name}"
        )


# ----------------------------------------------------------------------
# Command-line interface
# ----------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Parse SPPIDER-seq per-pair prediction files and generate "
            "per-query/per-model plots, partner cutoff CSV files, and/or "
            "three-line sequence tracks. PNG plots are generated by default "
            "when no output-type option is specified."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=(
            "Directory or ZIP archive containing SPPIDER-seq "
            "per-pair .txt prediction files."
        ),
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default="sppider_parsed_outputs",
        help=(
            "Directory for generated output files. "
            "Default: sppider_parsed_outputs"
        ),
    )

    output_group = parser.add_argument_group("output types")

    output_group.add_argument(
        "-p",
        "--png",
        action="store_true",
        help=(
            "Generate PNG probability-distribution plots. "
            "This is the default when none of --png, --csv, or --tracks "
            "is specified."
        ),
    )

    output_group.add_argument(
        "-c",
        "--csv",
        action="store_true",
        help=(
            "Generate per-query/per-model CSV files listing partners "
            "meeting the probability cutoffs."
        ),
    )

    output_group.add_argument(
        "-t",
        "--tracks",
        action="store_true",
        help=(
            "Generate three-line track files for individual query-partner "
            "predictions: header, amino-acid sequence, and probability vector."
        ),
    )

    output_group.add_argument(
        "-b",
        "--binary",
        action="store_true",
        help=(
            "Write track vectors as binary PPI calls instead of probabilities. "
            "Requires --tracks. Values >= --threshold become 1; others become 0."
        ),
    )

    output_group.add_argument(
        "-T",
        "--threshold",
        type=float,
        default=0.5,
        help=(
            "Probability threshold used with --binary. "
            "Default: 0.5"
        ),
    )

    mapping_group = parser.add_argument_group("UniProt-to-gene mapping")

    mapping_group.add_argument(
        "-g",
        "--convert-uniprot-to-gene",
        action="store_true",
        help=(
            "Convert UniProt identifiers in query and partner names "
            "to gene symbols using --mapping-file. "
            "Default: off."
        ),
    )

    mapping_group.add_argument(
        "-m",
        "--mapping-file",
        default=None,
        help=(
            "CSV or TSV file used for UniProt-to-gene conversion. "
            "Required with --convert-uniprot-to-gene."
        ),
    )

    mapping_group.add_argument(
        "-C",
        "--mapping-columns",
        nargs=2,
        metavar=("SOURCE_COLUMN", "TARGET_COLUMN"),
        default=["From", "To"],
        help=(
            "Mapping-file columns: first the UniProt identifier column, "
            "then the gene-symbol column. Default: From To"
        ),
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help=(
            "Print detailed information for every parsed prediction, "
            "including query, partner, model, and source file."
        ),
    )

    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()

    # PNG is the default only when no output type is selected explicitly.
    if not any((args.png, args.csv, args.tracks)):
        args.png = True

    if args.binary and not args.tracks:
        parser.error(
            "--binary requires --tracks."
        )

    if not 0.0 <= args.threshold <= 1.0:
        parser.error(
            "--threshold must be between 0.0 and 1.0."
        )

    if (
        args.convert_uniprot_to_gene
        and not args.mapping_file
    ):
        parser.error(
            "--mapping-file is required when "
            "--convert-uniprot-to-gene is specified."
        )

    if args.convert_uniprot_to_gene:

        source_column, target_column = args.mapping_columns

        id_to_gene = load_uniprot_to_gene_map(
            args.mapping_file,
            source_column=source_column,
            target_column=target_column,
        )

        if not id_to_gene:
            raise SystemExit(
                "The mapping file was loaded but contained "
                "no usable UniProt-to-gene mappings."
            )

    else:
        id_to_gene = None

    predictions = read_predictions(
        args.input
    )

    if not predictions:
        raise SystemExit(
            "No SPPIDER-seq prediction text files were found."
        )

    queries = sorted(
        {
            prediction.query
            for prediction in predictions
        }
    )

    partners = sorted(
        {
            prediction.partner
            for prediction in predictions
        }
    )

    written = write_outputs(
        predictions,
        args.output_dir,
        id_to_gene,
        write_png=args.png,
        write_csv_output=args.csv,
        write_track_output=args.tracks,
        binary_tracks=args.binary,
        binary_threshold=args.threshold,
    )

    models = sorted(
        {
            prediction.model
            for prediction in predictions
        }
    )

    print("\nSummary")
    print("-------")
    print(f"Prediction files parsed: {len(predictions)}")
    print(f"Unique queries:          {len(queries)}")
    print(f"Unique partners:         {len(partners)}")
    print(f"Models found:            {', '.join(models)}")

    output_types = []
    if args.png:
        output_types.append("PNG")
    if args.csv:
        output_types.append("CSV")
    if args.tracks:
        if args.binary:
            output_types.append(
                f"binary tracks (threshold={args.threshold:g})"
            )
        else:
            output_types.append("probability tracks")

    print(f"Output types:            {', '.join(output_types)}")
    print(f"Output files written:    {len(written)}")
    print(f"Output directory:        {args.output_dir}")

    print("\nSaved files:")
    for path in written:
        print(f"  {path}")

    if args.verbose:
        print_verbose_predictions(
            predictions,
            id_to_gene,
        )


if __name__ == "__main__":
    main()
