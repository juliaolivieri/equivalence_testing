import argparse
import io
from pathlib import Path
from urllib.request import urlopen

import pandas as pd


def make_unique(values):
    seen = {}
    out = []
    for value in values:
        key = str(value)
        if key not in seen:
            seen[key] = 0
            out.append(key)
        else:
            seen[key] += 1
            out.append(f"{key}__dup{seen[key]}")
    return out


def get_args():
    parser = argparse.ArgumentParser(
        description="Prepare a GEO count matrix and metadata file for batch equivalence analysis."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Local path or URL to the GEO count matrix.",
    )
    parser.add_argument(
        "--sample-sheet",
        required=True,
        help="CSV with columns sample,condition describing the comparison to run.",
    )
    parser.add_argument(
        "--output-counts",
        required=True,
        help="Output CSV path for the counts-only matrix.",
    )
    parser.add_argument(
        "--output-meta",
        required=True,
        help="Output CSV path for the metadata table.",
    )
    parser.add_argument(
        "--gene-column",
        default="gene_id",
        help="Column in the input matrix to use as the gene identifier.",
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Field separator for the input matrix. Default is tab.",
    )
    parser.add_argument(
        "--deduplicate-gene-column",
        action="store_true",
        help="If the gene identifier column contains duplicates, make them unique by appending suffixes.",
    )
    return parser.parse_args()


def load_table(path_or_url, sep):
    if path_or_url.startswith(("http://", "https://", "ftp://")):
        with urlopen(path_or_url) as handle:
            return pd.read_csv(handle, sep=sep, compression="infer")
    return pd.read_csv(path_or_url, sep=sep, compression="infer")


def main():
    args = get_args()
    sample_df = pd.read_csv(args.sample_sheet)
    if list(sample_df.columns) != ["sample", "condition"]:
        raise ValueError("sample sheet must have exactly two columns: sample,condition")

    counts_df = load_table(args.input, args.sep)
    needed = [args.gene_column] + sample_df["sample"].tolist()
    missing = [col for col in needed if col not in counts_df.columns]
    if missing:
        raise ValueError("missing columns in input matrix: {}".format(", ".join(missing)))

    out_counts = counts_df[needed].copy()
    if args.deduplicate_gene_column:
        gene_values = out_counts[args.gene_column].astype(str)
        if gene_values.duplicated().any():
            out_counts[args.gene_column] = make_unique(gene_values.tolist())

    out_counts_path = Path(args.output_counts)
    out_meta_path = Path(args.output_meta)
    out_counts_path.parent.mkdir(parents=True, exist_ok=True)
    out_meta_path.parent.mkdir(parents=True, exist_ok=True)
    out_counts.to_csv(out_counts_path, index=False)
    sample_df.to_csv(out_meta_path, index=False)

    print("Wrote counts to {}".format(out_counts_path))
    print("Wrote metadata to {}".format(out_meta_path))


if __name__ == "__main__":
    main()
