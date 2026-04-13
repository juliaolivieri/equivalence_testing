import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(
        description="Run the DESeq2 normalization + equivalence pipeline from an analysis manifest."
    )
    parser.add_argument(
        "--manifest",
        default="download_geo_datasets/analysis_manifest.csv",
        help="Path to the analysis manifest CSV.",
    )
    parser.add_argument(
        "--runs-dir",
        default="runs/batch_equiv",
        help="Base output directory for batch runs.",
    )
    parser.add_argument(
        "--analysis-ids",
        nargs="+",
        default=None,
        help="Optional subset of analysis_id values to run.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun analyses even if the result file already exists.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the analyses that would run without executing them.",
    )
    parser.add_argument(
        "--python-executable",
        default=sys.executable,
        help="Python interpreter to use for the equivalence script.",
    )
    parser.add_argument(
        "--rscript-executable",
        default="Rscript",
        help="Rscript executable to use for DESeq2 normalization.",
    )
    return parser.parse_args()


def run_command(cmd, cwd):
    subprocess.run(cmd, cwd=cwd, check=True)


def parse_bool(value):
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def format_delta(value):
    return str(float(value))


def summarize_result(result_path, counts_path, meta_path, condition_column, row):
    result_df = pd.read_csv(result_path, index_col=0)
    counts_df = pd.read_csv(counts_path, index_col=0)
    meta_df = pd.read_csv(meta_path, index_col=0)

    meta_df = meta_df.loc[counts_df.columns]
    group_sizes = meta_df[condition_column].value_counts().to_dict()
    library_sizes = counts_df.sum(axis=0)

    summary = {
        "analysis_id": row["analysis_id"],
        "gse_accession": row["gse_accession"],
        "description": row.get("description", ""),
        "species": row.get("species", ""),
        "comparison_name": row.get("comparison_name", ""),
        "comparison_description": row.get("comparison_description", ""),
        "condition_column": condition_column,
        "n_group1": min(group_sizes.values()) if group_sizes else 0,
        "n_group2": max(group_sizes.values()) if len(group_sizes) > 1 else 0,
        "total_n": int(sum(group_sizes.values())),
        "num_genes_tested": int(result_df.shape[0]),
        "num_different": int((result_df["category"] == "different").sum()),
        "num_equivalent": int((result_df["category"] == "equivalent").sum()),
        "num_inconclusive": int((result_df["category"] == "inconclusive").sum()),
        "num_low_information": int((result_df["category"] == "low_information").sum()),
        "frac_different": float((result_df["category"] == "different").mean()),
        "frac_equivalent": float((result_df["category"] == "equivalent").mean()),
        "frac_inconclusive": float((result_df["category"] == "inconclusive").mean()),
        "frac_low_information": float((result_df["category"] == "low_information").mean()),
        "mean_library_size": float(library_sizes.mean()),
        "median_library_size": float(library_sizes.median()),
        "delta": row["delta"],
        "low_count_threshold": row["low_count_threshold"],
        "min_nonzero_per_group": row["min_nonzero_per_group"],
        "status": "completed",
        "error_message": "",
        "result_path": str(result_path),
    }
    return summary


def main():
    args = get_args()
    repo_root = Path(__file__).resolve().parent.parent
    manifest_path = repo_root / args.manifest
    runs_dir = repo_root / args.runs_dir
    runs_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(manifest_path)
    manifest["selected_for_analysis"] = manifest["selected_for_analysis"].apply(parse_bool)
    manifest["analysis_ready"] = manifest["analysis_ready"].apply(parse_bool)

    queued = manifest[
        manifest["selected_for_analysis"] & manifest["analysis_ready"]
    ].copy()
    if args.analysis_ids:
        queued = queued[queued["analysis_id"].isin(args.analysis_ids)].copy()

    summary_rows = []
    for _, row in queued.iterrows():
        analysis_id = row["analysis_id"]
        analysis_dir = runs_dir / analysis_id
        analysis_dir.mkdir(parents=True, exist_ok=True)

        counts_path = repo_root / row["local_counts_path"]
        meta_path = repo_root / row["local_meta_path"]
        delta_str = format_delta(row["delta"])
        result_prefix = analysis_dir / f"{analysis_id}_{row['condition_column']}_{delta_str}"
        result_path = Path(str(result_prefix) + "_results.csv")

        if result_path.exists() and not args.force:
            summary_rows.append(
                summarize_result(
                    result_path,
                    counts_path,
                    meta_path,
                    row["condition_column"],
                    row,
                )
            )
            continue

        if args.dry_run:
            print(f"[dry-run] {analysis_id}: {counts_path} vs {meta_path}")
            continue

        try:
            normalized_prefix = analysis_dir / analysis_id
            run_command(
                [
                    args.rscript_executable,
                    "scripts/deseq_normalization.R",
                    str(counts_path),
                    str(meta_path),
                    row["condition_column"],
                    str(normalized_prefix),
                ],
                cwd=repo_root,
            )
            run_command(
                [
                    args.python_executable,
                    "scripts/equiv_test_vectorized.py",
                    "--savename",
                    analysis_id,
                    "--outpath",
                    str(analysis_dir) + "/",
                    "--infile",
                    str(normalized_prefix) + "_deseq2_normalized_counts.csv",
                    "--meta",
                    str(meta_path),
                    "--condition",
                    row["condition_column"],
                    "--delta",
                    delta_str,
                    "--low_count_threshold",
                    str(row["low_count_threshold"]),
                    "--min_nonzero_per_group",
                    str(int(row["min_nonzero_per_group"])),
                    "--deseq2_results",
                    str(normalized_prefix) + f"_{row['condition_column']}_deseq2_results.csv",
                ],
                cwd=repo_root,
            )
            summary_rows.append(
                summarize_result(
                    result_path,
                    counts_path,
                    meta_path,
                    row["condition_column"],
                    row,
                )
            )
        except Exception as exc:
            summary_rows.append(
                {
                    "analysis_id": analysis_id,
                    "gse_accession": row["gse_accession"],
                    "description": row.get("description", ""),
                    "species": row.get("species", ""),
                    "comparison_name": row.get("comparison_name", ""),
                    "comparison_description": row.get("comparison_description", ""),
                    "condition_column": row["condition_column"],
                    "n_group1": None,
                    "n_group2": None,
                    "total_n": None,
                    "num_genes_tested": None,
                    "num_different": None,
                    "num_equivalent": None,
                    "num_inconclusive": None,
                    "num_low_information": None,
                    "frac_different": None,
                    "frac_equivalent": None,
                    "frac_inconclusive": None,
                    "frac_low_information": None,
                    "mean_library_size": None,
                    "median_library_size": None,
                    "delta": row["delta"],
                    "low_count_threshold": row["low_count_threshold"],
                    "min_nonzero_per_group": row["min_nonzero_per_group"],
                    "status": "failed",
                    "error_message": str(exc),
                    "result_path": "",
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    summary_path = runs_dir / "batch_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Wrote batch summary to {summary_path}")


if __name__ == "__main__":
    main()
