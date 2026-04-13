import argparse
from pathlib import Path

import pandas as pd


DEFAULT_COLUMNS = [
    "analysis_id",
    "gse_accession",
    "selected_for_analysis",
    "analysis_ready",
    "local_counts_path",
    "local_meta_path",
    "condition_column",
    "control_label",
    "treatment_label",
    "comparison_name",
    "description",
    "species",
    "source_geo_url",
    "source_screening_status",
    "source_matrix_type_assessment",
    "source_metadata_sufficiency",
    "source_confidence",
    "comparison_description",
    "replicates_by_condition",
    "delta",
    "low_count_threshold",
    "min_nonzero_per_group",
    "notes",
]


def get_args():
    parser = argparse.ArgumentParser(
        description="Create an analysis manifest from GEO screening results."
    )
    parser.add_argument(
        "--screening-csv",
        default="download_geo_datasets/geo_bulk_rnaseq_screening.csv",
        help="Path to the GEO screening CSV.",
    )
    parser.add_argument(
        "--output",
        default="download_geo_datasets/analysis_manifest.csv",
        help="Where to write the analysis manifest CSV.",
    )
    parser.add_argument(
        "--include-status",
        nargs="+",
        default=["eligible"],
        help="Screening statuses to include in the manifest.",
    )
    parser.add_argument(
        "--append-smoke-test",
        action="store_true",
        help="Append a local smoke-test row using the intuitive test dataset.",
    )
    return parser.parse_args()


def build_manifest(screening_df):
    rows = []
    for _, row in screening_df.iterrows():
        rows.append(
            {
                "analysis_id": row["gse_accession"],
                "gse_accession": row["gse_accession"],
                "selected_for_analysis": False,
                "analysis_ready": False,
                "local_counts_path": "",
                "local_meta_path": "",
                "condition_column": "",
                "control_label": "",
                "treatment_label": "",
                "comparison_name": row.get("comparison_description", ""),
                "description": row.get("title", ""),
                "species": row.get("organism", ""),
                "source_geo_url": row.get("geo_url", ""),
                "source_screening_status": row.get("inclusion_status", ""),
                "source_matrix_type_assessment": row.get("matrix_type_assessment", ""),
                "source_metadata_sufficiency": row.get("metadata_sufficiency", ""),
                "source_confidence": row.get("confidence", ""),
                "comparison_description": row.get("comparison_description", ""),
                "replicates_by_condition": row.get("replicates_by_condition", ""),
                "delta": 1,
                "low_count_threshold": 5,
                "min_nonzero_per_group": 2,
                "notes": "",
            }
        )
    manifest = pd.DataFrame(rows, columns=DEFAULT_COLUMNS)
    return manifest


def smoke_test_row():
    return {
        "analysis_id": "TEST_INTUITIVE",
        "gse_accession": "TEST_INTUITIVE",
        "selected_for_analysis": True,
        "analysis_ready": True,
        "local_counts_path": "test_data/intuitive_data.csv",
        "local_meta_path": "test_data/intuitive_meta.csv",
        "condition_column": "condition",
        "control_label": "ctrl",
        "treatment_label": "trt",
        "comparison_name": "ctrl_vs_trt",
        "description": "Local smoke test using the intuitive toy dataset",
        "species": "Mus musculus",
        "source_geo_url": "",
        "source_screening_status": "local_test",
        "source_matrix_type_assessment": "raw_counts_confirmed",
        "source_metadata_sufficiency": "clear",
        "source_confidence": "high",
        "comparison_description": "Balanced toy control-versus-treatment comparison for pipeline validation",
        "replicates_by_condition": "ctrl=4; trt=4",
        "delta": 1,
        "low_count_threshold": 5,
        "min_nonzero_per_group": 2,
        "notes": "Useful for validating the manifest-driven runner before GEO downloads are prepared.",
    }


def main():
    args = get_args()
    screening_path = Path(args.screening_csv)
    output_path = Path(args.output)

    screening_df = pd.read_csv(screening_path)
    screening_df = screening_df[
        screening_df["inclusion_status"].isin(args.include_status)
    ].copy()
    manifest = build_manifest(screening_df)

    if args.append_smoke_test:
        manifest = pd.concat(
            [pd.DataFrame([smoke_test_row()], columns=DEFAULT_COLUMNS), manifest],
            ignore_index=True,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_path, index=False)
    print(
        "Wrote {} rows to {}".format(
            manifest.shape[0],
            output_path,
        )
    )


if __name__ == "__main__":
    main()
