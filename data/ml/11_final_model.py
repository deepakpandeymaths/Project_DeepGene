from pathlib import Path
import pandas as pd


# ============================================================
# PROJECT DEEPGENE
# STEP 11 — FINAL MODEL READINESS
# ============================================================

PROJECT_ROOT = Path(r"D:\Project_DeepGene")

FEATURE_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_ml_features_v2.csv"
)

MANIFEST_FILE = (
    PROJECT_ROOT
    / "data"
    / "ml_results"
    / "10_final_feature_manifest_v2.csv"
)

OUTPUT_FEATURE_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "deepgene_final_model_features_v1.csv"
)

REPORT_FILE = (
    PROJECT_ROOT
    / "data"
    / "ml_results"
    / "11_final_model.txt"
)


# Features removed during Step 10
REMOVED_CONSTANT_FEATURES = {
    "phenotype_available",
    "genomic_coordinates_present",
    "hgvs_present",
}

REMOVED_REDUNDANT_FEATURE = "conflict_flag"


def main():

    print("=" * 70)
    print("PROJECT DEEPGENE")
    print("STEP 11 — FINAL MODEL READINESS")
    print("=" * 70)

    # --------------------------------------------------------
    # 1. Load Step 10 feature dataset
    # --------------------------------------------------------
    if not FEATURE_FILE.exists():
        raise FileNotFoundError(
            f"Final feature dataset not found:\n{FEATURE_FILE}"
        )

    df = pd.read_csv(FEATURE_FILE)

    if df.empty:
        raise ValueError("Final feature dataset is empty.")

    print(f"Input rows: {len(df)}")
    print(f"Input columns: {len(df.columns)}")

    # --------------------------------------------------------
    # 2. Load Step 10 manifest
    # --------------------------------------------------------
    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(
            f"Feature manifest not found:\n{MANIFEST_FILE}"
        )

    manifest = pd.read_csv(MANIFEST_FILE)

    if "feature" not in manifest.columns:
        raise ValueError(
            "Feature manifest must contain a 'feature' column."
        )

    manifest_features = (
        manifest["feature"]
        .dropna()
        .astype(str)
        .tolist()
    )

    # --------------------------------------------------------
    # 3. Reconstruct the final retained feature set
    #
    # Step 10 established:
    #   18 original candidate features
    #   - 3 constant features
    #   - 1 redundant feature
    #   = 14 final features
    #
    # The manifest may contain the original candidate list,
    # so we explicitly apply the Step 10 removal decisions.
    # --------------------------------------------------------
    final_features = [
        feature
        for feature in manifest_features
        if feature not in REMOVED_CONSTANT_FEATURES
        and feature != REMOVED_REDUNDANT_FEATURE
    ]

    # Remove accidental duplicates while preserving order.
    final_features = list(dict.fromkeys(final_features))

    if len(final_features) != 14:
        raise ValueError(
            "Step 10 reconstruction did not produce the expected "
            f"14 final features. Found {len(final_features)}:\n"
            f"{final_features}"
        )

    # --------------------------------------------------------
    # 4. Validate feature availability
    # --------------------------------------------------------
    missing_features = [
        feature
        for feature in final_features
        if feature not in df.columns
    ]

    if missing_features:
        raise ValueError(
            "Final features missing from Step 10 dataset:\n"
            f"{missing_features}"
        )

    # --------------------------------------------------------
    # 5. Validate final feature matrix
    # --------------------------------------------------------
    X = df[final_features].copy()

    missing_values = int(X.isna().sum().sum())

    if missing_values != 0:
        raise ValueError(
            f"Final feature matrix contains {missing_values} missing values."
        )

    non_numeric_features = [
        feature
        for feature in final_features
        if not pd.api.types.is_numeric_dtype(X[feature])
    ]

    if non_numeric_features:
        raise ValueError(
            "Non-numeric final features detected:\n"
            f"{non_numeric_features}"
        )

    finite_matrix = bool(
        X.apply(
            lambda column: pd.Series(
                pd.api.types.is_number(value) and pd.notna(value)
                for value in column
            ),
            axis=0,
        ).all().all()
    )

    if not finite_matrix:
        raise ValueError(
            "Final feature matrix contains invalid numeric values."
        )

    # --------------------------------------------------------
    # 6. Validate variant IDs
    # --------------------------------------------------------
    if "variant_id" not in df.columns:
        raise ValueError(
            "variant_id is required for traceability."
        )

    variant_ids = df["variant_id"]

    if variant_ids.isna().any():
        raise ValueError("variant_id contains missing values.")

    if variant_ids.duplicated().any():
        raise ValueError("variant_id contains duplicates.")

    # --------------------------------------------------------
    # 7. Write final model feature dataset
    # --------------------------------------------------------
    output_df = X.copy()
    output_df.insert(0, "variant_id", variant_ids)

    output_df.to_csv(
        OUTPUT_FEATURE_FILE,
        index=False
    )

    # --------------------------------------------------------
    # 8. Independent target status
    # --------------------------------------------------------
    # No independent ground-truth target currently exists.
    target_available = False

    # --------------------------------------------------------
    # 9. Write report
    # --------------------------------------------------------
    report_lines = [
        "PROJECT DEEPGENE",
        "STEP 11 — FINAL MODEL READINESS",
        "",
        f"Input dataset: {FEATURE_FILE}",
        f"Variants: {len(df)}",
        f"Input columns: {len(df.columns)}",
        "",
        f"Manifest feature entries: {len(manifest_features)}",
        f"Final retained features: {len(final_features)}",
        "",
        "REMOVED CONSTANT FEATURES:",
    ]

    report_lines.extend(
        [
            f"  - {feature}"
            for feature in sorted(REMOVED_CONSTANT_FEATURES)
        ]
    )

    report_lines.extend(
        [
            "",
            "REMOVED REDUNDANT FEATURE:",
            f"  - {REMOVED_REDUNDANT_FEATURE}",
            "",
            "FINAL FEATURES:",
        ]
    )

    report_lines.extend(
        [
            f"  - {feature}"
            for feature in final_features
        ]
    )

    report_lines.extend(
        [
            "",
            f"Missing feature values: {missing_values}",
            f"Non-numeric features: {len(non_numeric_features)}",
            f"Final feature matrix finite: {finite_matrix}",
            "",
            "INDEPENDENT TARGET:",
            "NOT AVAILABLE",
            "",
            "ClinVar clinical_significance:",
            "NOT used as an independent ML target.",
            "",
            "FINAL MODEL STATUS:",
            "READY FOR SUPERVISED MODELING ONLY AFTER",
            "AN INDEPENDENT GROUND-TRUTH TARGET IS PROVIDED.",
            "",
            "Supervised model fitting:",
            "NOT PERFORMED.",
            "",
            "Model performance metrics:",
            "NOT AVAILABLE.",
            "",
            f"Final model feature dataset:",
            str(OUTPUT_FEATURE_FILE),
            "",
            f"Report:",
            str(REPORT_FILE),
            "",
            "STATUS: PASS",
        ]
    )

    REPORT_FILE.write_text(
        "\n".join(report_lines),
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # 10. Terminal summary
    # --------------------------------------------------------
    print(f"Manifest feature entries: {len(manifest_features)}")
    print(f"Final retained features: {len(final_features)}")

    print()
    print("Removed constant features:")
    for feature in sorted(REMOVED_CONSTANT_FEATURES):
        print(f"  - {feature}")

    print()
    print(f"Removed redundant feature: {REMOVED_REDUNDANT_FEATURE}")

    print()
    print("Final features:")
    for feature in final_features:
        print(f"  - {feature}")

    print()
    print(f"Missing feature values: {missing_values}")
    print(f"Final feature matrix finite: {finite_matrix}")

    print()
    print("Independent target: NOT AVAILABLE")
    print("Supervised model fitting: NOT PERFORMED")

    print()
    print("Final model feature dataset:")
    print(OUTPUT_FEATURE_FILE)

    print()
    print("Report:")
    print(REPORT_FILE)

    print()
    print("STATUS: PASS")
    print("=" * 70)


if __name__ == "__main__":
    main()