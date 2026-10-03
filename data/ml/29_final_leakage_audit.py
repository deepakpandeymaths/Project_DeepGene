from pathlib import Path

import pandas as pd


# ============================================================
# DEEPGENE ML STEP 29
# FINAL LEAKAGE & TARGET INDEPENDENCE AUDIT
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
REPORT_FILE = Path("data/analysis_results/29_final_leakage_audit.txt")

TARGET = "functional_target"

ALLOWED_FEATURES = [
    "protein_position",
    "same_amino_acid",
    "reference_hydrophobicity",
    "alternate_hydrophobicity",
    "delta_hydrophobicity",
    "absolute_delta_hydrophobicity",
    "reference_charge",
    "alternate_charge",
    "delta_charge",
    "absolute_delta_charge",
    "reference_volume",
    "alternate_volume",
    "delta_volume",
    "absolute_delta_volume",
]


# Columns intentionally present as dataset metadata but NOT
# used as model predictors.
ALLOWED_METADATA = [
    "variant_id",
    "functional_protein",
    "protein_variant",
    "reference_aa",
    "alternate_aa",
    "amino_acid_change_type",
    "functional_effect",
    "phenotype",
    "inheritance",
    "protein_position.1",
    "ml_target",
]


def main():

    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT_FILE)

    columns = list(df.columns)

    lines = []

    lines.append("=" * 70)
    lines.append("DEEPGENE ML STEP 29")
    lines.append("FINAL LEAKAGE & TARGET INDEPENDENCE AUDIT")
    lines.append("=" * 70)
    lines.append("")

    # --------------------------------------------------------
    # Dataset
    # --------------------------------------------------------

    lines.append("DATASET")
    lines.append("-" * 70)
    lines.append(f"Input: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Columns: {len(columns)}")
    lines.append(f"Target: {TARGET}")
    lines.append("")

    # --------------------------------------------------------
    # Target checks
    # --------------------------------------------------------

    target_exists = TARGET in columns

    lines.append("TARGET CHECK")
    lines.append("-" * 70)
    lines.append(f"Target exists: {target_exists}")

    if target_exists:
        target_missing = int(df[TARGET].isna().sum())
        target_unique = int(df[TARGET].nunique(dropna=True))

        lines.append(
            f"Target missing values: {target_missing}"
        )
        lines.append(
            f"Target unique classes: {target_unique}"
        )

        lines.append("Target distribution:")

        for cls, count in df[TARGET].value_counts().sort_index().items():
            lines.append(
                f"  {cls}: {count} "
                f"({count / len(df):.4f})"
            )
    else:
        target_missing = -1
        target_unique = 0

    lines.append("")

    # --------------------------------------------------------
    # Model feature-set audit
    # --------------------------------------------------------

    missing_features = [
        feature
        for feature in ALLOWED_FEATURES
        if feature not in columns
    ]

    lines.append("MODEL FEATURE SET CHECK")
    lines.append("-" * 70)
    lines.append(
        f"Expected model features: {len(ALLOWED_FEATURES)}"
    )
    lines.append(
        f"Expected metadata columns: {len(ALLOWED_METADATA)}"
    )

    lines.append(
        "Missing expected model features: "
        + (
            str(missing_features)
            if missing_features
            else "None"
        )
    )

    unexpected_columns = [
        column
        for column in columns
        if column != TARGET
        and column not in ALLOWED_FEATURES
        and column not in ALLOWED_METADATA
    ]

    lines.append(
        "Unexpected columns outside documented dataset structure: "
        + (
            str(unexpected_columns)
            if unexpected_columns
            else "None"
        )
    )

    lines.append("")

    # --------------------------------------------------------
    # Metadata audit
    # --------------------------------------------------------

    missing_metadata = [
        column
        for column in ALLOWED_METADATA
        if column not in columns
    ]

    lines.append("DOCUMENTED METADATA CHECK")
    lines.append("-" * 70)
    lines.append(
        "Expected metadata columns present: "
        + str(
            [
                column
                for column in ALLOWED_METADATA
                if column in columns
            ]
        )
    )

    lines.append(
        "Missing documented metadata columns: "
        + (
            str(missing_metadata)
            if missing_metadata
            else "None"
        )
    )

    lines.append("")

    # --------------------------------------------------------
    # Verify model predictor list exactly
    # --------------------------------------------------------

    model_feature_duplicates = (
        len(ALLOWED_FEATURES)
        != len(set(ALLOWED_FEATURES))
    )

    lines.append("MODEL PREDICTOR DEFINITION")
    lines.append("-" * 70)
    lines.append(
        "The supervised models are restricted to exactly "
        f"{len(ALLOWED_FEATURES)} biological features."
    )

    for feature in ALLOWED_FEATURES:
        lines.append(f"  {feature}")

    lines.append(
        f"Duplicate feature names in predictor list: "
        f"{model_feature_duplicates}"
    )

    lines.append("")

    # --------------------------------------------------------
    # Target duplication
    # --------------------------------------------------------

    exact_duplicates = []

    if target_exists:

        for feature in ALLOWED_FEATURES:

            if feature in df.columns:

                if df[feature].astype(str).equals(
                    df[TARGET].astype(str)
                ):
                    exact_duplicates.append(feature)

    lines.append("TARGET DUPLICATION AUDIT")
    lines.append("-" * 70)
    lines.append(
        "Model features exactly equal to target: "
        + (
            str(exact_duplicates)
            if exact_duplicates
            else "None"
        )
    )

    lines.append("")

    # --------------------------------------------------------
    # Functional-target metadata separation
    # --------------------------------------------------------

    target_related_metadata = [
        column
        for column in ALLOWED_METADATA
        if column in columns
        and (
            "functional" in column.lower()
            or column.lower() in {
                "ml_target",
                "protein_variant",
                "protein_position.1",
            }
        )
    ]

    lines.append("TARGET-RELATED METADATA SEPARATION")
    lines.append("-" * 70)
    lines.append(
        "Target-related metadata present in dataset: "
        + str(target_related_metadata)
    )
    lines.append(
        "These columns are NOT included in ALLOWED_FEATURES "
        "and therefore are not model predictors."
    )

    lines.append("")

    # --------------------------------------------------------
    # Numeric model feature audit
    # --------------------------------------------------------

    non_numeric_features = []

    for feature in ALLOWED_FEATURES:

        if feature in df.columns:

            if not pd.api.types.is_numeric_dtype(
                df[feature]
            ):
                non_numeric_features.append(feature)

    lines.append("MODEL FEATURE TYPE AUDIT")
    lines.append("-" * 70)
    lines.append(
        "Non-numeric model features: "
        + (
            str(non_numeric_features)
            if non_numeric_features
            else "None"
        )
    )

    lines.append("")

    # --------------------------------------------------------
    # Missing model feature values
    # --------------------------------------------------------

    missing_values = {}

    for feature in ALLOWED_FEATURES:

        if feature in df.columns:

            count = int(df[feature].isna().sum())

            if count > 0:
                missing_values[feature] = count

    lines.append("MODEL FEATURE MISSING-VALUE AUDIT")
    lines.append("-" * 70)

    if missing_values:
        for feature, count in missing_values.items():
            lines.append(
                f"{feature}: {count}"
            )
    else:
        lines.append(
            "No missing values detected in model features."
        )

    lines.append("")

    # --------------------------------------------------------
    # Identifier exclusion audit
    # --------------------------------------------------------

    identifier_columns = [
        "variant_id",
        "variant",
        "rsid",
        "hgvs",
        "gene",
        "chrom",
        "chromosome",
    ]

    present_identifiers = [
        column
        for column in identifier_columns
        if column in columns
    ]

    model_identifier_overlap = [
        column
        for column in ALLOWED_FEATURES
        if column in present_identifiers
    ]

    lines.append("IDENTIFIER EXCLUSION AUDIT")
    lines.append("-" * 70)
    lines.append(
        "Identifier-like dataset columns present: "
        + (
            str(present_identifiers)
            if present_identifiers
            else "None"
        )
    )
    lines.append(
        "Identifier columns included as model predictors: "
        + (
            str(model_identifier_overlap)
            if model_identifier_overlap
            else "None"
        )
    )

    lines.append("")

    # --------------------------------------------------------
    # Clinical / ClinVar feature exclusion
    # --------------------------------------------------------

    clinical_patterns = [
        "clinvar",
        "clinical",
        "pathogenic",
        "benign",
        "significance",
        "review",
        "evidence",
        "submitter",
        "phenotype",
        "inheritance",
        "disease",
        "condition",
        "classification",
        "conflict",
    ]

    clinical_columns = []

    for column in columns:

        lower = column.lower()

        if any(
            pattern in lower
            for pattern in clinical_patterns
        ):
            clinical_columns.append(column)

    clinical_predictor_overlap = [
        column
        for column in ALLOWED_FEATURES
        if column in clinical_columns
    ]

    lines.append("CLINICAL / CLINVAR FEATURE EXCLUSION")
    lines.append("-" * 70)
    lines.append(
        "Clinical/ClinVar-related dataset columns detected: "
        + (
            str(clinical_columns)
            if clinical_columns
            else "None"
        )
    )
    lines.append(
        "Clinical/ClinVar-related model predictors: "
        + (
            str(clinical_predictor_overlap)
            if clinical_predictor_overlap
            else "None"
        )
    )

    lines.append("")

    # --------------------------------------------------------
    # Final decision
    # --------------------------------------------------------

    failures = []

    if not target_exists:
        failures.append("Target column is missing.")

    if target_exists and target_missing != 0:
        failures.append("Target contains missing values.")

    if missing_features:
        failures.append(
            "One or more expected biological model features are missing."
        )

    if unexpected_columns:
        failures.append(
            "Unexpected columns exist outside the documented dataset structure."
        )

    if missing_metadata:
        failures.append(
            "Documented metadata columns are missing."
        )

    if model_feature_duplicates:
        failures.append(
            "Duplicate feature names exist in predictor definition."
        )

    if exact_duplicates:
        failures.append(
            "A model feature exactly duplicates the target."
        )

    if non_numeric_features:
        failures.append(
            "A model feature is non-numeric."
        )

    if missing_values:
        failures.append(
            "Missing values exist in model predictors."
        )

    if model_identifier_overlap:
        failures.append(
            "An identifier column is being used as a model predictor."
        )

    if clinical_predictor_overlap:
        failures.append(
            "A clinical/ClinVar-related column is being used as a model predictor."
        )

    lines.append("=" * 70)
    lines.append("FINAL LEAKAGE DECISION")
    lines.append("=" * 70)

    if failures:

        lines.append("STATUS: FAIL")
        lines.append("")
        lines.append("Failure flags:")

        for failure in failures:
            lines.append(
                f"- {failure}"
            )

    else:

        lines.append(
            "No structural leakage indicators detected."
        )
        lines.append("")
        lines.append(
            "The functional target is kept separate from the "
            "14 biological predictor features."
        )
        lines.append(
            "Metadata columns are present in the dataset but are "
            "excluded from the model predictor list."
        )
        lines.append(
            "Identifier columns are not model predictors."
        )
        lines.append(
            "Clinical/ClinVar-derived variables are not model predictors."
        )
        lines.append(
            "Phenotype and inheritance variables are not model predictors."
        )
        lines.append(
            "No model feature exactly duplicates the target."
        )
        lines.append("")
        lines.append("STATUS: PASS")

    lines.append("")
    lines.append("=" * 70)
    lines.append("LIMITATION")
    lines.append("=" * 70)
    lines.append(
        "This audit verifies the documented predictor set and structural "
        "separation of target, metadata, identifiers, and clinical variables."
    )
    lines.append(
        "It does not establish clinical validity or biological causality."
    )
    lines.append(
        "No model selection is performed in Step 29."
    )
    lines.append("=" * 70)

    REPORT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("\n".join(lines))


if __name__ == "__main__":
    main()