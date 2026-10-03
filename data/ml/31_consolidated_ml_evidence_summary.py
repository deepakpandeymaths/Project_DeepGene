from pathlib import Path

import pandas as pd


# ============================================================
# DEEPGENE ML STEP 31
# CONSOLIDATED ML EVIDENCE SUMMARY
# ============================================================

INPUT_FILE = Path(
    "data/processed/deepgene_ml_ready_v1.csv"
)

REPORT_FILE = Path(
    "data/analysis_results/31_consolidated_ml_evidence_summary.txt"
)

TARGET = "functional_target"

FEATURES = [
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


def main():

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    df = pd.read_csv(INPUT_FILE)

    lines = []

    # ========================================================
    # HEADER
    # ========================================================

    lines.append("=" * 78)
    lines.append("DEEPGENE ML STEP 31")
    lines.append("CONSOLIDATED ML EVIDENCE SUMMARY")
    lines.append("=" * 78)
    lines.append("")

    # ========================================================
    # PURPOSE
    # ========================================================

    lines.append("PURPOSE")
    lines.append("-" * 78)
    lines.append(
        "Consolidate the evidence generated across the DeepGene "
        "exploratory machine-learning workflow."
    )
    lines.append(
        "The purpose is to determine whether the current dataset "
        "supports selection of a final predictive model."
    )
    lines.append(
        "This step does not train a new model."
    )
    lines.append(
        "This step does not perform hyperparameter tuning."
    )
    lines.append(
        "This step does not perform feature selection."
    )
    lines.append(
        "This step does not provide clinical interpretation."
    )
    lines.append("")

    # ========================================================
    # DATASET
    # ========================================================

    lines.append("DATASET")
    lines.append("-" * 78)
    lines.append(
        f"Input file: {INPUT_FILE}"
    )
    lines.append(
        f"Rows: {len(df)}"
    )
    lines.append(
        f"Columns: {len(df.columns)}"
    )
    lines.append(
        f"Target: {TARGET}"
    )
    lines.append(
        f"Model features: {len(FEATURES)}"
    )
    lines.append("")

    # ========================================================
    # TARGET DISTRIBUTION
    # ========================================================

    counts = (
        df[TARGET]
        .value_counts()
        .sort_index()
    )

    lines.append("TARGET DISTRIBUTION")
    lines.append("-" * 78)

    for target_class, count in counts.items():

        fraction = count / len(df)

        lines.append(
            f"{target_class}: "
            f"{count} "
            f"({fraction:.4f})"
        )

    lines.append("")

    lines.append(
        "Smallest target class: "
        f"{counts.min()}"
    )

    lines.append(
        "Largest target class: "
        f"{counts.max()}"
    )

    lines.append("")

    # ========================================================
    # FEATURE AUDIT
    # ========================================================

    lines.append("MODEL FEATURE AUDIT")
    lines.append("-" * 78)

    missing_features = [
        feature
        for feature in FEATURES
        if feature not in df.columns
    ]

    if missing_features:

        lines.append(
            "ERROR: Missing model features:"
        )

        for feature in missing_features:
            lines.append(
                f"  - {feature}"
            )

    else:

        lines.append(
            "All 14 predefined biological model features are present."
        )

    numeric_failures = []

    for feature in FEATURES:

        if not pd.api.types.is_numeric_dtype(
            df[feature]
        ):
            numeric_failures.append(feature)

    if numeric_failures:

        lines.append(
            "ERROR: Non-numeric model features:"
        )

        for feature in numeric_failures:
            lines.append(
                f"  - {feature}"
            )

    else:

        lines.append(
            "All model features are numeric."
        )

    missing_values = (
        df[FEATURES]
        .isna()
        .sum()
        .sum()
    )

    lines.append(
        f"Missing model-feature values: "
        f"{missing_values}"
    )

    lines.append("")

    # ========================================================
    # EVIDENCE FROM PREVIOUS STEPS
    # ========================================================

    lines.append("CONSOLIDATED EVIDENCE")
    lines.append("-" * 78)

    lines.append(
        "Step 13 — Majority-class baseline:"
    )
    lines.append(
        "  Balanced accuracy: 0.3333"
    )
    lines.append(
        "  Macro F1: 0.2713"
    )
    lines.append(
        "  Majority class: Loss"
    )
    lines.append("")

    lines.append(
        "Step 15 — Logistic Regression:"
    )
    lines.append(
        "  Repeated-CV balanced accuracy: 0.2750 ± 0.0740"
    )
    lines.append(
        "  Repeated-CV macro F1: 0.2472 ± 0.0537"
    )
    lines.append("")

    lines.append(
        "Step 16 — Random Forest:"
    )
    lines.append(
        "  Repeated-CV balanced accuracy: 0.3032 ± 0.1663"
    )
    lines.append(
        "  Repeated-CV macro F1: 0.2718 ± 0.1274"
    )
    lines.append("")

    lines.append(
        "Step 17 — Initial model comparison:"
    )
    lines.append(
        "  Neither supervised model showed a robust advantage "
        "over the majority baseline."
    )
    lines.append("")

    lines.append(
        "Step 18 — Repeated stratified CV:"
    )
    lines.append(
        "  Logistic Regression remained below the majority baseline."
    )
    lines.append(
        "  Random Forest remained below the majority baseline "
        "on average."
    )
    lines.append(
        "  Results were unstable across fold assignments."
    )
    lines.append("")

    lines.append(
        "Step 19 — Feature class-separability audit:"
    )
    lines.append(
        "  delta_hydrophobicity showed an exploratory "
        "class-distribution difference."
    )
    lines.append(
        "  This did not establish predictive utility."
    )
    lines.append("")

    lines.append(
        "Step 20 — Feature-importance stability:"
    )
    lines.append(
        "  Feature importance was unstable across repeated "
        "validation folds."
    )
    lines.append(
        "  No feature was selected or removed."
    )
    lines.append("")

    lines.append(
        "Step 21 — Learning-curve analysis:"
    )
    lines.append(
        "  No clear pattern showed improved performance "
        "with additional training data."
    )
    lines.append(
        "  Dataset size remains a major limitation."
    )
    lines.append("")

    lines.append(
        "Step 22 — Controlled uncertainty analysis:"
    )
    lines.append(
        "  Logistic Regression mean balanced accuracy: 0.2750"
    )
    lines.append(
        "  Random Forest mean balanced accuracy: 0.3032"
    )
    lines.append(
        "  Both showed substantial fold-level uncertainty."
    )
    lines.append("")

    lines.append(
        "Step 23 — Out-of-fold prediction audit:"
    )
    lines.append(
        "  Predictions showed some repeated stability, "
        "but stable predictions were not equivalent to correct predictions."
    )
    lines.append(
        "  Logistic Regression repeatedly misclassified Gain variants."
    )
    lines.append(
        "  Random Forest showed substantial Gain/Mixed confusion."
    )
    lines.append("")

    lines.append(
        "Step 24 — Error-pattern analysis:"
    )
    lines.append(
        "  Major repeated error directions included "
        "Gain→Loss, Loss→Gain, Loss→Mixed, and Mixed→Loss."
    )
    lines.append(
        "  Error persistence was high for minority classes."
    )
    lines.append("")

    lines.append(
        "Step 25 — Feature correlation audit:"
    )
    lines.append(
        "  No pairwise feature correlation reached |r| >= 0.80."
    )
    lines.append(
        "  No feature was removed."
    )
    lines.append("")

    lines.append(
        "Step 26 — Feature variance audit:"
    )
    lines.append(
        "  same_amino_acid was constant."
    )
    lines.append(
        "  No additional near-zero-variance feature "
        "was identified."
    )
    lines.append("")

    lines.append(
        "Step 27 — Feature scale audit:"
    )
    lines.append(
        "  Biological features have substantially different numerical scales."
    )
    lines.append(
        "  Logistic Regression appropriately uses "
        "fold-contained standardization."
    )
    lines.append("")

    lines.append(
        "Step 28 — Class-weight sensitivity:"
    )
    lines.append(
        "  Logistic Regression performance decreased "
        "under class weighting."
    )
    lines.append(
        "  Random Forest mean performance was similar "
        "with and without class weighting, but balanced "
        "Random Forest remained highly variable."
    )
    lines.append("")

    lines.append(
        "Step 29 — Final leakage and target-independence audit:"
    )
    lines.append(
        "  STATUS PASS."
    )
    lines.append(
        "  The documented predictor set is structurally "
        "separated from the functional target."
    )
    lines.append("")

    lines.append(
        "Step 30 — Target permutation sanity test:"
    )
    lines.append(
        "  Logistic Regression real-label balanced accuracy: 0.2750"
    )
    lines.append(
        "  Logistic Regression permutation mean: 0.3394"
    )
    lines.append(
        "  95/100 randomized permutations reached or exceeded "
        "the real-label score."
    )
    lines.append(
        "  Random Forest real-label balanced accuracy: 0.3012"
    )
    lines.append(
        "  Random Forest permutation mean: 0.3390"
    )
    lines.append(
        "  78/100 randomized permutations reached or exceeded "
        "the real-label score."
    )
    lines.append(
        "  The observed scores were therefore not unusually "
        "high relative to randomized-target performance."
    )
    lines.append("")

    # ========================================================
    # OVERALL ASSESSMENT
    # ========================================================

    lines.append("=" * 78)
    lines.append("OVERALL ML ASSESSMENT")
    lines.append("=" * 78)

    lines.append(
        "1. The independent functional-effect dataset is suitable "
        "for exploratory supervised analysis."
    )

    lines.append(
        "2. Structural leakage checks passed."
    )

    lines.append(
        "3. The current biological feature set does not produce "
        "robust predictive performance above the majority-class "
        "baseline."
    )

    lines.append(
        "4. Repeated cross-validation demonstrates substantial "
        "performance variability."
    )

    lines.append(
        "5. Feature importance is unstable."
    )

    lines.append(
        "6. Error analysis shows persistent difficulty with "
        "the minority functional classes."
    )

    lines.append(
        "7. The target-permutation test does not provide evidence "
        "that the observed model scores are unusually high relative "
        "to randomized target labels."
    )

    lines.append(
        "8. The Gain class contains only 4 variants, creating a "
        "major limitation for supervised classification."
    )

    lines.append(
        "9. No model has demonstrated sufficient evidence for "
        "selection as a final DeepGene predictive model."
    )

    lines.append("")

    # ========================================================
    # MODEL DECISION
    # ========================================================

    lines.append("=" * 78)
    lines.append("MODEL DECISION")
    lines.append("=" * 78)

    lines.append(
        "FINAL MODEL SELECTED: NO"
    )

    lines.append(
        "REASON:"
    )

    lines.append(
        "The available evidence does not support selecting "
        "Logistic Regression or Random Forest as a final "
        "predictive model."
    )

    lines.append(
        "Both models perform around or below the majority-class "
        "balanced-accuracy baseline, show substantial uncertainty, "
        "and do not outperform randomized-target performance "
        "convincingly."
    )

    lines.append("")

    # ========================================================
    # WHAT THIS MEANS
    # ========================================================

    lines.append("=" * 78)
    lines.append("SCIENTIFIC INTERPRETATION")
    lines.append("=" * 78)

    lines.append(
        "The current result should be interpreted as an exploratory "
        "negative/insufficient predictive finding."
    )

    lines.append(
        "It does not demonstrate that the biological features "
        "have no relationship with SCN1A functional effect."
    )

    lines.append(
        "It demonstrates that this 51-variant dataset and the "
        "current 14-feature representation do not provide sufficient "
        "evidence for a reliable three-class predictive model."
    )

    lines.append(
        "The result should not be interpreted as clinical prediction, "
        "pathogenicity classification, or disease-risk prediction."
    )

    lines.append("")

    # ========================================================
    # LIMITATIONS
    # ========================================================

    lines.append("=" * 78)
    lines.append("KEY LIMITATIONS")
    lines.append("=" * 78)

    limitations = [
        "Only 51 variants are available for the primary supervised task.",
        "The Gain class contains only 4 variants.",
        "The functional dataset originates from a limited experimental source.",
        "Repeated cross-validation remains unstable.",
        "Feature importance is not stable across validation assignments.",
        "The current feature representation is limited to amino-acid-level biological properties and protein position.",
        "No external independent functional benchmark has been evaluated in this workflow.",
        "No hyperparameter tuning was used to rescue model performance.",
        "No clinical validity or clinical utility has been established.",
    ]

    for number, limitation in enumerate(
        limitations,
        start=1,
    ):

        lines.append(
            f"{number}. {limitation}"
        )

    lines.append("")

    # ========================================================
    # NEXT-STAGE DECISION
    # ========================================================

    lines.append("=" * 78)
    lines.append("NEXT-STAGE DECISION")
    lines.append("=" * 78)

    lines.append(
        "The project should not proceed to deployment or clinical "
        "prediction using the current supervised models."
    )

    lines.append(
        "The next stage should focus on documenting the negative/limited "
        "ML finding and determining whether additional independent "
        "functional data or scientifically justified features are "
        "available."
    )

    lines.append(
        "Any future model development should preserve the independent "
        "functional target and the leakage controls established "
        "in Steps 29 and 30."
    )

    lines.append("")

    # ========================================================
    # STATUS
    # ========================================================

    lines.append("=" * 78)
    lines.append("STATUS: PASS")
    lines.append("=" * 78)

    REPORT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("\n".join(lines))


if __name__ == "__main__":
    main()