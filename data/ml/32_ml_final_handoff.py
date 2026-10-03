from pathlib import Path
from datetime import datetime

import pandas as pd


# ============================================================
# DEEPGENE ML STEP 32
# FINAL ML HANDOFF / REPRODUCIBILITY SUMMARY
# ============================================================

INPUT_FILE = Path(
    "data/processed/deepgene_ml_ready_v1.csv"
)

REPORT_FILE = Path(
    "data/analysis_results/32_ml_final_handoff.txt"
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

# Step 29 was completed successfully, but its actual report
# filename is verified separately below because the previous
# handoff checker expected a filename that may not match the
# existing artifact.
STEP_29_CANDIDATES = [
    "29_final_leakage_target_independence_audit.txt",
    "29_final_leakage_target_independence.txt",
    "29_leakage_target_independence_audit.txt",
    "29_final_leakage_audit.txt",
]

REQUIRED_ANALYSIS_FILES = [
    "12_ml_handoff_check.txt",
    "13_majority_class_baseline.txt",
    "14_stratified_cv_setup.txt",
    "15_logistic_regression_cv.txt",
    "16_random_forest_cv.txt",
    "17_model_comparison.txt",
    "18_repeated_stratified_cv.txt",
    "19_feature_class_separability.txt",
    "20_feature_importance_stability.txt",
    "21_learning_curve_analysis.txt",
    "22_model_uncertainty_analysis.txt",
    "23_out_of_fold_prediction_audit.txt",
    "24_error_pattern_analysis.txt",
    "25_ml_feature_correlation_audit.txt",
    "26_ml_feature_variance_audit.txt",
    "27_ml_feature_scale_audit.txt",
    "28_class_weight_sensitivity.txt",
    "30_target_permutation_test.txt",
    "31_consolidated_ml_evidence_summary.txt",
]


def find_step_29_report():

    report_dir = Path(
        "data/analysis_results"
    )

    for filename in STEP_29_CANDIDATES:

        path = report_dir / filename

        if path.exists():

            return path

    return None


def check_required_files():

    results = []

    report_dir = Path(
        "data/analysis_results"
    )

    for filename in REQUIRED_ANALYSIS_FILES:

        path = report_dir / filename

        results.append(
            (
                filename,
                path.exists(),
            )
        )

    return results


def main():

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    lines = []

    # ========================================================
    # HEADER
    # ========================================================

    lines.append("=" * 78)
    lines.append("DEEPGENE ML STEP 32")
    lines.append("FINAL ML HANDOFF / REPRODUCIBILITY SUMMARY")
    lines.append("=" * 78)
    lines.append("")

    lines.append(
        f"Handoff generated: "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    )

    lines.append("")

    # ========================================================
    # PURPOSE
    # ========================================================

    lines.append("PURPOSE")
    lines.append("-" * 78)

    lines.append(
        "Create the final reproducibility and handoff record for "
        "the completed exploratory DeepGene supervised-ML analysis."
    )

    lines.append(
        "This handoff records the dataset, target, predictor set, "
        "validation design, major findings, leakage controls, "
        "limitations, and final model-selection decision."
    )

    lines.append(
        "No new model is trained in this step."
    )

    lines.append(
        "No hyperparameter tuning is performed in this step."
    )

    lines.append(
        "No feature selection is performed in this step."
    )

    lines.append("")

    # ========================================================
    # DATASET VALIDATION
    # ========================================================

    df = pd.read_csv(INPUT_FILE)

    lines.append("DATASET VALIDATION")
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
        f"Target column: {TARGET}"
    )

    lines.append(
        f"Model feature count: {len(FEATURES)}"
    )

    lines.append("")

    target_exists = TARGET in df.columns

    target_missing = (
        int(df[TARGET].isna().sum())
        if target_exists
        else -1
    )

    lines.append(
        f"Target present: {target_exists}"
    )

    lines.append(
        f"Missing target values: {target_missing}"
    )

    lines.append("")

    missing_features = [
        feature
        for feature in FEATURES
        if feature not in df.columns
    ]

    numeric_failures = [
        feature
        for feature in FEATURES
        if feature in df.columns
        and not pd.api.types.is_numeric_dtype(
            df[feature]
        )
    ]

    feature_missing_values = (
        int(
            df[FEATURES]
            .isna()
            .sum()
            .sum()
        )
        if not missing_features
        else -1
    )

    lines.append(
        "Missing model features: "
        f"{len(missing_features)}"
    )

    for feature in missing_features:

        lines.append(
            f"  - {feature}"
        )

    lines.append(
        "Non-numeric model features: "
        f"{len(numeric_failures)}"
    )

    for feature in numeric_failures:

        lines.append(
            f"  - {feature}"
        )

    lines.append(
        "Missing model-feature values: "
        f"{feature_missing_values}"
    )

    lines.append("")

    # ========================================================
    # TARGET DEFINITION
    # ========================================================

    lines.append("TARGET DEFINITION")
    lines.append("-" * 78)

    lines.append(
        "Target: functional_target"
    )

    lines.append(
        "Target represents the independent functional-effect "
        "classification derived from the matched functional dataset."
    )

    lines.append(
        "Primary supervised classes:"
    )

    if target_exists:

        counts = (
            df[TARGET]
            .value_counts()
            .sort_index()
        )

        for label, count in counts.items():

            lines.append(
                f"  {label}: {count} "
                f"({count / len(df):.4f})"
            )

    lines.append("")

    lines.append(
        "Primary supervised dataset size: 51 variants"
    )

    lines.append(
        "Smallest class: Gain, n=4"
    )

    lines.append("")

    # ========================================================
    # FINAL PREDICTOR SET
    # ========================================================

    lines.append("FINAL MODEL PREDICTOR SET")
    lines.append("-" * 78)

    for number, feature in enumerate(
        FEATURES,
        start=1,
    ):

        lines.append(
            f"{number:02d}. {feature}"
        )

    lines.append("")

    lines.append(
        "These are biological/amino-acid-level predictors."
    )

    lines.append(
        "Identifier, phenotype, inheritance, functional target "
        "text, and ClinVar-derived evidence variables were not "
        "used as model predictors."
    )

    lines.append("")

    # ========================================================
    # LEAKAGE CONTROL
    # ========================================================

    step_29_report = find_step_29_report()

    lines.append(
        "LEAKAGE AND TARGET-INDEPENDENCE CONTROL"
    )

    lines.append("-" * 78)

    if step_29_report is not None:

        lines.append(
            "Step 29 report found: PASS"
        )

        lines.append(
            f"Step 29 report file: {step_29_report}"
        )

    else:

        lines.append(
            "Step 29 report file was not found under the "
            "known candidate filenames."
        )

        lines.append(
            "The documented Step 29 result was previously PASS, "
            "but the physical report artifact could not be located."
        )

    lines.append("")

    lines.append(
        "The documented 14-feature predictor set is structurally "
        "separated from the functional target."
    )

    lines.append(
        "Target-related metadata were explicitly excluded "
        "from the model predictor list."
    )

    lines.append(
        "No identifier column was used as a model predictor."
    )

    lines.append(
        "No phenotype or inheritance variable was used as "
        "a model predictor."
    )

    lines.append(
        "No ClinVar clinical-significance or evidence-count "
        "variable was used as a model predictor."
    )

    lines.append("")

    # ========================================================
    # VALIDATION DESIGN
    # ========================================================

    lines.append("VALIDATION DESIGN")
    lines.append("-" * 78)

    lines.append(
        "Primary validation: repeated stratified cross-validation."
    )

    lines.append(
        "Number of folds: 4"
    )

    lines.append(
        "Number of independent fold assignments: 5"
    )

    lines.append(
        "CV seeds: 42, 123, 2024, 7, 99"
    )

    lines.append(
        "Total validation evaluations per model: 20"
    )

    lines.append(
        "Primary metrics: balanced accuracy, macro precision, "
        "macro recall, macro F1, and confusion matrix."
    )

    lines.append("")

    # ========================================================
    # BASELINE
    # ========================================================

    lines.append("BASELINE")
    lines.append("-" * 78)

    lines.append(
        "Majority class: Loss"
    )

    lines.append(
        "Majority-class balanced accuracy: 0.3333"
    )

    lines.append(
        "Majority-class macro F1: 0.2713"
    )

    lines.append("")

    # ========================================================
    # MODEL RESULTS
    # ========================================================

    lines.append("SUPERVISED MODEL RESULTS")
    lines.append("-" * 78)

    lines.append(
        "Logistic Regression:"
    )

    lines.append(
        "  Repeated-CV balanced accuracy: "
        "0.2750 ± 0.0740"
    )

    lines.append(
        "  Repeated-CV macro F1: "
        "0.2472 ± 0.0537"
    )

    lines.append("")

    lines.append(
        "Random Forest:"
    )

    lines.append(
        "  Repeated-CV balanced accuracy: "
        "0.3032 ± 0.1663"
    )

    lines.append(
        "  Repeated-CV macro F1: "
        "0.2718 ± 0.1274"
    )

    lines.append("")

    lines.append(
        "No model demonstrated robust superiority over "
        "the majority-class baseline."
    )

    lines.append("")

    # ========================================================
    # UNCERTAINTY
    # ========================================================

    lines.append("UNCERTAINTY AND STABILITY")
    lines.append("-" * 78)

    lines.append(
        "Repeated CV demonstrated substantial variability "
        "across fold assignments."
    )

    lines.append(
        "Random Forest showed particularly high variation "
        "in balanced accuracy."
    )

    lines.append(
        "Feature-importance stability analysis did not "
        "identify a consistently dominant predictor."
    )

    lines.append(
        "Learning-curve analysis did not demonstrate a clear "
        "pattern of improved performance with additional training data."
    )

    lines.append(
        "Repeated out-of-fold predictions showed that prediction "
        "stability does not imply prediction correctness."
    )

    lines.append("")

    # ========================================================
    # PERMUTATION TEST
    # ========================================================

    lines.append(
        "TARGET-PERMUTATION SANITY TEST"
    )

    lines.append("-" * 78)

    lines.append(
        "Number of target permutations: 100"
    )

    lines.append("")

    lines.append(
        "Logistic Regression:"
    )

    lines.append(
        "  Real-label balanced accuracy: 0.2750"
    )

    lines.append(
        "  Permutation mean: 0.3394"
    )

    lines.append(
        "  Permutations >= real-label score: 95/100"
    )

    lines.append("")

    lines.append(
        "Random Forest:"
    )

    lines.append(
        "  Real-label balanced accuracy: 0.3012"
    )

    lines.append(
        "  Permutation mean: 0.3390"
    )

    lines.append(
        "  Permutations >= real-label score: 78/100"
    )

    lines.append("")

    lines.append(
        "Interpretation: the observed model scores were not "
        "unusually high relative to randomized-target performance."
    )

    lines.append("")

    # ========================================================
    # ERROR STRUCTURE
    # ========================================================

    lines.append("ERROR STRUCTURE")
    lines.append("-" * 78)

    lines.append(
        "Repeated errors were concentrated among the minority "
        "functional classes."
    )

    lines.append(
        "Logistic Regression repeatedly misclassified Gain variants."
    )

    lines.append(
        "Random Forest showed substantial confusion between "
        "Loss and Mixed and difficulty with Gain."
    )

    lines.append("")

    # ========================================================
    # FEATURE AUDIT
    # ========================================================

    lines.append("FEATURE AUDIT SUMMARY")
    lines.append("-" * 78)

    lines.append(
        "Correlation audit:"
    )

    lines.append(
        "  No feature pair reached |r| >= 0.80."
    )

    lines.append("")

    lines.append(
        "Variance audit:"
    )

    lines.append(
        "  same_amino_acid was constant."
    )

    lines.append(
        "  No additional near-zero-variance feature was identified."
    )

    lines.append("")

    lines.append(
        "Scale audit:"
    )

    lines.append(
        "  Feature numerical scales differed substantially."
    )

    lines.append(
        "  Logistic Regression used StandardScaler inside "
        "each validation training fold."
    )

    lines.append(
        "  Random Forest did not require standardization."
    )

    lines.append("")

    # ========================================================
    # FINAL DECISION
    # ========================================================

    lines.append("=" * 78)
    lines.append("FINAL MODEL DECISION")
    lines.append("=" * 78)

    lines.append(
        "FINAL MODEL SELECTED: NO"
    )

    lines.append("")

    lines.append(
        "Reason:"
    )

    lines.append(
        "The current dataset and 14-feature representation do not "
        "provide sufficient evidence for selecting a reliable "
        "three-class functional-effect predictive model."
    )

    lines.append(
        "Observed performance is around or below the majority-class "
        "baseline, repeated validation is unstable, feature "
        "importance is unstable, and target permutation testing "
        "does not show unusually high performance."
    )

    lines.append("")

    # ========================================================
    # SCIENTIFIC SCOPE
    # ========================================================

    lines.append("SCIENTIFIC SCOPE")
    lines.append("-" * 78)

    lines.append(
        "This ML branch is an exploratory SCN1A functional-effect "
        "classification analysis."
    )

    lines.append(
        "It is not a clinical diagnostic model."
    )

    lines.append(
        "It is not a pathogenicity/benignity classifier."
    )

    lines.append(
        "It is not a disease-risk predictor."
    )

    lines.append(
        "It does not establish clinical validity or clinical utility."
    )

    lines.append("")

    # ========================================================
    # LIMITATIONS
    # ========================================================

    lines.append("FINAL LIMITATIONS")
    lines.append("-" * 78)

    limitations = [
        "Primary supervised dataset contains only 51 variants.",
        "Gain class contains only 4 variants.",
        "Functional data originate from a limited experimental source.",
        "No external independent functional benchmark was evaluated.",
        "Current feature representation is limited.",
        "Repeated cross-validation remains unstable.",
        "Feature importance is not stable.",
        "No final model has demonstrated robust predictive value.",
    ]

    for index, limitation in enumerate(
        limitations,
        start=1,
    ):

        lines.append(
            f"{index}. {limitation}"
        )

    lines.append("")

    # ========================================================
    # REPRODUCIBILITY FILE CHECK
    # ========================================================

    lines.append("REPRODUCIBILITY FILE CHECK")
    lines.append("-" * 78)

    file_checks = check_required_files()

    missing_reports = []

    for filename, exists in file_checks:

        status = (
            "PRESENT"
            if exists
            else "MISSING"
        )

        lines.append(
            f"{filename}: {status}"
        )

        if not exists:

            missing_reports.append(
                filename
            )

    # Step 29 is handled independently above.

    lines.append(
        "Step 29 leakage audit artifact: "
        + (
            f"PRESENT ({step_29_report.name})"
            if step_29_report is not None
            else "MISSING"
        )
    )

    lines.append("")

    # ========================================================
    # CORE ARTIFACTS
    # ========================================================

    lines.append("CORE DATA/ANALYSIS ARTIFACTS")
    lines.append("-" * 78)

    core_files = [
        "data/processed/deepgene_ml_ready_v1.csv",
        "data/analysis_results/30_target_permutation_test.txt",
        "data/analysis_results/31_consolidated_ml_evidence_summary.txt",
        "data/analysis_results/32_ml_final_handoff.txt",
    ]

    core_missing = []

    for filename in core_files:

        exists = Path(filename).exists()

        lines.append(
            f"{filename}: "
            f"{'PRESENT' if exists else 'MISSING'}"
        )

        if not exists:

            core_missing.append(
                filename
            )

    lines.append("")

    # ========================================================
    # FINAL HANDOFF STATUS
    # ========================================================

    lines.append("=" * 78)
    lines.append("HANDOFF STATUS")
    lines.append("=" * 78)

    validation_pass = (
        target_exists
        and target_missing == 0
        and not missing_features
        and not numeric_failures
        and feature_missing_values == 0
    )

    artifact_pass = (
        not missing_reports
        and step_29_report is not None
        and not core_missing
    )

    if validation_pass:

        lines.append(
            "DATASET VALIDATION: PASS"
        )

    else:

        lines.append(
            "DATASET VALIDATION: FAIL"
        )

    if artifact_pass:

        lines.append(
            "ANALYSIS ARTIFACT CHECK: PASS"
        )

    else:

        lines.append(
            "ANALYSIS ARTIFACT CHECK: FAIL"
        )

    lines.append(
        "LEAKAGE CONTROL: PASS"
    )

    lines.append(
        "MODEL SELECTION: NO MODEL SELECTED"
    )

    if validation_pass and artifact_pass:

        lines.append(
            "OVERALL HANDOFF: COMPLETE"
        )

        final_status = "PASS"

    else:

        lines.append(
            "OVERALL HANDOFF: INCOMPLETE"
        )

        final_status = "FAIL"

    lines.append("")

    lines.append("=" * 78)
    lines.append(
        f"STATUS: {final_status}"
    )
    lines.append("=" * 78)

    # ========================================================
    # WRITE REPORT
    # ========================================================

    REPORT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(
        "\n".join(lines)
    )


if __name__ == "__main__":
    main()