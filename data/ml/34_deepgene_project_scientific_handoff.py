from pathlib import Path
from datetime import datetime

import pandas as pd


# ============================================================
# DEEPGENE STEP 34
# PROJECT SCIENTIFIC HANDOFF
# ============================================================

# ------------------------------------------------------------
# CORE DATASET
# ------------------------------------------------------------

V1_DATASET = Path(
    "data/processed/deepgene_v1_final.csv"
)

ML_DATASET = Path(
    "data/processed/deepgene_ml_ready_v1.csv"
)

# ------------------------------------------------------------
# ANALYSIS REPORTS
# ------------------------------------------------------------

ANALYSIS_DIR = Path(
    "data/analysis_results"
)

# Major completed project reports
REQUIRED_REPORTS = [
    "42_final_ml_design.txt",
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
    "29_final_leakage_audit.txt",
    "30_target_permutation_test.txt",
    "31_consolidated_ml_evidence_summary.txt",
    "32_ml_final_handoff.txt",
    "33_ml_project_integration_handoff.txt",
]

OUTPUT_FILE = (
    ANALYSIS_DIR
    / "34_deepgene_project_scientific_handoff.txt"
)


# ============================================================
# HELPER
# ============================================================

def status(path):
    return "PRESENT" if path.exists() else "MISSING"


# ============================================================
# MAIN
# ============================================================

def main():

    ANALYSIS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------------
    # LOAD CORE DATASETS
    # --------------------------------------------------------

    v1_df = pd.read_csv(V1_DATASET)
    ml_df = pd.read_csv(ML_DATASET)

    # --------------------------------------------------------
    # BASIC DATASET INFORMATION
    # --------------------------------------------------------

    v1_rows = len(v1_df)
    v1_columns = len(v1_df.columns)

    ml_rows = len(ml_df)
    ml_columns = len(ml_df.columns)

    target = "functional_target"

    ml_features = [
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

    # --------------------------------------------------------
    # TARGET CHECK
    # --------------------------------------------------------

    target_present = target in ml_df.columns

    if target_present:

        target_missing = int(
            ml_df[target].isna().sum()
        )

        target_counts = (
            ml_df[target]
            .value_counts()
            .sort_index()
        )

    else:

        target_missing = -1
        target_counts = pd.Series(dtype=int)

    # --------------------------------------------------------
    # FEATURE CHECK
    # --------------------------------------------------------

    missing_features = [
        feature
        for feature in ml_features
        if feature not in ml_df.columns
    ]

    if not missing_features:

        missing_feature_values = int(
            ml_df[ml_features]
            .isna()
            .sum()
            .sum()
        )

        non_numeric_features = [
            feature
            for feature in ml_features
            if not pd.api.types.is_numeric_dtype(
                ml_df[feature]
            )
        ]

    else:

        missing_feature_values = -1
        non_numeric_features = []

    # --------------------------------------------------------
    # REPORT AVAILABILITY
    # --------------------------------------------------------

    report_status = {
        report: status(
            ANALYSIS_DIR / report
        )
        for report in REQUIRED_REPORTS
    }

    all_reports_present = all(
        value == "PRESENT"
        for value in report_status.values()
    )

    dataset_status = (
        V1_DATASET.exists()
        and ML_DATASET.exists()
    )

    ml_validation = (
        target_present
        and target_missing == 0
        and len(missing_features) == 0
        and missing_feature_values == 0
        and len(non_numeric_features) == 0
    )

    # --------------------------------------------------------
    # BUILD REPORT
    # --------------------------------------------------------

    lines = []

    lines.append("=" * 80)
    lines.append(
        "DEEPGENE STEP 34"
    )
    lines.append(
        "PROJECT SCIENTIFIC HANDOFF"
    )
    lines.append("=" * 80)
    lines.append("")

    lines.append(
        "Generated: "
        + datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
    )

    lines.append("")

    # ========================================================
    # PURPOSE
    # ========================================================

    lines.append("PURPOSE")
    lines.append("-" * 80)

    lines.append(
        "Create the final project-level scientific handoff "
        "for the completed DeepGene V1 analysis and its "
        "exploratory machine-learning branch."
    )

    lines.append(
        "This step does not train a model."
    )

    lines.append(
        "This step does not perform model tuning."
    )

    lines.append(
        "This step does not perform feature selection."
    )

    lines.append(
        "This step does not create a clinical prediction system."
    )

    lines.append("")

    # ========================================================
    # PROJECT STATUS
    # ========================================================

    lines.append("=" * 80)
    lines.append("PROJECT STATUS")
    lines.append("=" * 80)

    lines.append(
        "DeepGene V1 dataset construction: COMPLETE"
    )

    lines.append(
        "Evidence consistency auditing: COMPLETE"
    )

    lines.append(
        "Phenotype normalization/auditing: COMPLETE"
    )

    lines.append(
        "Feature redundancy auditing: COMPLETE"
    )

    lines.append(
        "Independent functional-data investigation: COMPLETE"
    )

    lines.append(
        "Independent functional ML dataset construction: COMPLETE"
    )

    lines.append(
        "Exploratory supervised ML analysis: COMPLETE"
    )

    lines.append(
        "Leakage/target-independence auditing: COMPLETE"
    )

    lines.append(
        "Permutation sanity testing: COMPLETE"
    )

    lines.append(
        "ML reproducibility handoff: COMPLETE"
    )

    lines.append(
        "Project scientific handoff: CURRENT STEP"
    )

    lines.append("")

    # ========================================================
    # V1 DATASET
    # ========================================================

    lines.append("DEEPGENE V1 DATASET")
    lines.append("-" * 80)

    lines.append(
        f"Dataset: {V1_DATASET}"
    )

    lines.append(
        f"Rows: {v1_rows}"
    )

    lines.append(
        f"Columns: {v1_columns}"
    )

    lines.append(
        "The V1 dataset is the canonical project-level "
        "variant dataset used for the completed preprocessing "
        "and evidence-analysis workflow."
    )

    lines.append("")

    # ========================================================
    # DATA QUALITY / EVIDENCE
    # ========================================================

    lines.append("DATA AND EVIDENCE AUDIT STATUS")
    lines.append("-" * 80)

    lines.append(
        "Evidence consistency audit: PASS"
    )

    lines.append(
        "Evidence-count consistency audit: PASS"
    )

    lines.append(
        "Phenotype evidence audit: PASS"
    )

    lines.append(
        "Phenotype normalization validation: PASS"
    )

    lines.append(
        "Feature redundancy audit: PASS"
    )

    lines.append("")

    lines.append(
        "The completed audits did not identify the documented "
        "structural inconsistencies that were specifically tested."
    )

    lines.append("")

    # ========================================================
    # INDEPENDENT FUNCTIONAL DATA
    # ========================================================

    lines.append(
        "INDEPENDENT FUNCTIONAL-EFFECT ANALYSIS"
    )
    lines.append("-" * 80)

    lines.append(
        "The functional-effect analysis used an independent "
        "functional dataset rather than using ClinVar "
        "pathogenicity labels as the supervised target."
    )

    lines.append(
        "The matched functional dataset contained "
        "56 matched functional variants."
    )

    lines.append(
        "The primary supervised ML subset contained "
        "51 variants after excluding functional classes "
        "that were not suitable for the primary three-class "
        "supervised analysis."
    )

    lines.append(
        "Primary supervised classes:"
    )

    lines.append(
        "  Loss: 35"
    )

    lines.append(
        "  Mixed: 12"
    )

    lines.append(
        "  Gain: 4"
    )

    lines.append("")

    # ========================================================
    # ML DATASET
    # ========================================================

    lines.append(
        "FINAL EXPLORATORY ML DATASET"
    )
    lines.append("-" * 80)

    lines.append(
        f"Dataset: {ML_DATASET}"
    )

    lines.append(
        f"Rows: {ml_rows}"
    )

    lines.append(
        f"Columns: {ml_columns}"
    )

    lines.append(
        f"Target: {target}"
    )

    lines.append(
        f"Predictor count: {len(ml_features)}"
    )

    lines.append("")

    if target_present:

        lines.append(
            "Target distribution:"
        )

        for label, count in target_counts.items():

            lines.append(
                f"  {label}: {count} "
                f"({count / len(ml_df):.4f})"
            )

    lines.append("")

    lines.append(
        "Smallest supervised class: Gain, n=4"
    )

    lines.append("")

    # ========================================================
    # ML FEATURE SET
    # ========================================================

    lines.append(
        "FINAL BIOLOGICAL PREDICTOR SET"
    )
    lines.append("-" * 80)

    for i, feature in enumerate(
        ml_features,
        start=1,
    ):

        lines.append(
            f"{i:02d}. {feature}"
        )

    lines.append("")

    lines.append(
        "The final predictor representation consists of "
        "protein-position and amino-acid physicochemical "
        "features."
    )

    lines.append(
        "Identifier, phenotype, inheritance, functional target "
        "text, and ClinVar-derived evidence variables were "
        "excluded from the model predictor set."
    )

    lines.append("")

    # ========================================================
    # ML VALIDATION
    # ========================================================

    lines.append(
        "ML DATASET VALIDATION"
    )
    lines.append("-" * 80)

    lines.append(
        f"Target present: {target_present}"
    )

    lines.append(
        f"Missing target values: {target_missing}"
    )

    lines.append(
        f"Missing predictor columns: "
        f"{len(missing_features)}"
    )

    lines.append(
        f"Missing predictor values: "
        f"{missing_feature_values}"
    )

    lines.append(
        f"Non-numeric predictors: "
        f"{len(non_numeric_features)}"
    )

    lines.append("")

    # ========================================================
    # LEAKAGE CONTROL
    # ========================================================

    lines.append(
        "LEAKAGE AND TARGET-INDEPENDENCE"
    )
    lines.append("-" * 80)

    lines.append(
        "Step 29 leakage audit: PASS"
    )

    lines.append(
        "The documented predictor set is structurally "
        "separated from the functional target."
    )

    lines.append(
        "Target-derived functional metadata were excluded "
        "from model predictors."
    )

    lines.append(
        "Identifier variables were not used as predictors."
    )

    lines.append(
        "Phenotype and inheritance were not used as predictors."
    )

    lines.append(
        "ClinVar clinical-significance and evidence-count "
        "variables were not used as predictors."
    )

    lines.append("")

    # ========================================================
    # ML VALIDATION DESIGN
    # ========================================================

    lines.append(
        "ML VALIDATION DESIGN"
    )
    lines.append("-" * 80)

    lines.append(
        "Validation strategy: repeated stratified "
        "cross-validation"
    )

    lines.append(
        "Folds per repeat: 4"
    )

    lines.append(
        "Independent fold assignments: 5"
    )

    lines.append(
        "Total validation evaluations per model: 20"
    )

    lines.append(
        "Seeds: 42, 123, 2024, 7, 99"
    )

    lines.append(
        "Primary metrics: balanced accuracy, macro precision, "
        "macro recall, macro F1, confusion matrix"
    )

    lines.append("")

    # ========================================================
    # ML RESULTS
    # ========================================================

    lines.append(
        "ML RESULTS"
    )
    lines.append("-" * 80)

    lines.append(
        "Majority-class baseline:"
    )

    lines.append(
        "  Balanced accuracy: 0.3333"
    )

    lines.append(
        "  Macro F1: 0.2713"
    )

    lines.append("")

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
        "No final predictive model was selected."
    )

    lines.append("")

    # ========================================================
    # UNCERTAINTY
    # ========================================================

    lines.append(
        "UNCERTAINTY AND STABILITY"
    )
    lines.append("-" * 80)

    lines.append(
        "Repeated CV showed substantial variability."
    )

    lines.append(
        "Random Forest showed particularly high variability "
        "in balanced accuracy."
    )

    lines.append(
        "Feature-importance stability did not identify a "
        "consistently dominant predictor."
    )

    lines.append(
        "Learning-curve analysis did not demonstrate a clear "
        "improvement pattern with additional training data."
    )

    lines.append(
        "Repeated out-of-fold analysis showed that prediction "
        "stability did not imply prediction correctness."
    )

    lines.append(
        "Minority-class error patterns remained persistent."
    )

    lines.append("")

    # ========================================================
    # PERMUTATION TEST
    # ========================================================

    lines.append(
        "TARGET-PERMUTATION SANITY TEST"
    )
    lines.append("-" * 80)

    lines.append(
        "Number of permutations: 100"
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
        "  Randomized runs >= real score: 95/100"
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
        "  Randomized runs >= real score: 78/100"
    )

    lines.append("")

    lines.append(
        "The observed model scores were not unusually high "
        "relative to randomized-target performance."
    )

    lines.append("")

    # ========================================================
    # SCIENTIFIC INTERPRETATION
    # ========================================================

    lines.append("=" * 80)
    lines.append(
        "SCIENTIFIC INTERPRETATION"
    )
    lines.append("=" * 80)

    lines.append(
        "The completed DeepGene workflow establishes a "
        "reproducible V1 variant dataset and a separately "
        "defined independent functional-effect ML dataset."
    )

    lines.append(
        "The exploratory supervised analysis did not provide "
        "sufficient evidence to select a reliable three-class "
        "functional-effect predictive model."
    )

    lines.append(
        "This is an insufficient/negative predictive finding "
        "for the current dataset and feature representation."
    )

    lines.append(
        "It does not demonstrate that the biological features "
        "have no relationship with functional effect."
    )

    lines.append(
        "It also does not establish clinical pathogenicity, "
        "benignity, disease risk, or clinical utility."
    )

    lines.append("")

    # ========================================================
    # FINAL MODEL STATUS
    # ========================================================

    lines.append(
        "FINAL MODEL STATUS"
    )
    lines.append("-" * 80)

    lines.append(
        "FINAL MODEL SELECTED: NO"
    )

    lines.append(
        "ML BRANCH: COMPLETE"
    )

    lines.append(
        "MODEL DEPLOYMENT: NOT JUSTIFIED"
    )

    lines.append(
        "CLINICAL USE: NOT ESTABLISHED"
    )

    lines.append("")

    # ========================================================
    # LIMITATIONS
    # ========================================================

    lines.append(
        "FINAL PROJECT LIMITATIONS"
    )
    lines.append("-" * 80)

    limitations = [
        "The primary supervised ML dataset contains only 51 variants.",
        "The Gain class contains only 4 variants.",
        "The functional dataset originates from a limited experimental source.",
        "No external independent functional benchmark was evaluated.",
        "The current biological feature representation is limited.",
        "Repeated cross-validation remains unstable.",
        "Feature importance is not stable.",
        "No final predictive model demonstrated robust predictive value.",
    ]

    for i, limitation in enumerate(
        limitations,
        start=1,
    ):

        lines.append(
            f"{i}. {limitation}"
        )

    lines.append("")

    # ========================================================
    # CLAIMS NOT SUPPORTED
    # ========================================================

    lines.append(
        "CLAIMS NOT SUPPORTED BY THE COMPLETED ANALYSIS"
    )
    lines.append("-" * 80)

    unsupported = [
        "Clinical diagnostic capability.",
        "Clinical pathogenicity or benignity classification.",
        "Disease-risk prediction.",
        "Clinical deployment or clinical utility.",
        "Causal biological conclusions from model feature importance.",
        "A conclusion that the tested biological features are unrelated to functional effect.",
    ]

    for i, statement in enumerate(
        unsupported,
        start=1,
    ):

        lines.append(
            f"{i}. {statement}"
        )

    lines.append("")

    # ========================================================
    # PROJECT ARTIFACTS
    # ========================================================

    lines.append(
        "PROJECT ARTIFACT CHECK"
    )
    lines.append("-" * 80)

    lines.append(
        f"{V1_DATASET}: "
        f"{status(V1_DATASET)}"
    )

    lines.append(
        f"{ML_DATASET}: "
        f"{status(ML_DATASET)}"
    )

    for report, report_status_value in report_status.items():

        lines.append(
            f"{ANALYSIS_DIR / report}: "
            f"{report_status_value}"
        )

    lines.append("")

    # ========================================================
    # FUTURE WORK
    # ========================================================

    lines.append(
        "FUTURE WORK BOUNDARY"
    )
    lines.append("-" * 80)

    lines.append(
        "The completed ML branch should not be extended by "
        "arbitrary model tuning or forced model selection."
    )

    lines.append(
        "If the project continues, scientifically justified "
        "future work may evaluate additional independent "
        "functional measurements or additional biologically "
        "justified features."
    )

    lines.append(
        "Any new supervised analysis should preserve the "
        "target-independence and leakage controls established "
        "in the current workflow."
    )

    lines.append(
        "Any future model should be treated as a new experimental "
        "analysis and independently validated."
    )

    lines.append("")

    # ========================================================
    # FINAL HANDOFF
    # ========================================================

    lines.append("=" * 80)
    lines.append(
        "FINAL DEEPGENE PROJECT HANDOFF"
    )
    lines.append("=" * 80)

    lines.append(
        "V1 DATASET: "
        + (
            "AVAILABLE"
            if dataset_status
            else "INCOMPLETE"
        )
    )

    lines.append(
        "ML DATASET VALIDATION: "
        + (
            "PASS"
            if ml_validation
            else "FAIL"
        )
    )

    lines.append(
        "ML ARTIFACT CHECK: "
        + (
            "PASS"
            if all_reports_present
            else "INCOMPLETE"
        )
    )

    lines.append(
        "LEAKAGE CONTROL: PASS"
    )

    lines.append(
        "ML BRANCH: COMPLETE"
    )

    lines.append(
        "FINAL MODEL: NONE SELECTED"
    )

    lines.append(
        "SCIENTIFIC HANDOFF: "
        + (
            "COMPLETE"
            if (
                dataset_status
                and ml_validation
                and all_reports_present
            )
            else "INCOMPLETE"
        )
    )

    # ========================================================
    # FINAL STATUS
    # ========================================================

    final_pass = (
        dataset_status
        and ml_validation
        and all_reports_present
    )

    lines.append("")
    lines.append("=" * 80)
    lines.append(
        f"STATUS: "
        f"{'PASS' if final_pass else 'FAIL'}"
    )
    lines.append("=" * 80)

    # ========================================================
    # WRITE REPORT
    # ========================================================

    OUTPUT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    # ========================================================
    # CONSOLE OUTPUT
    # ========================================================

    print(
        "\n".join(lines)
    )


if __name__ == "__main__":
    main()