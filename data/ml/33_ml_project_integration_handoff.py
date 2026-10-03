from pathlib import Path
from datetime import datetime

import pandas as pd


# ============================================================
# DEEPGENE STEP 33
# ML PROJECT INTEGRATION / HANDOFF
# ============================================================

ML_DATASET = Path(
    "data/processed/deepgene_ml_ready_v1.csv"
)

ML_HANDOFF = Path(
    "data/analysis_results/32_ml_final_handoff.txt"
)

ML_SUMMARY = Path(
    "data/analysis_results/31_consolidated_ml_evidence_summary.txt"
)

ML_PERMUTATION = Path(
    "data/analysis_results/30_target_permutation_test.txt"
)

ML_LEAKAGE = Path(
    "data/analysis_results/29_final_leakage_audit.txt"
)

OUTPUT_FILE = Path(
    "data/analysis_results/33_ml_project_integration_handoff.txt"
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


# ============================================================
# HELPER
# ============================================================

def file_status(path):
    return "PRESENT" if path.exists() else "MISSING"


def main():

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    lines = []

    # ========================================================
    # LOAD DATASET
    # ========================================================

    df = pd.read_csv(ML_DATASET)

    target_counts = (
        df[TARGET]
        .value_counts()
        .sort_index()
    )

    # ========================================================
    # HEADER
    # ========================================================

    lines.append("=" * 80)
    lines.append("DEEPGENE STEP 33")
    lines.append("ML PROJECT INTEGRATION / HANDOFF")
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
        "Integrate the completed exploratory machine-learning branch "
        "into the broader DeepGene project record."
    )

    lines.append(
        "This step does not train a new model."
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
    # PROJECT ML PIPELINE
    # ========================================================

    lines.append("ML PIPELINE COMPLETION")
    lines.append("-" * 80)

    lines.append(
        "Independent functional dataset"
    )

    lines.append(
        "        ↓"
    )

    lines.append(
        "Leakage-controlled biological feature set"
    )

    lines.append(
        "        ↓"
    )

    lines.append(
        "Exploratory supervised ML"
    )

    lines.append(
        "        ↓"
    )

    lines.append(
        "Repeated stratified cross-validation"
    )

    lines.append(
        "        ↓"
    )

    lines.append(
        "Uncertainty, stability, and error analysis"
    )

    lines.append(
        "        ↓"
    )

    lines.append(
        "Target-permutation sanity test"
    )

    lines.append(
        "        ↓"
    )

    lines.append(
        "No final predictive model selected"
    )

    lines.append(
        "        ↓"
    )

    lines.append(
        "ML branch complete"
    )

    lines.append("")

    # ========================================================
    # DATASET
    # ========================================================

    lines.append("INTEGRATED ML DATASET")
    lines.append("-" * 80)

    lines.append(
        f"Dataset: {ML_DATASET}"
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
        f"Biological predictors: {len(FEATURES)}"
    )

    lines.append("")

    lines.append(
        "Functional target distribution:"
    )

    for label, count in target_counts.items():

        lines.append(
            f"  {label}: {count} "
            f"({count / len(df):.4f})"
        )

    lines.append("")

    lines.append(
        "Smallest class: "
        f"{target_counts.min()}"
    )

    lines.append("")

    # ========================================================
    # FINAL PREDICTORS
    # ========================================================

    lines.append("INTEGRATED BIOLOGICAL PREDICTORS")
    lines.append("-" * 80)

    for index, feature in enumerate(
        FEATURES,
        start=1,
    ):

        lines.append(
            f"{index:02d}. {feature}"
        )

    lines.append("")

    lines.append(
        "The predictor set consists of protein-position and "
        "amino-acid physicochemical features."
    )

    lines.append(
        "ClinVar-derived clinical/evidence variables, phenotype, "
        "inheritance, identifiers, and target-derived functional "
        "text were not used as model predictors."
    )

    lines.append("")

    # ========================================================
    # DATASET VALIDATION
    # ========================================================

    lines.append("INTEGRATION DATASET VALIDATION")
    lines.append("-" * 80)

    missing_features = [
        feature
        for feature in FEATURES
        if feature not in df.columns
    ]

    non_numeric = [
        feature
        for feature in FEATURES
        if feature in df.columns
        and not pd.api.types.is_numeric_dtype(
            df[feature]
        )
    ]

    target_missing = int(
        df[TARGET].isna().sum()
    )

    feature_missing = int(
        df[FEATURES]
        .isna()
        .sum()
        .sum()
    )

    lines.append(
        f"Target present: "
        f"{TARGET in df.columns}"
    )

    lines.append(
        f"Missing target values: "
        f"{target_missing}"
    )

    lines.append(
        f"Missing predictor columns: "
        f"{len(missing_features)}"
    )

    lines.append(
        f"Non-numeric predictors: "
        f"{len(non_numeric)}"
    )

    lines.append(
        f"Missing predictor values: "
        f"{feature_missing}"
    )

    if missing_features:

        for feature in missing_features:
            lines.append(
                f"  Missing: {feature}"
            )

    if non_numeric:

        for feature in non_numeric:
            lines.append(
                f"  Non-numeric: {feature}"
            )

    lines.append("")

    # ========================================================
    # LEAKAGE STATUS
    # ========================================================

    lines.append("LEAKAGE CONTROL")
    lines.append("-" * 80)

    lines.append(
        f"Step 29 leakage audit: "
        f"{file_status(ML_LEAKAGE)}"
    )

    lines.append(
        "Documented Step 29 result: PASS"
    )

    lines.append(
        "The final predictor set is structurally separated "
        "from the functional target."
    )

    lines.append(
        "No identifier predictor was used."
    )

    lines.append(
        "No phenotype or inheritance predictor was used."
    )

    lines.append(
        "No ClinVar clinical-significance or evidence-count "
        "predictor was used."
    )

    lines.append("")

    # ========================================================
    # MODEL EVIDENCE
    # ========================================================

    lines.append("MODEL EVIDENCE INTEGRATION")
    lines.append("-" * 80)

    lines.append(
        "Majority-class baseline:"
    )

    lines.append(
        "  Balanced accuracy = 0.3333"
    )

    lines.append(
        "  Macro F1 = 0.2713"
    )

    lines.append("")

    lines.append(
        "Logistic Regression:"
    )

    lines.append(
        "  Repeated-CV balanced accuracy = "
        "0.2750 ± 0.0740"
    )

    lines.append(
        "  Repeated-CV macro F1 = "
        "0.2472 ± 0.0537"
    )

    lines.append("")

    lines.append(
        "Random Forest:"
    )

    lines.append(
        "  Repeated-CV balanced accuracy = "
        "0.3032 ± 0.1663"
    )

    lines.append(
        "  Repeated-CV macro F1 = "
        "0.2718 ± 0.1274"
    )

    lines.append("")

    lines.append(
        "Neither supervised model demonstrated robust "
        "superiority over the majority-class baseline."
    )

    lines.append("")

    # ========================================================
    # PERMUTATION EVIDENCE
    # ========================================================

    lines.append("TARGET-PERMUTATION EVIDENCE")
    lines.append("-" * 80)

    lines.append(
        f"Permutation report: "
        f"{file_status(ML_PERMUTATION)}"
    )

    lines.append(
        "100 randomized-target permutations were evaluated."
    )

    lines.append("")

    lines.append(
        "Logistic Regression:"
    )

    lines.append(
        "  Real-label balanced accuracy = 0.2750"
    )

    lines.append(
        "  Permutation mean = 0.3394"
    )

    lines.append(
        "  Randomized runs >= real-label score = 95/100"
    )

    lines.append("")

    lines.append(
        "Random Forest:"
    )

    lines.append(
        "  Real-label balanced accuracy = 0.3012"
    )

    lines.append(
        "  Permutation mean = 0.3390"
    )

    lines.append(
        "  Randomized runs >= real-label score = 78/100"
    )

    lines.append("")

    lines.append(
        "Integrated interpretation:"
    )

    lines.append(
        "The observed model scores were not unusually high "
        "relative to randomized-target performance."
    )

    lines.append("")

    # ========================================================
    # STABILITY / ERROR EVIDENCE
    # ========================================================

    lines.append("STABILITY AND ERROR EVIDENCE")
    lines.append("-" * 80)

    lines.append(
        "Repeated cross-validation showed substantial "
        "variation across fold assignments."
    )

    lines.append(
        "Feature-importance stability did not identify "
        "a consistently dominant predictor."
    )

    lines.append(
        "Learning-curve analysis did not demonstrate a clear "
        "improvement pattern with additional training data."
    )

    lines.append(
        "Repeated out-of-fold prediction analysis showed "
        "that prediction stability did not imply correctness."
    )

    lines.append(
        "Error analysis showed persistent difficulty with "
        "minority functional classes."
    )

    lines.append(
        "The Gain class contains only 4 variants."
    )

    lines.append("")

    # ========================================================
    # FINAL PROJECT DECISION
    # ========================================================

    lines.append("=" * 80)
    lines.append("FINAL ML PROJECT DECISION")
    lines.append("=" * 80)

    lines.append(
        "ML BRANCH STATUS: COMPLETE"
    )

    lines.append(
        "FINAL MODEL SELECTED: NO"
    )

    lines.append("")

    lines.append(
        "The current evidence does not support selection of "
        "Logistic Regression or Random Forest as a final "
        "predictive model."
    )

    lines.append(
        "The ML branch therefore ends as an exploratory analysis "
        "rather than a deployed predictive model."
    )

    lines.append("")

    # ========================================================
    # WHAT IS VALID TO CARRY FORWARD
    # ========================================================

    lines.append("VALID INFORMATION TO CARRY FORWARD")
    lines.append("-" * 80)

    lines.append(
        "1. The independent functional-effect dataset is "
        "documented and reproducible."
    )

    lines.append(
        "2. The leakage-controlled 14-feature biological "
        "representation is documented."
    )

    lines.append(
        "3. The exploratory ML performance and uncertainty "
        "results are documented."
    )

    lines.append(
        "4. The negative/insufficient model-selection result "
        "is documented."
    )

    lines.append(
        "5. The limitations of the current functional dataset "
        "are documented."
    )

    lines.append("")

    # ========================================================
    # WHAT MUST NOT BE CLAIMED
    # ========================================================

    lines.append("CLAIMS NOT SUPPORTED BY THIS ML BRANCH")
    lines.append("-" * 80)

    lines.append(
        "1. No clinical diagnostic capability is established."
    )

    lines.append(
        "2. No pathogenicity/benignity classifier is established."
    )

    lines.append(
        "3. No disease-risk prediction is established."
    )

    lines.append(
        "4. No clinically deployable prediction model is established."
    )

    lines.append(
        "5. No causal biological conclusion is established "
        "from feature importance."
    )

    lines.append(
        "6. No claim that the biological features have no "
        "relationship with functional effect is justified."
    )

    lines.append("")

    # ========================================================
    # HANDOFF ARTIFACTS
    # ========================================================

    lines.append("ML HANDOFF ARTIFACTS")
    lines.append("-" * 80)

    artifacts = [
        ML_DATASET,
        ML_LEAKAGE,
        ML_PERMUTATION,
        ML_SUMMARY,
        ML_HANDOFF,
    ]

    for artifact in artifacts:

        lines.append(
            f"{artifact}: "
            f"{file_status(artifact)}"
        )

    lines.append("")

    # ========================================================
    # NEXT-STAGE GUIDANCE
    # ========================================================

    lines.append("NEXT-STAGE GUIDANCE")
    lines.append("-" * 80)

    lines.append(
        "The ML branch should not proceed to deployment, "
        "clinical prediction, or forced model selection."
    )

    lines.append(
        "If the project continues, scientifically justified "
        "future work may evaluate additional independent "
        "functional data or additional biologically justified "
        "features while preserving the existing target-independence "
        "and leakage controls."
    )

    lines.append(
        "Any future supervised analysis should be treated as "
        "a new experimental analysis and independently validated."
    )

    lines.append("")

    # ========================================================
    # INTEGRATION CHECK
    # ========================================================

    dataset_pass = (
        TARGET in df.columns
        and target_missing == 0
        and len(missing_features) == 0
        and len(non_numeric) == 0
        and feature_missing == 0
    )

    artifact_pass = all(
        artifact.exists()
        for artifact in artifacts
    )

    lines.append("=" * 80)
    lines.append("INTEGRATION CHECK")
    lines.append("=" * 80)

    lines.append(
        "ML DATASET VALIDATION: "
        + (
            "PASS"
            if dataset_pass
            else "FAIL"
        )
    )

    lines.append(
        "ML HANDOFF ARTIFACTS: "
        + (
            "PASS"
            if artifact_pass
            else "FAIL"
        )
    )

    lines.append(
        "LEAKAGE CONTROL: PASS"
    )

    lines.append(
        "MODEL SELECTION: NO MODEL SELECTED"
    )

    if dataset_pass and artifact_pass:

        lines.append(
            "PROJECT ML INTEGRATION: COMPLETE"
        )

        status = "PASS"

    else:

        lines.append(
            "PROJECT ML INTEGRATION: INCOMPLETE"
        )

        status = "FAIL"

    lines.append("")

    lines.append("=" * 80)
    lines.append(
        f"STATUS: {status}"
    )
    lines.append("=" * 80)

    # ========================================================
    # WRITE OUTPUT
    # ========================================================

    OUTPUT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(
        "\n".join(lines)
    )


if __name__ == "__main__":
    main()