from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# ============================================================
# ML STEP 22
# Controlled uncertainty analysis using repeated stratified CV
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
OUTPUT_FILE = Path("data/analysis_results/22_model_uncertainty_analysis.txt")

RANDOM_SEEDS = [42, 123, 2024, 7, 99]
N_SPLITS = 4

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


def make_models():
    return {
        "Logistic Regression": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "model",
                    LogisticRegression(
                        max_iter=5000,
                        random_state=42,
                    ),
                ),
            ]
        ),
        "Random Forest": RandomForestClassifier(
            n_estimators=300,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
    }


def percentile_interval(values, lower=2.5, upper=97.5):
    values = np.asarray(values, dtype=float)
    return (
        float(np.percentile(values, lower)),
        float(np.percentile(values, upper)),
    )


def summarize(values):
    values = np.asarray(values, dtype=float)
    low, high = percentile_interval(values)

    return {
        "n": len(values),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)),
        "median": float(np.median(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "ci_low": low,
        "ci_high": high,
    }


def main():
    df = pd.read_csv(INPUT_FILE)

    # --------------------------------------------------------
    # Basic validation
    # --------------------------------------------------------
    required_columns = [TARGET] + FEATURES
    missing = [c for c in required_columns if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if df[required_columns].isnull().any().any():
        raise ValueError("Missing values detected in target/features.")

    X = df[FEATURES].copy()
    y = df[TARGET].copy()

    class_counts = y.value_counts()

    if class_counts.min() < N_SPLITS:
        raise ValueError(
            f"Smallest class has {class_counts.min()} samples, "
            f"but {N_SPLITS} folds are required."
        )

    # --------------------------------------------------------
    # Majority-class baseline
    # --------------------------------------------------------
    majority_class = class_counts.idxmax()

    baseline_predictions = np.full(
        shape=len(y),
        fill_value=majority_class,
        dtype=object,
    )

    baseline_ba = balanced_accuracy_score(y, baseline_predictions)

    # --------------------------------------------------------
    # Repeated stratified CV
    # --------------------------------------------------------
    results = {
        "Logistic Regression": [],
        "Random Forest": [],
    }

    model_fold_records = {
        "Logistic Regression": [],
        "Random Forest": [],
    }

    for seed in RANDOM_SEEDS:
        cv = StratifiedKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=seed,
        )

        for fold_number, (train_idx, valid_idx) in enumerate(
            cv.split(X, y),
            start=1,
        ):
            X_train = X.iloc[train_idx]
            X_valid = X.iloc[valid_idx]
            y_train = y.iloc[train_idx]
            y_valid = y.iloc[valid_idx]

            models = make_models()

            for model_name, model in models.items():
                model.fit(X_train, y_train)

                predictions = model.predict(X_valid)

                ba = balanced_accuracy_score(
                    y_valid,
                    predictions,
                )

                results[model_name].append(float(ba))

                model_fold_records[model_name].append(
                    {
                        "seed": seed,
                        "fold": fold_number,
                        "training_size": len(train_idx),
                        "validation_size": len(valid_idx),
                        "balanced_accuracy": float(ba),
                    }
                )

    # --------------------------------------------------------
    # Model-vs-baseline differences
    #
    # The baseline is fixed at 0.3333. We do NOT claim
    # statistical significance. This is descriptive uncertainty.
    # --------------------------------------------------------
    difference_results = {}

    for model_name, values in results.items():
        values = np.asarray(values)

        differences = values - baseline_ba

        difference_results[model_name] = differences

    # --------------------------------------------------------
    # Write report
    # --------------------------------------------------------
    lines = []

    lines.append("DEEPGENE ML STEP 22 — CONTROLLED UNCERTAINTY ANALYSIS")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Purpose:")
    lines.append(
        "Quantify the observed variability of Logistic Regression and "
        "Random Forest using the same repeated stratified cross-validation "
        "design used in Steps 18–21."
    )
    lines.append(
        "This is a descriptive uncertainty analysis, not a claim of "
        "statistical significance or clinical validity."
    )
    lines.append("")

    lines.append("INPUT")
    lines.append("-" * 72)
    lines.append(f"File: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Features: {len(FEATURES)}")
    lines.append(f"Target: {TARGET}")
    lines.append("")

    lines.append("TARGET DISTRIBUTION")
    lines.append("-" * 72)
    for class_name, count in class_counts.items():
        lines.append(
            f"{class_name}: {count} "
            f"({count / len(y):.4f})"
        )
    lines.append("")

    lines.append("BASELINE")
    lines.append("-" * 72)
    lines.append(f"Majority class: {majority_class}")
    lines.append(f"Balanced accuracy: {baseline_ba:.6f}")
    lines.append("")

    lines.append("CROSS-VALIDATION DESIGN")
    lines.append("-" * 72)
    lines.append(f"Seeds: {RANDOM_SEEDS}")
    lines.append(f"Folds per seed: {N_SPLITS}")
    lines.append(
        f"Total validation evaluations per model: "
        f"{len(RANDOM_SEEDS) * N_SPLITS}"
    )
    lines.append(
        "Metric: balanced accuracy"
    )
    lines.append(
        "No hyperparameter tuning performed."
    )
    lines.append("")

    # --------------------------------------------------------
    # Summary statistics
    # --------------------------------------------------------
    for model_name, values in results.items():
        summary = summarize(values)

        lines.append(model_name.upper())
        lines.append("-" * 72)
        lines.append(f"Evaluations: {summary['n']}")
        lines.append(f"Mean: {summary['mean']:.6f}")
        lines.append(f"SD: {summary['sd']:.6f}")
        lines.append(f"Median: {summary['median']:.6f}")
        lines.append(f"Minimum: {summary['min']:.6f}")
        lines.append(f"Maximum: {summary['max']:.6f}")
        lines.append(
            "Empirical 95% interval across repeated CV evaluations: "
            f"{summary['ci_low']:.6f} to {summary['ci_high']:.6f}"
        )

        diff_summary = summarize(difference_results[model_name])

        lines.append(
            f"Mean difference vs majority baseline: "
            f"{diff_summary['mean']:.6f}"
        )
        lines.append(
            "Empirical 95% interval for fold-level difference vs baseline: "
            f"{diff_summary['ci_low']:.6f} to "
            f"{diff_summary['ci_high']:.6f}"
        )

        lines.append("")

    # --------------------------------------------------------
    # Fold-level observations
    # --------------------------------------------------------
    lines.append("FOLD-LEVEL RESULTS")
    lines.append("-" * 72)

    for model_name, records in model_fold_records.items():
        lines.append("")
        lines.append(model_name)

        for record in records:
            lines.append(
                f"seed={record['seed']}, "
                f"fold={record['fold']}, "
                f"train={record['training_size']}, "
                f"valid={record['validation_size']}, "
                f"balanced_accuracy="
                f"{record['balanced_accuracy']:.6f}"
            )

    lines.append("")

    # --------------------------------------------------------
    # Explicit interpretation without selecting a model
    # --------------------------------------------------------
    lines.append("INTERPRETATION")
    lines.append("-" * 72)

    for model_name, values in results.items():
        values = np.asarray(values)

        below_or_equal = int(np.sum(values <= baseline_ba))
        above = int(np.sum(values > baseline_ba))

        lines.append(
            f"{model_name}: {above}/{len(values)} repeated validation "
            f"evaluations exceeded the majority-class balanced accuracy "
            f"of {baseline_ba:.6f}; "
            f"{below_or_equal}/{len(values)} were at or below it."
        )

    lines.append("")
    lines.append(
        "The repeated-validation distributions should be interpreted "
        "as evidence of substantial uncertainty in model performance "
        "on this 51-variant dataset."
    )
    lines.append(
        "The empirical intervals describe variation across the chosen "
        "validation evaluations; they are not formal confidence intervals "
        "for a population-level model performance parameter."
    )
    lines.append(
        "Because the Gain class contains only 4 variants, estimates "
        "involving that class remain particularly unstable."
    )
    lines.append(
        "No final model is selected from this analysis."
    )
    lines.append(
        "No hyperparameter tuning is performed."
    )
    lines.append(
        "No clinical interpretation or pathogenicity/benignity conclusion "
        "is made."
    )
    lines.append("")

    lines.append("STATUS")
    lines.append("-" * 72)
    lines.append("PASS")
    lines.append(
        "Controlled uncertainty analysis completed successfully."
    )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text("\n".join(lines), encoding="utf-8")

    print("\n".join(lines))


if __name__ == "__main__":
    main()