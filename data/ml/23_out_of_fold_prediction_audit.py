from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
OUTPUT_FILE = Path("data/analysis_results/23_out_of_fold_prediction_audit.txt")

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

SEEDS = [42, 123, 2024, 7, 99]
N_SPLITS = 4

CLASSES = ["Loss", "Mixed", "Gain"]


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


def main():
    df = pd.read_csv(INPUT_FILE)

    required = FEATURES + [TARGET]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if df[required].isnull().any().any():
        raise ValueError("Missing values detected.")

    X = df[FEATURES]
    y = df[TARGET]

    if sorted(y.unique()) != sorted(CLASSES):
        raise ValueError(
            f"Unexpected target classes: {sorted(y.unique())}"
        )

    # --------------------------------------------------------
    # Collect repeated out-of-fold predictions
    # --------------------------------------------------------
    prediction_records = []

    for seed in SEEDS:
        cv = StratifiedKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=seed,
        )

        for fold, (train_idx, valid_idx) in enumerate(
            cv.split(X, y),
            start=1,
        ):
            X_train = X.iloc[train_idx]
            X_valid = X.iloc[valid_idx]
            y_train = y.iloc[train_idx]

            models = make_models()

            for model_name, model in models.items():
                model.fit(X_train, y_train)

                predictions = model.predict(X_valid)

                for idx, prediction in zip(
                    valid_idx,
                    predictions,
                ):
                    prediction_records.append(
                        {
                            "row_index": int(idx),
                            "seed": seed,
                            "fold": fold,
                            "model": model_name,
                            "actual": y.iloc[idx],
                            "predicted": prediction,
                        }
                    )

    pred_df = pd.DataFrame(prediction_records)

    # --------------------------------------------------------
    # Per-model prediction stability
    # --------------------------------------------------------
    lines = []

    lines.append(
        "DEEPGENE ML STEP 23 — OUT-OF-FOLD PREDICTION AUDIT"
    )
    lines.append("=" * 72)
    lines.append("")
    lines.append(
        "Purpose:"
    )
    lines.append(
        "Audit repeated out-of-fold predictions at the individual "
        "variant level to determine how consistently each model "
        "assigns functional-effect classes."
    )
    lines.append(
        "This is a descriptive stability analysis and does not "
        "establish clinical validity."
    )
    lines.append("")

    lines.append("INPUT")
    lines.append("-" * 72)
    lines.append(f"File: {INPUT_FILE}")
    lines.append(f"Rows: {len(df)}")
    lines.append(f"Features: {len(FEATURES)}")
    lines.append(f"Repeated CV seeds: {SEEDS}")
    lines.append(f"Folds per seed: {N_SPLITS}")
    lines.append("")

    # --------------------------------------------------------
    # Overall prediction distributions
    # --------------------------------------------------------
    for model_name in ["Logistic Regression", "Random Forest"]:
        model_df = pred_df[pred_df["model"] == model_name]

        lines.append(model_name.upper())
        lines.append("-" * 72)

        actual_counts = model_df["actual"].value_counts()

        lines.append("Actual class counts across validation predictions:")
        for cls in CLASSES:
            lines.append(
                f"{cls}: {actual_counts.get(cls, 0)}"
            )

        predicted_counts = model_df["predicted"].value_counts()

        lines.append("")
        lines.append("Predicted class counts:")
        for cls in CLASSES:
            lines.append(
                f"{cls}: {predicted_counts.get(cls, 0)}"
            )

        lines.append("")

        cm = confusion_matrix(
            model_df["actual"],
            model_df["predicted"],
            labels=CLASSES,
        )

        lines.append("Aggregate repeated-CV confusion matrix:")
        lines.append("Rows = actual; columns = predicted")
        lines.append("             " + "  ".join(f"{c:>7}" for c in CLASSES))

        for cls, row in zip(CLASSES, cm):
            lines.append(
                f"{cls:>12} "
                + "  ".join(f"{int(v):7d}" for v in row)
            )

        lines.append("")

    # --------------------------------------------------------
    # Variant-level prediction stability
    # --------------------------------------------------------
    lines.append("VARIANT-LEVEL PREDICTION STABILITY")
    lines.append("-" * 72)

    for model_name in ["Logistic Regression", "Random Forest"]:
        model_df = pred_df[pred_df["model"] == model_name]

        lines.append("")
        lines.append(model_name)

        grouped = (
            model_df
            .groupby(["row_index", "actual"])["predicted"]
            .agg(list)
            .reset_index()
        )

        stability_values = []

        for _, row in grouped.iterrows():
            predictions = row["predicted"]

            counts = pd.Series(predictions).value_counts()

            dominant_class = counts.index[0]
            dominant_count = counts.iloc[0]

            stability = dominant_count / len(predictions)

            stability_values.append(stability)

            lines.append(
                f"row={int(row['row_index'])}, "
                f"actual={row['actual']}, "
                f"dominant_prediction={dominant_class}, "
                f"dominant_fraction={stability:.3f}, "
                f"predictions="
                f"{dict(counts)}"
            )

        stability_values = np.asarray(stability_values)

        lines.append("")
        lines.append(
            f"Mean variant prediction stability: "
            f"{np.mean(stability_values):.4f}"
        )
        lines.append(
            f"Median variant prediction stability: "
            f"{np.median(stability_values):.4f}"
        )
        lines.append(
            f"Variants with 100% identical predictions: "
            f"{np.sum(stability_values == 1.0)}/{len(stability_values)}"
        )
        lines.append(
            f"Variants with prediction stability < 0.75: "
            f"{np.sum(stability_values < 0.75)}/{len(stability_values)}"
        )

    lines.append("")

    # --------------------------------------------------------
    # Actual-class stability
    # --------------------------------------------------------
    lines.append("STABILITY BY ACTUAL FUNCTIONAL CLASS")
    lines.append("-" * 72)

    for model_name in ["Logistic Regression", "Random Forest"]:
        model_df = pred_df[pred_df["model"] == model_name]

        lines.append("")
        lines.append(model_name)

        for cls in CLASSES:
            cls_df = model_df[model_df["actual"] == cls]

            if cls_df.empty:
                continue

            stability_values = []

            for row_index in sorted(cls_df["row_index"].unique()):
                predictions = cls_df[
                    cls_df["row_index"] == row_index
                ]["predicted"]

                counts = predictions.value_counts()
                stability_values.append(
                    counts.iloc[0] / len(predictions)
                )

            lines.append(
                f"{cls}: "
                f"variants={len(stability_values)}, "
                f"mean_stability="
                f"{np.mean(stability_values):.4f}, "
                f"median_stability="
                f"{np.median(stability_values):.4f}"
            )

    lines.append("")

    # --------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------
    lines.append("INTERPRETATION")
    lines.append("-" * 72)
    lines.append(
        "Repeated out-of-fold predictions are used here to assess "
        "prediction stability across different stratified validation "
        "assignments."
    )
    lines.append(
        "Low variant-level stability indicates that predictions depend "
        "strongly on which variants are included in the training fold."
    )
    lines.append(
        "Because the dataset contains only 4 Gain variants, prediction "
        "behavior for Gain is expected to be particularly sensitive "
        "to fold composition."
    )
    lines.append(
        "This analysis does not establish that a model is clinically "
        "useful and does not establish pathogenicity or benignity."
    )
    lines.append(
        "No final model is selected."
    )
    lines.append("")

    lines.append("STATUS")
    lines.append("-" * 72)
    lines.append("PASS")
    lines.append(
        "Out-of-fold prediction stability audit completed successfully."
    )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("\n".join(lines))


if __name__ == "__main__":
    main()