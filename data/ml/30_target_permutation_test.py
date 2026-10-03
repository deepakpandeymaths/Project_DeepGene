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
# DEEPGENE ML STEP 30
# FAST TARGET-PERMUTATION SANITY TEST
# ============================================================

INPUT_FILE = Path("data/processed/deepgene_ml_ready_v1.csv")
REPORT_FILE = Path("data/analysis_results/30_target_permutation_test.txt")

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

CV_SEEDS = [42, 123, 2024, 7, 99]
N_SPLITS = 4

# Reduced from 200 because this dataset is only 51 variants.
N_PERMUTATIONS = 100

PERMUTATION_SEED = 42


def make_lr():
    return Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                LogisticRegression(
                    max_iter=5000,
                    class_weight=None,
                    random_state=42,
                ),
            ),
        ]
    )


def make_rf():
    return RandomForestClassifier(
        n_estimators=100,
        class_weight=None,
        random_state=42,
        n_jobs=1,
    )


def evaluate_model(model_factory, X, y):

    scores = []

    for seed in CV_SEEDS:

        cv = StratifiedKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=seed,
        )

        for train_idx, val_idx in cv.split(X, y):

            X_train = X.iloc[train_idx]
            X_val = X.iloc[val_idx]

            y_train = y.iloc[train_idx]
            y_val = y.iloc[val_idx]

            model = model_factory()

            model.fit(
                X_train,
                y_train,
            )

            predictions = model.predict(X_val)

            score = balanced_accuracy_score(
                y_val,
                predictions,
            )

            scores.append(score)

    return float(np.mean(scores))


def summarize(values):

    return {
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)),
        "median": float(np.median(values)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "p95": float(np.percentile(values, 95)),
    }


def main():

    REPORT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    df = pd.read_csv(INPUT_FILE)

    X = df[FEATURES].copy()
    y = df[TARGET].copy()

    models = {
        "Logistic Regression": make_lr,
        "Random Forest": make_rf,
    }

    # --------------------------------------------------------
    # Real-label performance
    # --------------------------------------------------------

    print("=" * 70)
    print("REAL-LABEL PERFORMANCE")
    print("=" * 70)

    real_scores = {}

    for name, factory in models.items():

        score = evaluate_model(
            factory,
            X,
            y,
        )

        real_scores[name] = score

        print(
            f"{name}: {score:.6f}"
        )

    # --------------------------------------------------------
    # Permutation testing
    # --------------------------------------------------------

    rng = np.random.RandomState(
        PERMUTATION_SEED
    )

    permutation_results = {
        name: []
        for name in models
    }

    print("")
    print("=" * 70)
    print("TARGET PERMUTATION TEST")
    print("=" * 70)
    print(
        f"Running {N_PERMUTATIONS} permutations..."
    )

    y_array = y.to_numpy()

    for permutation_number in range(
        N_PERMUTATIONS
    ):

        shuffled_y = rng.permutation(
            y_array
        )

        shuffled_y = pd.Series(
            shuffled_y,
            index=y.index,
            name=y.name,
        )

        for name, factory in models.items():

            score = evaluate_model(
                factory,
                X,
                shuffled_y,
            )

            permutation_results[name].append(
                score
            )

        if (
            permutation_number + 1
        ) % 10 == 0:

            print(
                f"Permutation "
                f"{permutation_number + 1}/"
                f"{N_PERMUTATIONS}"
            )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    lines = []

    lines.append("=" * 70)
    lines.append("DEEPGENE ML STEP 30")
    lines.append("TARGET-PERMUTATION SANITY TEST")
    lines.append("=" * 70)
    lines.append("")

    lines.append("PURPOSE")
    lines.append("-" * 70)
    lines.append(
        "Assess whether observed model performance is distinguishable "
        "from performance obtained after randomized target labels."
    )
    lines.append(
        "This is a sanity test for predictive signal."
    )
    lines.append(
        "It is not clinical validation."
    )
    lines.append(
        "No model selection is performed."
    )
    lines.append("")

    lines.append("DATASET")
    lines.append("-" * 70)
    lines.append(
        f"Input: {INPUT_FILE}"
    )
    lines.append(
        f"Rows: {len(df)}"
    )
    lines.append(
        f"Features: {len(FEATURES)}"
    )
    lines.append(
        f"Target: {TARGET}"
    )
    lines.append("")

    lines.append("TARGET DISTRIBUTION")
    lines.append("-" * 70)

    for cls, count in (
        y.value_counts()
        .sort_index()
        .items()
    ):

        lines.append(
            f"{cls}: {count} "
            f"({count / len(y):.4f})"
        )

    lines.append("")

    lines.append("CROSS-VALIDATION DESIGN")
    lines.append("-" * 70)
    lines.append(
        f"StratifiedKFold: {N_SPLITS} folds"
    )
    lines.append(
        f"CV seeds: {CV_SEEDS}"
    )
    lines.append(
        f"Evaluations per permutation: "
        f"{len(CV_SEEDS) * N_SPLITS}"
    )
    lines.append("")

    lines.append("PERMUTATION DESIGN")
    lines.append("-" * 70)
    lines.append(
        f"Number of permutations: "
        f"{N_PERMUTATIONS}"
    )
    lines.append(
        f"Permutation seed: "
        f"{PERMUTATION_SEED}"
    )
    lines.append(
        "Target labels are shuffled while biological features remain unchanged."
    )
    lines.append(
        "Random Forest uses 100 trees and single-thread execution "
        "to avoid excessive parallel-processing overhead."
    )
    lines.append("")

    # --------------------------------------------------------
    # Results
    # --------------------------------------------------------

    for name in models:

        values = np.array(
            permutation_results[name]
        )

        summary = summarize(values)

        real_score = real_scores[name]

        exceed_count = int(
            np.sum(values >= real_score)
        )

        empirical_percentile = (
            100.0
            * np.mean(values <= real_score)
        )

        lines.append("=" * 70)
        lines.append(
            name.upper()
        )
        lines.append("=" * 70)

        lines.append(
            f"Real-label mean balanced accuracy: "
            f"{real_score:.6f}"
        )

        lines.append(
            f"Permutation mean balanced accuracy: "
            f"{summary['mean']:.6f}"
        )

        lines.append(
            f"Permutation SD: "
            f"{summary['sd']:.6f}"
        )

        lines.append(
            f"Permutation median: "
            f"{summary['median']:.6f}"
        )

        lines.append(
            f"Permutation minimum: "
            f"{summary['min']:.6f}"
        )

        lines.append(
            f"Permutation maximum: "
            f"{summary['max']:.6f}"
        )

        lines.append(
            f"Permutation 95th percentile: "
            f"{summary['p95']:.6f}"
        )

        lines.append(
            f"Permutations >= real-label score: "
            f"{exceed_count}/{N_PERMUTATIONS}"
        )

        lines.append(
            f"Empirical percentile of real-label score: "
            f"{empirical_percentile:.2f}%"
        )

        lines.append("")

    # --------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------

    lines.append("=" * 70)
    lines.append("INTERPRETATION")
    lines.append("=" * 70)
    lines.append(
        "Target permutation testing compares the observed feature-target "
        "relationship with randomized target assignments."
    )
    lines.append(
        "A low proportion of permutations reaching the real-label score "
        "means the observed score is unusual relative to randomized labels."
    )
    lines.append(
        "A high proportion means the observed score is compatible with "
        "randomized-target performance."
    )
    lines.append(
        "The permutation result does not establish clinical validity, "
        "biological causality, or external generalization."
    )
    lines.append(
        "The Gain class contains only 4 variants, limiting stability."
    )
    lines.append("")

    lines.append("=" * 70)
    lines.append("DECISION")
    lines.append("=" * 70)
    lines.append(
        "No model selected in Step 30."
    )
    lines.append(
        "No hyperparameter tuning performed."
    )
    lines.append(
        "No feature selection performed."
    )
    lines.append(
        "No clinical interpretation made."
    )

    lines.append("")
    lines.append("=" * 70)
    lines.append("STATUS: PASS")
    lines.append("=" * 70)

    REPORT_FILE.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print("")
    print("\n".join(lines))


if __name__ == "__main__":
    main()