from pathlib import Path


# ============================================================
# ML STEP 17 — MODEL COMPARISON
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

REPORT_FILE = (
    PROJECT_ROOT
    / "data"
    / "analysis_results"
    / "17_model_comparison.txt"
)


# ------------------------------------------------------------
# Results from completed ML steps
# ------------------------------------------------------------

models = [
    {
        "name": "Majority-class baseline",
        "balanced_accuracy": 0.3333,
        "macro_precision": 0.2288,
        "macro_recall": 0.3333,
        "macro_f1": 0.2713,
        "std_balanced_accuracy": None,
        "std_macro_f1": None,
    },
    {
        "name": "Logistic Regression",
        "balanced_accuracy": 0.3299,
        "macro_precision": 0.2771,
        "macro_recall": 0.3299,
        "macro_f1": 0.2968,
        "std_balanced_accuracy": 0.0736,
        "std_macro_f1": 0.0657,
    },
    {
        "name": "Random Forest",
        "balanced_accuracy": 0.3530,
        "macro_precision": 0.2796,
        "macro_recall": 0.3530,
        "macro_f1": 0.2838,
        "std_balanced_accuracy": 0.2099,
        "std_macro_f1": 0.1519,
    },
]


# ------------------------------------------------------------
# Calculate differences from majority baseline
# ------------------------------------------------------------

baseline = models[0]

for model in models:
    model["delta_balanced_accuracy"] = (
        model["balanced_accuracy"]
        - baseline["balanced_accuracy"]
    )

    model["delta_macro_f1"] = (
        model["macro_f1"]
        - baseline["macro_f1"]
    )


# ------------------------------------------------------------
# Write report
# ------------------------------------------------------------

lines = []

lines.append("PROJECT DEEPGENE — ML STEP 17")
lines.append("MODEL COMPARISON")
lines.append("=" * 60)
lines.append("")

lines.append("MODELS INCLUDED")
lines.append("-" * 60)
lines.append(
    "1. Majority-class baseline"
)
lines.append(
    "2. Logistic Regression"
)
lines.append(
    "3. Random Forest"
)
lines.append("")

lines.append("PERFORMANCE SUMMARY")
lines.append("-" * 60)

lines.append(
    f"{'Model':<28}"
    f"{'Bal.Acc.':>12}"
    f"{'Macro Prec.':>14}"
    f"{'Macro Recall':>14}"
    f"{'Macro F1':>12}"
)

lines.append("-" * 80)

for model in models:

    lines.append(
        f"{model['name']:<28}"
        f"{model['balanced_accuracy']:>12.4f}"
        f"{model['macro_precision']:>14.4f}"
        f"{model['macro_recall']:>14.4f}"
        f"{model['macro_f1']:>12.4f}"
    )

lines.append("")

lines.append("CROSS-VALIDATION VARIABILITY")
lines.append("-" * 60)

lines.append(
    "Logistic Regression:"
)
lines.append(
    f"  Balanced accuracy SD: "
    f"{models[1]['std_balanced_accuracy']:.4f}"
)
lines.append(
    f"  Macro F1 SD: "
    f"{models[1]['std_macro_f1']:.4f}"
)
lines.append("")

lines.append(
    "Random Forest:"
)
lines.append(
    f"  Balanced accuracy SD: "
    f"{models[2]['std_balanced_accuracy']:.4f}"
)
lines.append(
    f"  Macro F1 SD: "
    f"{models[2]['std_macro_f1']:.4f}"
)
lines.append("")

lines.append("CHANGE FROM MAJORITY BASELINE")
lines.append("-" * 60)

for model in models[1:]:

    lines.append(
        f"{model['name']}:"
    )

    lines.append(
        f"  Balanced accuracy change: "
        f"{model['delta_balanced_accuracy']:+.4f}"
    )

    lines.append(
        f"  Macro F1 change: "
        f"{model['delta_macro_f1']:+.4f}"
    )

    lines.append("")


lines.append("INTERPRETATION")
lines.append("-" * 60)

lines.append(
    "The majority-class baseline provides the reference "
    "performance without biological features."
)

lines.append(
    "Logistic Regression has balanced accuracy of 0.3299 "
    "and macro F1 of 0.2968."
)

lines.append(
    "Random Forest has balanced accuracy of 0.3530 "
    "and macro F1 of 0.2838."
)

lines.append(
    "The Random Forest balanced-accuracy variability is "
    "substantial (SD=0.2099)."
)

lines.append(
    "The small number of Gain variants limits the stability "
    "and precision of model comparison."
)

lines.append(
    "These results are exploratory and do not establish "
    "clinical validity or clinical predictive performance."
)

lines.append(
    "No model is selected as a final model in this step."
)

lines.append("")

lines.append("STATUS: PASS")


REPORT_FILE.parent.mkdir(
    parents=True,
    exist_ok=True,
)

REPORT_FILE.write_text(
    "\n".join(lines),
    encoding="utf-8",
)


print("\n".join(lines))