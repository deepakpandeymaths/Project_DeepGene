# DeepGene exploratory functional predictor

Build the research-only artifact from the independent functional dataset:

```powershell
python data/ml/35_build_exploratory_predictor.py
```

Predict a missense substitution:

```powershell
python data/ml/36_predict_functional_effect.py --variant p.Arg1234Gly
```

The output is a three-class estimate (`Loss`, `Mixed`, `Gain`) with class
probabilities, normalized entropy, and an abstention flag. The artifact is
explicitly exploratory because the 51-row Generation-1 evaluation did not
establish performance above the majority/permutation baselines. It is not a
pathogenicity classifier, disease-risk predictor, diagnostic tool, or treatment
recommendation system.
