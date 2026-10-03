# 🧬 DeepGene — SCN1A Variant Evidence & Prioritization

DeepGene is a research-oriented bioinformatics/data-science project for **systematic analysis and prioritization of genetic variants**, initially focused on **SCN1A** and SCN1A-related disorders.

The project starts with ClinVar evidence and is being developed as a transparent, reproducible pipeline before introducing machine-learning models.

> **Scope:** DeepGene V1 is an evidence representation and prioritization research project. It is **not a clinical diagnostic system** and its outputs should not be used as a medical diagnosis or treatment decision.

---

## 1. Project Vision

The long-term goal of DeepGene is to build an interpretable computational framework that can combine different types of variant evidence and eventually prioritize variants for further scientific investigation.

The planned progression is:

```text
Raw genetic variant data
        ↓
Database + data quality
        ↓
Evidence representation
        ↓
Transparent prioritization
        ↓
Independent benchmarking
        ↓
Machine learning
        ↓
Phenotype-aware models
        ↓
Advanced models
       ├── Gradient Boosted Decision Trees
       ├── Multilayer Perceptron
       └── Graph Neural Network
```

The central principle is:

> **Build a scientifically transparent baseline first, then add model complexity only when the data and evaluation design justify it.**

---

# 2. Initial Biological Focus

### Gene

**SCN1A — Sodium Voltage-Gated Channel Alpha Subunit 1**

### Initial disease context

SCN1A-related disorders, including disorders associated with developmental and epileptic encephalopathy and related epilepsy phenotypes.

The initial version focuses on **SCN1A variants available in ClinVar**.

---

# 3. V1 Data Strategy

DeepGene V1 intentionally starts with:

```text
ClinVar only
```

External resources are not being added at this stage.

Future versions may incorporate additional evidence such as:

- Population frequency
- Functional evidence
- Conservation
- Splice prediction
- Protein-impact prediction
- Protein structure
- Phenotype ontologies
- Independent literature evidence
- Patient-level phenotype/genotype information

These features are **not invented when unavailable**.

A critical distinction is maintained between:

```text
0  = evidence was assessed and absent

NA = evidence is currently unavailable
```

This prevents missing information from being incorrectly interpreted as negative evidence.

---

# 4. Project Development History

## Phase 1 — Data Acquisition

The project began by extracting SCN1A-specific records from the ClinVar `variant_summary` dataset.

The raw ClinVar dataset contained:

```text
10,696 SCN1A ClinVar records
```

These records became the raw evidence source for DeepGene V1.

Raw data location:

```text
data/raw/SCN1A_clinvar.tsv
```

---

# 5. Phase 2 — SQLite Database

A SQLite database was created to make the dataset queryable and reproducible.

Database:

```text
data/database/scn1a.db
```

The database contains four main tables:

```text
genes
   │
   ▼
variants
   ├── variant_representations
   └── clinvar_records
```

### Database summary

```text
Genes:                  1
Unique variants:        5,381
Genomic representations: 10,696
Raw ClinVar records:    10,696
```

### Main schema

#### `genes`

Stores gene identity.

#### `variants`

Stores the biological variant identity:

```text
variant_id
gene_id
allele_id
variation_id
hgvs_name
```

#### `clinvar_records`

Stores ClinVar evidence:

```text
variant_id
clinical_significance
review_status
number_submitters
phenotypes
```

#### `variant_representations`

Stores genomic representations:

```text
variant_id
assembly
chromosome
start
stop
reference_allele
alternate_allele
dbsnp_id
```

This separation is important because one biological variant can have multiple genomic representations.

---

# 6. Phase 3 — Exploratory Data Analysis

The project then performed six analysis stages:

```text
08_variant_statistics.py
09_clinical_significance.py
10_evidence_quality.py
11_phenotype_analysis.py
12_variant_representation.py
13_prioritization_signals.py
```

These analyses established the structure and limitations of the dataset before feature engineering.

---

# 7. Current ClinVar Landscape

After evidence deduplication, the 5,381 unique variants have the following primary classifications:

| Clinical significance | Variants |
|---|---:|
| Uncertain significance | 1,710 |
| Pathogenic | 1,380 |
| Likely pathogenic | 637 |
| Conflicting classifications | 386 |
| Pathogenic/Likely pathogenic | 213 |
| Likely benign | 788 |
| Benign | 123 |
| Benign/Likely benign | 62 |
| Not provided | 59 |
| `-` | 22 |
| Drug response | 1 |

These classifications are descriptive representations of the ClinVar data; they are not independent DeepGene predictions.

---

# 8. Evidence Reliability

DeepGene separately represents evidence reliability/context.

Current dataset:

```text
Multiple submitters:       1,592
Expert panel reviewed:        22
Conflicting classifications: 386
```

Important:

> A larger submitter count is not interpreted as a numerical measure of pathogenicity.

Instead, submitter count provides evidence context.

Relevant features include:

```text
submitter_count
multiple_submitters
expert_panel_review
conflict_flag
review_status
```

---

# 9. Phenotype Information

Phenotype information is preserved rather than discarded.

Current dataset:

```text
Variants with phenotype information: 5,381 / 5,381
```

Examples of phenotype descriptions encountered include SCN1A-related epilepsy and developmental/epileptic encephalopathy terminology.

For V1, phenotype information is retained as evidence/context.

The planned future direction is:

```text
Phenotype description
        ↓
HPO representation
        ↓
Phenotype similarity
        ↓
DeepGene feature
```

The presence of a phenotype description is **not itself treated as proof of pathogenicity**.

---

# 10. Variant Identity and QC

The project also preserves identity and representation features.

Current coverage:

```text
GRCh37 present:                  5,369 / 5,381
GRCh38 present:                  5,317 / 5,381
dbSNP present:                   4,400 / 5,381
Genomic coordinates present:     5,381 / 5,381
HGVS present:                    5,381 / 5,381
```

These are primarily used for:

- Variant identity
- Normalization
- Database integration
- Quality control
- Future joining with external datasets

They are **not automatically treated as pathogenicity predictors**.

---

# 11. Feature Engineering

The first engineered evidence dataset is:

```text
data/processed/deepgene_evidence_v1.csv
```

Dimensions:

```text
5,381 variants × 26 columns
```

The dataset contains four broad groups.

### Clinical evidence

```text
clinical_significance
pathogenic_count
pathogenic_likely_pathogenic_count
likely_pathogenic_count
benign_count
benign_likely_benign_count
likely_benign_count
vus_count
conflicting_count
```

### Reliability/context

```text
submitter_count
multiple_submitters
expert_panel_review
conflict_flag
review_status
```

### Phenotype

```text
phenotypes
phenotype_available
```

### Identity/QC

```text
variant_id
gene_id
allele_id
variation_id
hgvs_name
grch37_present
grch38_present
dbsnp_present
genomic_coordinates_present
hgvs_present
```

---

# 12. Important Deduplication Decision

The original database contains:

```text
10,696 genomic/evidence rows
```

but multiple rows can represent the same underlying evidence because a variant may have multiple genomic representations.

Therefore, DeepGene does **not** simply treat every genomic representation as an independent ClinVar assertion.

Evidence is deduplicated using:

```text
variant_id
clinical_significance
review_status
number_submitters
phenotypes
```

After deduplication:

```text
5,381 unique variants
5,381 unique evidence rows in the current database representation
```

This means:

> `total_clinvar_records` became constant at 1 in the current engineered representation and was removed from the downstream feature dataset.

The raw database remains available, so this information is not destroyed.

---

# 13. Transparent Baseline

The first DeepGene baseline is intentionally **rule-based and interpretable**.

Script:

```text
data/data_analysis/16_transparent_baseline.py
```

Output:

```text
data/processed/deepgene_baseline_v1.csv
```

The baseline does not assign arbitrary values such as:

```text
Pathogenic = 10
Likely pathogenic = 8
VUS = 3
```

Instead, it preserves the evidence interpretation directly.

Example categories:

```text
Pathogenic evidence
Likely pathogenic evidence
Pathogenic/Likely pathogenic evidence
Uncertain evidence
Conflicting evidence
Benign evidence
Likely benign evidence
Benign/Likely benign evidence
Other / insufficiently classified
```

It also reports:

```text
reliability_context
requires_evidence_review
review_group
evidence_explanation
```

This makes the baseline explainable.

---

# 14. Leakage Problem

One of the most important methodological issues in DeepGene is **target leakage**.

For example, suppose we use:

```text
clinical_significance
pathogenic_count
likely_pathogenic_count
```

as model inputs while asking the model to predict:

```text
Pathogenic vs Benign
```

The model could simply learn the ClinVar classification that it was already given.

That would not demonstrate independent predictive ability.

Therefore, DeepGene separates:

```text
Evidence used for prioritization
              ≠
Labels used for evaluation/training
```

This separation is fundamental to the research design.

---

# 15. Three-Dataset Architecture

DeepGene is designed around three distinct datasets.

## Dataset A — Evidence Dataset

```text
deepgene_evidence_v1.csv
```

Contains:

```text
all 5,381 variants
+
ClinVar evidence
+
reliability/context
+
phenotype information
+
identity/QC information
```

---

## Dataset B — Benchmark / Label Dataset

Planned:

```text
deepgene_labels_v1.csv
```

Conceptually:

```text
variant_id
label
label_confidence
label_source
```

However, the current V1 ClinVar-only dataset does **not** contain an independent external label source.

Therefore, no artificial independent labels are being created.

Current status:

```text
Independent benchmark labels: NOT AVAILABLE IN V1
```

---

## Dataset C — Final Model Matrix

Planned:

```text
deepgene_model_v1.csv
```

This will be constructed only after the evidence/label relationship and leakage rules are clearly defined.

Architecture:

```text
Evidence Dataset A
        +
Benchmark Dataset B
        ↓
Leakage-controlled feature selection
        ↓
Model Matrix C
```

---

# 16. Current Project Status

```text
Data acquisition                         ✅
SQLite database                          ✅
Exploratory analysis                     ✅
Clinical evidence analysis               ✅
Evidence reliability analysis            ✅
Phenotype analysis                       ✅
Variant representation analysis          ✅
Evidence feature engineering              ✅
Feature distribution validation           ✅
Transparent baseline                      ✅
Benchmark feasibility check              ✅
Independent benchmark labels             ⏳
Final model matrix                       ⏳
ML baseline                              ⏳
Phenotype similarity                     ⏳
GBDT                                     ⏳
MLP                                      ⏳
GNN                                      ⏳
External evidence integration            ⏳
Advanced validation                      ⏳
```

---

# 17. Complete DeepGene Roadmap

## Phase A — Foundation

### A1. Data acquisition

```text
ClinVar
   ↓
SCN1A extraction
   ↓
Raw TSV
```

### A2. Database construction

```text
Raw ClinVar
   ↓
SQLite
   ↓
genes
variants
clinvar_records
variant_representations
```

### A3. Data validation

Check:

- Number of variants
- Duplicate representations
- Clinical classifications
- Review status
- Submitter information
- Phenotypes
- Genomic representations

---

# Phase B — Evidence Engineering

### B1. Exploratory analysis

Scripts:

```text
08
09
10
11
12
13
```

### B2. Feature engineering

```text
14_feature_engineering.py
```

Output:

```text
deepgene_evidence_v1.csv
```

### B3. Feature validation

```text
15_transparent_prioritization.py
```

Current role:

```text
Feature distribution / consistency validation
```

### B4. Transparent baseline

```text
16_transparent_baseline.py
```

Output:

```text
deepgene_baseline_v1.csv
```

---

# Phase C — Benchmark Design

### C1. Label feasibility

```text
17_label_benchmark_check.py
```

Current result:

```text
Independent labels unavailable in ClinVar-only V1
```

### C2. Obtain an appropriate independent benchmark

Future work may add an independently defined benchmark source.

The source, population, time period, inclusion criteria, and confidence rules must be documented.

### C3. Create

```text
deepgene_labels_v1.csv
```

---

# Phase D — Leakage-Controlled Model Matrix

Create:

```text
deepgene_model_v1.csv
```

The model matrix should be generated from the evidence dataset and benchmark labels according to explicitly documented leakage rules.

Potential model inputs may include:

```text
Reliability features
Phenotype-derived features
Variant representation/QC features
Other future evidence sources
```

Features that directly encode the target must be handled carefully or excluded from the corresponding benchmark.

---

# Phase E — Phenotype Similarity

A future phenotype module will transform textual/ontology phenotype information into computational representations.

Planned architecture:

```text
Variant phenotype information
          ↓
Phenotype normalization
          ↓
HPO representation
          ↓
Disease/phenotype representation
          ↓
Similarity calculation
          ↓
Phenotype similarity feature
```

This should be evaluated separately before being incorporated into a larger model.

---

# Phase F — Machine Learning

The planned progression is deliberately incremental.

## F1. Simple baseline

Start with transparent models appropriate to the final label definition.

## F2. Gradient Boosted Decision Trees

```text
Evidence/features
      ↓
GBDT
      ↓
Ranking / classification
```

Evaluate using appropriate held-out data.

## F3. Multilayer Perceptron

```text
Feature vector
      ↓
MLP
      ↓
Prediction
```

## F4. Graph Neural Network

Potential graph structure:

```text
Variants
   │
   ├── phenotype relationships
   ├── genomic relationships
   ├── gene relationships
   └── future biological relationships
```

Then:

```text
Graph
  ↓
GNN
  ↓
Variant representation
  ↓
Prioritization
```

The GNN should only be introduced once the graph definition and benchmark are scientifically justified.

---

# Phase G — Model Evaluation

Models should be compared using the same properly defined benchmark.

Potential evaluation metrics include:

```text
ROC-AUC
PR-AUC
Precision
Recall
F1
Calibration
Ranking metrics
```

The appropriate metrics depend on the final task and class distribution.

The project should also inspect:

```text
False positives
False negatives
Conflicting variants
Uncertain variants
Review-status subsets
```

rather than relying on one aggregate metric.

---

# Phase H — Interpretability

DeepGene should explain why a variant received its prioritization.

Possible future tools:

```text
Feature importance
Permutation importance
SHAP
Error analysis
Phenotype contribution
Evidence contribution
```

The final system should make it possible to trace:

```text
Variant
  ↓
Evidence
  ↓
Features
  ↓
Model
  ↓
Prioritization
  ↓
Explanation
```

---

# Phase I — External Evidence Expansion

Only after the V1 pipeline is stable should additional data sources be introduced.

Possible future evidence layers:

```text
Population evidence
Functional assays
Conservation
Splicing predictions
Protein impact
Protein structure
Literature
Phenotype ontologies
Patient-level data
```

Each source should be added as a separate, documented evidence layer.

---

# Phase J — Research Validation

Future validation should include:

- Independent benchmark evaluation
- Leakage testing
- Ablation studies
- Feature importance analysis
- Error analysis
- Robustness checks
- Subgroup analysis
- Reproducibility checks
- Versioned datasets
- Reproducible preprocessing

The objective is not merely to obtain a high score, but to determine **what information actually contributes to variant prioritization**.

---

# 18. Planned Final Architecture

The long-term DeepGene architecture is:

```text
                         ┌─────────────────────┐
                         │   SCN1A Variants    │
                         └──────────┬──────────┘
                                    │
                 ┌──────────────────┼──────────────────┐
                 │                  │                  │
                 ▼                  ▼                  ▼
          ClinVar Evidence     Phenotypes       Variant Identity
                 │                  │                  │
                 ▼                  ▼                  ▼
           Evidence Layer     HPO/Semantics       QC Layer
                 │                  │                  │
                 └──────────────────┼──────────────────┘
                                    │
                                    ▼
                          Feature Representation
                                    │
                       ┌────────────┴────────────┐
                       │                         │
                       ▼                         ▼
                Transparent Baseline       ML Models
                       │               ┌────────┼────────┐
                       │               ▼        ▼        ▼
                       │              GBDT      MLP      GNN
                       │
                       └────────────┬────────────┘
                                    ▼
                            Variant Prioritization
                                    │
                                    ▼
                              Explanation
```

---

# 19. Planned Repository Structure

Current/recommended structure:

```text
Project_DeepGene/
│
├── LICENSE
├── README.md
├── .gitignore
│
└── data/
    │
    ├── raw/
    │   └── SCN1A_clinvar.tsv
    │
    ├── database/
    │   ├── 05_create_database.py
    │   ├── 06_load_clinvar.py
    │   ├── 07_query_database.py
    │   └── scn1a.db
    │
    ├── data_analysis/
    │   ├── 08_variant_statistics.py
    │   ├── 09_clinical_significance.py
    │   ├── 10_evidence_quality.py
    │   ├── 11_phenotype_analysis.py
    │   ├── 12_variant_representation.py
    │   ├── 13_prioritization_signals.py
    │   ├── 14_feature_engineering.py
    │   ├── 15_transparent_prioritization.py
    │   ├── 16_transparent_baseline.py
    │   └── 17_label_benchmark_check.py
    │
    ├── analysis_results/
    │   ├── 08_variant_statistics.txt
    │   ├── 09_clinical_significance.txt
    │   ├── 10_evidence_quality.txt
    │   ├── 11_phenotype_analysis.txt
    │   ├── 12_variant_representation.txt
    │   ├── 13_prioritization_signals.txt
    │   └── 16_transparent_baseline.txt
    │
    └── processed/
        ├── SCN1A_clinvar_clean.csv
        ├── deepgene_evidence_v1.csv
        └── deepgene_baseline_v1.csv
```

Future additions may include:

```text
models/
notebooks/
src/
tests/
docs/
```

as the project grows.

---

# 20. Scientific Design Principles

DeepGene follows these principles:

### 1. Transparency first

Every major transformation should be understandable.

### 2. No invented evidence

Unavailable evidence is not silently converted into zero.

### 3. Preserve raw information

The raw database remains available even when downstream representations aggregate information.

### 4. Separate evidence from labels

This is essential for preventing leakage.

### 5. No arbitrary scores without justification

A numerical score should have a documented statistical or scientific basis.

### 6. Complexity should be earned

A GNN is not automatically better than a transparent baseline. More complex models should only be introduced when justified by the data and benchmark.

### 7. Reproducibility

Every dataset and transformation should be traceable to its source and processing step.

### 8. V1 is not a clinical system

DeepGene is a research and computational prioritization project.

---

# 21. Current Key Outputs

The most important current artifacts are:

```text
data/database/scn1a.db
data/processed/deepgene_evidence_v1.csv
data/processed/deepgene_baseline_v1.csv
```

Together they represent the current DeepGene V1 evidence pipeline.

---

# 22. Current Milestone

### DeepGene V1 Evidence Layer — Completed

```text
ClinVar extraction                    ✅
SQLite database                       ✅
SCN1A variant representation          ✅
Exploratory analysis                  ✅
Evidence feature engineering          ✅
Feature validation                    ✅
Transparent baseline                  ✅
Benchmark feasibility assessment      ✅
```

### Next milestone

```text
Independent benchmark design
        ↓
Leakage-controlled model matrix
        ↓
Phenotype similarity
        ↓
ML baseline
        ↓
GBDT / MLP
        ↓
GNN
        ↓
Independent validation
```

---

## DeepGene in one sentence

> **DeepGene is an interpretable, reproducible research framework that starts from SCN1A ClinVar evidence and progressively develops toward phenotype-aware and machine-learning-based genetic variant prioritization while explicitly controlling for missing evidence, data quality, and target leakage.**

---

# 23. Generation-2 Functional Model Status

The Generation-2 functional branch now uses the published SCION cross-channel
functional dataset covering 375 variants across nine sodium-channel genes.
After removing 40 exact protein-level overlaps with the DeepGene V1 functional
dataset, 335 variants remain in the current cross-channel training matrix.

Implemented artifacts:

```text
data/raw/scion_clean_tbl.csv
data/raw/generation2_source_manifest.json
data/processed/deepgene_gen2_scion_harmonized_v1.csv
data/processed/deepgene_gen2_ml_ready_v1.csv
data/model/deepgene_generation2_gbm.pkl
```

The current Generation-2 candidate predicts functional `LOF` versus `GOF`.
Its primary validation is leave-one-gene-out validation. The current Gradient
Boosting candidate achieved balanced accuracy of approximately `0.6368 +/-
0.1465`; this is a research benchmark, not clinical validation or an
independent external test.

Run a Generation-2 functional prediction with:

```powershell
python data/ml/39_predict_generation2_effect.py --gene SCN1A --variant p.Arg1234Gly
```

# 24. Dravet-Syndrome Model Status

The Dravet model is a separate, not-yet-trained clinical-outcome branch.
Functional `LOF/GOF` labels cannot be used as a Dravet target because SCN1A
loss-of-function is associated with a spectrum of phenotypes, not Dravet
syndrome alone. A Dravet model requires clinically curated variant-level labels
and phenotype information.

The safe workflow is:

```text
Reviewed Dravet/non-Dravet labels
              ↓
Leakage-controlled feature matrix
              ↓
Patient/cohort-aware validation
              ↓
Independent clinical test set
              ↓
Research-only prediction model
```

Create the required label template:

```powershell
python data/ml/45_create_dravet_label_template.py
```

Then populate:

```text
data/raw/dravet_clinical_labels.csv
```

Required fields are `variant_key`, `dravet_label`, `label_source`,
`source_accession`, `cohort_id`, and `reviewer_status`. Unknown or unresolved
cases must not be assigned a negative label. The preparation and training
commands are:

```powershell
python data/ml/46_prepare_dravet_dataset.py
python data/ml/47_train_dravet_model.py
```

Until reviewed labels are supplied, the Dravet branch intentionally stops and
does not train. No current DeepGene output should be interpreted as a Dravet
diagnosis, pathogenicity classification, disease-risk probability, or treatment
recommendation.

If a variant has no reviewed Dravet label, the non-crashing status is expected:

```powershell
python data/ml/47_train_dravet_model.py
```

For documented phenotype mentions in the local ClinVar-derived dataset, use:

```powershell
python data/ml/48_scn1a_syndrome_evidence.py --variant p.Arg1648His
```

This reports observed syndrome terms and whether Dravet-related wording is
present. It intentionally reports no percentage when a validated penetrance or
clinical-outcome dataset is unavailable.

---

# 25. Local Web Dashboard

DeepGene now includes a local two-page HTML/CSS/JavaScript dashboard.

Start it from the project root:

```powershell
python app/server.py
```

Open:

```text
http://127.0.0.1:8765
```

Pages:

```text
/             AI analysis and local project statistics
/upload.html  dataset validation, SQLite import, and functional retraining
```

The dashboard works without an external AI service. It answers basic project
questions locally from SQLite and current CSV artifacts. Optional natural-
language answers can be enabled by setting an `OPENROUTER_API_KEY`; the key is
read only by the local server and is never sent to the browser. The server uses
OpenRouter's `openrouter/free` router by default, subject to provider
availability and rate limits.

Supported upload workflows:

```text
ClinVar CSV/TSV
    -> column validation
    -> transactional SQLite update
    -> upload audit log

Functional CSV/TSV
    -> LOF/GOF validation
    -> cleaned feature construction
    -> research-model retraining
    -> upload audit log
```

The dashboard does not automatically create Dravet labels. Dravet model
training remains blocked until clinically reviewed labels are supplied through
the dedicated label template.
