# DeepGene local dashboard

From the project root:

```powershell
python app/server.py
```

Open <http://127.0.0.1:8765>.

The app has two pages:

- `/` — local statistics and optional AI question answering.
- `/upload.html` — CSV/TSV validation, SQLite ClinVar import, and functional-model retraining.

Optional language answers use OpenRouter only when `OPENROUTER_API_KEY` is set.
Without a key, the local analyzer remains available.

## Conversation memory and learning

Each chat turn is saved locally in `data/memory/conversations.jsonl` and the
latest remembered turns are supplied as context to optional external AI calls.
This is retrieval-style conversation memory; it does not update the weights of
the OpenRouter model. The memory file is ignored by Git because it may contain
research notes or sensitive text.

Functional model retraining is separate and requires a labelled CSV/TSV with
`gene, aa1, aa2, pos, y`, where `y` is `LOF` or `GOF`. Uploaded rows are
validated and combined with the existing training data to create a research
candidate model. Ordinary chat messages are never used as training labels.

## SCN1A benchmark refresh

The project includes a public supplementary SCN1A functional dataset from
Brunklaus et al. (2020), downloaded from UCL Discovery. The reproducible
preparation and training commands are:

```powershell
python data/ml/49_prepare_scn1a_supplement.py
python data/ml/49_train_scn1a_supplement.py
```

The resulting model is
`data/model/deepgene_scn1a_supplement_candidate.pkl`. It is a research-only
LOF/GOF candidate and excludes mixed or uncertain functional labels.
