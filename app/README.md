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
