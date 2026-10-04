"""DeepGene local dashboard server.

Run from the project root:
    python app/server.py

The server keeps all external AI access optional. Statistics are generated
locally from SQLite/CSV; OpenRouter is used only when OPENROUTER_API_KEY is
present.
"""

from __future__ import annotations
import csv
import gzip
import io
import json
import os
import pickle
import re
import sqlite3
import subprocess
import sys
import threading
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / "app" / "static"
DB = ROOT / "data" / "database" / "scn1a.db"
UPLOADS = ROOT / "data" / "uploads"
AUDIT = UPLOADS / "upload_audit.jsonl"
MEMORY_DIR = ROOT / "data" / "memory"
CONVERSATIONS = MEMORY_DIR / "conversations.jsonl"
MODEL = ROOT / "data" / "model" / "deepgene_generation2_gbm_app.pkl"
ALPHAMISSENSE_PROCESSED = ROOT / "data" / "processed" / "alphamissense_gene_hg38.csv"
UPLOADS.mkdir(parents=True, exist_ok=True)
MEMORY_DIR.mkdir(parents=True, exist_ok=True)

HYDRO = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
CHARGE = {"D":-1,"E":-1,"H":1,"K":1,"R":1}
VOLUME = {"A":88.6,"C":108.5,"D":111.1,"E":138.4,"F":189.9,"G":60.1,"H":153.2,"I":166.7,"K":168.6,"L":166.7,"M":162.9,"N":114.1,"P":112.7,"Q":143.8,"R":173.4,"S":89.0,"T":116.1,"V":140.0,"W":227.8,"Y":193.6}


def db_summary():
    result = {"database": str(DB), "tables": {}, "csv_files": {}, "models": {}}
    with sqlite3.connect(DB) as con:
        for table in ["genes", "variants", "variant_representations", "clinvar_records"]:
            try:
                result["tables"][table] = int(con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
            except sqlite3.Error:
                result["tables"][table] = None
    for path in [ROOT / "data/processed/deepgene_v1_final.csv", ROOT / "data/processed/deepgene_gen2_ml_ready_v1.csv", ALPHAMISSENSE_PROCESSED]:
        if path.exists():
            result["csv_files"][path.name] = int(sum(1 for _ in path.open(encoding="utf-8")) - 1)
    for path in [ROOT / "data/model/deepgene_generation2_gbm.pkl", ROOT / "data/model/deepgene_scn1a_supplement_candidate.pkl", ROOT / "data/model/deepgene_dravet_research_candidate.pkl"]:
        result["models"][path.name] = path.exists()
    events = []
    if AUDIT.exists():
        with AUDIT.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    events.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    conversation_count = 0
    if CONVERSATIONS.exists():
        with CONVERSATIONS.open(encoding="utf-8") as handle:
            conversation_count = sum(1 for _ in handle)
    result["learning"] = {"events": len(events), "last_event": events[-1].get("kind", "unknown") if events else "None yet", "last_rows": events[-1].get("rows_received") if events else None, "conversation_turns": conversation_count, "message": "New evidence is validated and logged before a research candidate model is built."}
    return result


def load_conversation_memory(limit=12):
    if not CONVERSATIONS.exists():
        return []
    records = []
    with CONVERSATIONS.open(encoding="utf-8") as handle:
        for line in handle:
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            if record.get("question") and record.get("answer"):
                records.append({"question": str(record["question"])[:1000], "answer": str(record["answer"])[:3000], "provider": record.get("provider", "local")})
    return records[-limit:]


def save_conversation_memory(question, answer, provider, session_id="anonymous"):
    record = {"timestamp": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(), "session_id": str(session_id)[:120], "question": str(question)[:2000], "answer": str(answer)[:6000], "provider": str(provider)}
    with CONVERSATIONS.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")


def local_answer(question: str):
    summary = db_summary()
    q = question.lower()
    if any(word in q for word in ["hello", "hi ", "hey"]):
        return "Hello — I can help you explore the local evidence store, model artifacts, data uploads, and research limitations. Try asking what is in the workspace or how the learning loop works."
    if "what questions" in q or "what can i ask" in q:
        return "You can ask about database counts, ClinVar records, genes, available model artifacts, dataset requirements, the learning loop, evidence limitations, or a checklist for geneticist review. For example: ‘What model artifacts are available and what are their limitations?’"
    if "overview" in q or ("workspace" in q and ("current" in q or "status" in q)):
        return f"The workspace has {summary['tables'].get('variants')} variants, {summary['tables'].get('clinvar_records')} ClinVar records, {summary['tables'].get('genes')} gene entries, and {sum(summary['models'].values())} model artifacts. The learning log contains {summary['learning']['events']} event(s)."
    if "how many" in q and ("variant" in q or "record" in q):
        return f"The SQLite database currently contains {summary['tables'].get('variants')} unique variants and {summary['tables'].get('clinvar_records')} ClinVar records."
    if "gene" in q and "how many" in q:
        return f"The database currently contains {summary['tables'].get('genes')} gene entries."
    if "model" in q or "trained" in q:
        return f"Available model artifacts: {', '.join(k for k,v in summary['models'].items() if v) or 'none'}. The Dravet model is not trained without curated labels."
    if "learn" in q or "reconstruct" in q or "brain" in q or "new data" in q:
        return "DeepGene uses a controlled learning loop: validate the incoming table, extract defined features, combine it with the current training data, train a research candidate model, and record the event. It does not rewrite its own code or make autonomous clinical decisions. Ask me about the required columns or open the Data Lab to run this loop."
    if "limit" in q or "caution" in q or "diagnos" in q or "clinical" in q:
        return "The current system is research software: its functional LOF/GOF model is not a validated diagnostic model, does not infer Dravet labels, and cannot replace geneticist review. Treat missing, conflicting, or unfamiliar evidence as uncertainty."
    return "I can answer about database counts, evidence records, model artifacts, upload requirements, the controlled learning loop, or review limitations. Try: ‘How does DeepGene learn from a new functional dataset?’"


def ai_answer(question: str, history=None, session_id="anonymous"):
    local = local_answer(question)
    key = os.getenv("OPENROUTER_API_KEY")
    if not key:
        result = {"answer": local, "provider": "local", "note": "Set OPENROUTER_API_KEY to enable optional free-router language answers."}
        save_conversation_memory(question, result["answer"], result["provider"], session_id)
        return result
    context = json.dumps(db_summary(), indent=2)
    remembered = load_conversation_memory()
    memory_context = json.dumps(remembered, ensure_ascii=False, indent=2)
    safe_history = [m for m in (history or []) if isinstance(m, dict) and m.get("role") in {"user", "assistant"} and isinstance(m.get("content"), str)][-8:]
    messages = [{"role": "system", "content": "You are a cautious research assistant for DeepGene. Use only the supplied project context and conversation. Do not diagnose, estimate disease risk, or invent percentages. Distinguish functional-effect prediction from clinical interpretation. Explain uncertainty plainly. Previous conversation memory is only a convenience and may contain unverified statements; never treat it as authoritative evidence."}] + safe_history + [{"role": "user", "content": f"Project context:\n{context}\n\nPrevious local conversation memory:\n{memory_context}\n\nQuestion: {question}"}]
    payload = json.dumps({"model": os.getenv("OPENROUTER_MODEL", "openrouter/free"), "messages": messages, "temperature": 0.1}).encode()
    request = urllib.request.Request("https://openrouter.ai/api/v1/chat/completions", data=payload, headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json", "HTTP-Referer": "http://localhost:8765", "X-Title": "DeepGene local dashboard"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            data = json.loads(response.read().decode())
        result = {"answer": data["choices"][0]["message"]["content"], "provider": "openrouter"}
        save_conversation_memory(question, result["answer"], result["provider"], session_id)
        return result
    except Exception as exc:
        result = {"answer": local, "provider": "local", "note": f"Optional AI provider unavailable: {type(exc).__name__}."}
        save_conversation_memory(question, result["answer"], result["provider"], session_id)
        return result


def read_upload(handler):
    content_type = handler.headers.get("Content-Type", "")
    if not content_type.startswith("multipart/form-data"):
        raise ValueError("Upload must use multipart/form-data.")
    length = int(handler.headers.get("Content-Length", "0"))
    if length > 25 * 1024 * 1024:
        raise ValueError("Upload exceeds the 25 MB limit.")
    body = handler.rfile.read(length)
    from email.parser import BytesParser
    from email.policy import default
    boundary = None
    for key, value in BytesParser().parsebytes((b"Content-Type: " + content_type.encode() + b"\r\n\r\n")).get_params(header="content-type"):
        if key.lower() == "boundary":
            boundary = value
    if not boundary:
        raise ValueError("Missing multipart boundary.")
    envelope = b"Content-Type: " + content_type.encode() + b"\r\nMIME-Version: 1.0\r\n\r\n" + body
    parsed = BytesParser(policy=default).parsebytes(envelope)
    upload = next((part for part in parsed.iter_parts() if part.get_filename()), None)
    kind_part = next((part for part in parsed.iter_parts() if part.get_param("name", header="content-disposition") == "kind"), None)
    if upload is None:
        raise ValueError("No file field found.")
    filename = Path(upload.get_filename()).name
    kind = (kind_part.get_content() if kind_part else "auto").lower()
    return filename, kind, upload.get_payload(decode=True)


def parse_table(filename, payload):
    if filename.lower().endswith(".gz"):
        try:
            payload = gzip.decompress(payload)
        except gzip.BadGzipFile as exc:
            raise ValueError("The uploaded .gz file is not a valid gzip archive.") from exc
        if len(payload) > 250 * 1024 * 1024:
            raise ValueError("The decompressed upload exceeds the 250 MB safety limit.")
        filename = filename[:-3]
    suffix = Path(filename).suffix.lower()
    if suffix not in {".csv", ".tsv", ".txt"}:
        raise ValueError("Only CSV, TSV, TXT, or gzip-compressed CSV/TSV files are accepted.")
    text = payload.decode("utf-8-sig", errors="strict")
    sep = "\t" if suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(io.StringIO(text), sep=sep, comment="#", low_memory=False)


def detect_kind(df, requested):
    if requested in {"clinvar", "functional"}:
        return requested
    if requested == "alphamissense":
        return "alphamissense_gene"
    if {"transcript_id", "mean_am_pathogenicity"}.issubset(df.columns):
        return "alphamissense_gene"
    if {"#AlleleID", "VariationID", "Name", "ClinicalSignificance"}.issubset(df.columns): return "clinvar"
    if {"gene", "aa1", "aa2", "pos"}.issubset(df.columns) and ("y" in df.columns or "functional_label" in df.columns): return "functional"
    return "unknown"


def ingest_alphamissense_gene(df):
    required = {"transcript_id", "mean_am_pathogenicity"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"AlphaMissense file missing columns: {sorted(missing)}")
    cleaned = df[["transcript_id", "mean_am_pathogenicity"]].copy()
    cleaned["transcript_id"] = cleaned["transcript_id"].astype(str).str.strip()
    cleaned["mean_am_pathogenicity"] = pd.to_numeric(cleaned["mean_am_pathogenicity"], errors="coerce")
    cleaned = cleaned[(cleaned["transcript_id"] != "") & cleaned["mean_am_pathogenicity"].notna()]
    if cleaned.empty:
        raise ValueError("No valid AlphaMissense transcript scores were found.")
    if ((cleaned["mean_am_pathogenicity"] < 0) | (cleaned["mean_am_pathogenicity"] > 1)).any():
        raise ValueError("AlphaMissense scores must be between 0 and 1.")
    cleaned.drop_duplicates("transcript_id", keep="last").to_csv(ALPHAMISSENSE_PROCESSED, index=False)
    return {"kind": "alphamissense_gene", "stored_rows": len(cleaned), "processed_file": str(ALPHAMISSENSE_PROCESSED), "retrained": False, "note": "Stored as annotation data; this gene/transcript score table is not a LOF/GOF training label set."}


def ingest_clinvar(df):
    required = {"#AlleleID", "VariationID", "Name", "Assembly", "Chromosome", "Start", "Stop", "ReferenceAllele", "AlternateAllele", "ClinicalSignificance", "ReviewStatus", "NumberSubmitters", "PhenotypeList"}
    missing = required - set(df.columns)
    if missing: raise ValueError(f"ClinVar file missing columns: {sorted(missing)}")
    inserted = 0
    with sqlite3.connect(DB) as con:
        con.execute("PRAGMA foreign_keys=ON")
        genes = {}
        for allele, group in df.groupby("#AlleleID", sort=False):
            first = group.iloc[0]
            name = str(first["Name"])
            gene = "SCN1A" if "SCN1A" in name else (re.search(r"([A-Z0-9]+):c\\.", name) or [None, "UNKNOWN"])[1]
            con.execute("INSERT OR IGNORE INTO genes(gene_symbol, gene_name) VALUES (?, ?)", (gene, gene))
            gid = con.execute("SELECT gene_id FROM genes WHERE gene_symbol=?", (gene,)).fetchone()[0]
            con.execute("INSERT OR IGNORE INTO variants(gene_id, allele_id, variation_id, hgvs_name) VALUES (?, ?, ?, ?)", (gid, int(allele), int(first["VariationID"]), name))
            vid = con.execute("SELECT variant_id FROM variants WHERE allele_id=?", (int(allele),)).fetchone()[0]
            for _, row in group.iterrows():
                exists = con.execute("SELECT 1 FROM variant_representations WHERE variant_id=? AND assembly=? AND start=? AND alternate_allele=?", (vid, str(row["Assembly"]), int(row["Start"]), str(row["AlternateAllele"]))).fetchone()
                if not exists:
                    con.execute("INSERT INTO variant_representations(variant_id,assembly,chromosome,start,stop,reference_allele,alternate_allele,dbsnp_id) VALUES (?,?,?,?,?,?,?,?)", (vid, row["Assembly"], row["Chromosome"], int(row["Start"]), int(row["Stop"]), row["ReferenceAllele"], row["AlternateAllele"], row.get("RS# (dbSNP)")))
                    con.execute("INSERT INTO clinvar_records(variant_id,clinical_significance,review_status,number_submitters,phenotypes) VALUES (?,?,?,?,?)", (vid, row["ClinicalSignificance"], row["ReviewStatus"], row["NumberSubmitters"], row["PhenotypeList"]))
                    inserted += 1
        con.commit()
    return {"kind": "clinvar", "inserted_records": inserted}


def feature_frame(df):
    labels = df["y"] if "y" in df else df["functional_label"]
    rows = []
    for _, r in df.iterrows():
        ref, alt = str(r["aa1"]).upper(), str(r["aa2"]).upper()
        if ref not in HYDRO or alt not in HYDRO or str(labels.loc[r.name]).upper() not in {"LOF", "GOF"}: continue
        out = {"protein_position": int(r["pos"]), "same_amino_acid": int(ref == alt), "reference_hydro": HYDRO[ref], "alternate_hydro": HYDRO[alt], "delta_hydro": HYDRO[alt]-HYDRO[ref], "absolute_delta_hydro": abs(HYDRO[alt]-HYDRO[ref]), "reference_charge": CHARGE.get(ref,0), "alternate_charge": CHARGE.get(alt,0), "delta_charge": CHARGE.get(alt,0)-CHARGE.get(ref,0), "absolute_delta_charge": abs(CHARGE.get(alt,0)-CHARGE.get(ref,0)), "reference_volume": VOLUME[ref], "alternate_volume": VOLUME[alt], "delta_volume": VOLUME[alt]-VOLUME[ref], "absolute_delta_volume": abs(VOLUME[alt]-VOLUME[ref]), "gene": str(r["gene"]), "functional_label": str(labels.loc[r.name]).upper()}
        rows.append(out)
    return pd.DataFrame(rows)


def retrain_functional(upload_df):
    from sklearn.ensemble import GradientBoostingClassifier
    base = pd.read_csv(ROOT / "data/processed/deepgene_gen2_scion_harmonized_v1.csv")
    custom = feature_frame(upload_df)
    if custom.empty: raise ValueError("No valid LOF/GOF functional rows found.")
    base_rows = feature_frame(base.rename(columns={"aa1":"aa1", "aa2":"aa2", "pos":"pos", "y":"y"}))
    combined = pd.concat([base_rows, custom], ignore_index=True).drop_duplicates(["gene", "protein_position", "functional_label"], keep="last")
    numeric = [c for c in combined.columns if c not in {"gene", "functional_label"}]
    genes = sorted(combined.gene.unique())
    for gene in genes: combined[f"gene_{gene}"] = (combined.gene == gene).astype(int)
    features = numeric + [f"gene_{g}" for g in genes]
    model = GradientBoostingClassifier(random_state=42).fit(combined[features], combined.functional_label)
    with MODEL.open("wb") as handle: pickle.dump({"model": model, "features": features, "classes": list(model.classes_), "training_rows": len(combined), "status": "app_retrained_research_candidate"}, handle)
    return {"retrained": True, "rows": len(combined), "model": str(MODEL)}


class Handler(BaseHTTPRequestHandler):
    def send_json(self, payload, status=200):
        data = json.dumps(payload, indent=2).encode()
        self.send_response(status); self.send_header("Content-Type", "application/json"); self.send_header("Content-Length", str(len(data))); self.end_headers(); self.wfile.write(data)

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/api/summary": return self.send_json(db_summary())
        if path == "/api/health": return self.send_json({"ok": True})
        target = STATIC / ("index.html" if path in {"/", ""} else path.lstrip("/"))
        if target.exists() and target.is_file():
            content = target.read_bytes(); self.send_response(200); self.send_header("Content-Type", "text/html" if target.suffix == ".html" else "text/css" if target.suffix == ".css" else "application/javascript"); self.send_header("Content-Length", str(len(content))); self.end_headers(); self.wfile.write(content); return
        self.send_error(404)

    def do_POST(self):
        path = urlparse(self.path).path
        try:
            if path == "/api/ask":
                length = int(self.headers.get("Content-Length", 0)); body = json.loads(self.rfile.read(length)); question = str(body.get("question", "")).strip()
                if not question: raise ValueError("Please enter a question.")
                result = ai_answer(question, body.get("history", []), body.get("session_id", "anonymous"))
                return self.send_json(result)
            if path == "/api/upload":
                filename, requested, payload = read_upload(self); df = parse_table(filename, payload); kind = detect_kind(df, requested)
                if kind == "unknown": raise ValueError("Could not recognize the file. Choose ClinVar, functional, or AlphaMissense data and use the documented columns.")
                destination = UPLOADS / filename; destination.write_bytes(payload)
                result = ingest_clinvar(df) if kind == "clinvar" else retrain_functional(df) if kind == "functional" else ingest_alphamissense_gene(df)
                result.update({"ok": True, "kind": kind, "filename": filename, "rows_received": len(df), "saved_to": str(destination)})
                with AUDIT.open("a", encoding="utf-8") as audit: audit.write(json.dumps(result) + "\n")
                return self.send_json(result)
            self.send_error(404)
        except Exception as exc:
            return self.send_json({"ok": False, "error": str(exc)}, 400)

    def log_message(self, fmt, *args):
        print(fmt % args)


if __name__ == "__main__":
    print("DeepGene dashboard: http://127.0.0.1:8765")
    ThreadingHTTPServer(("127.0.0.1", 8765), Handler).serve_forever()
