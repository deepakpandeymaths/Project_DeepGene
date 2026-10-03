"""DeepGene local dashboard server.

Run from the project root:
    python app/server.py

The server keeps all external AI access optional. Statistics are generated
locally from SQLite/CSV; OpenRouter is used only when OPENROUTER_API_KEY is
present.
"""

from __future__ import annotations
import csv
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
MODEL = ROOT / "data" / "model" / "deepgene_generation2_gbm_app.pkl"
UPLOADS.mkdir(parents=True, exist_ok=True)

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
    for path in [ROOT / "data/processed/deepgene_v1_final.csv", ROOT / "data/processed/deepgene_gen2_ml_ready_v1.csv"]:
        if path.exists():
            result["csv_files"][path.name] = int(sum(1 for _ in path.open(encoding="utf-8")) - 1)
    for path in [ROOT / "data/model/deepgene_generation2_gbm.pkl", ROOT / "data/model/deepgene_dravet_research_candidate.pkl"]:
        result["models"][path.name] = path.exists()
    return result


def local_answer(question: str):
    summary = db_summary()
    q = question.lower()
    if "how many" in q and ("variant" in q or "record" in q):
        return f"The SQLite database currently contains {summary['tables'].get('variants')} unique variants and {summary['tables'].get('clinvar_records')} ClinVar records."
    if "gene" in q and "how many" in q:
        return f"The database currently contains {summary['tables'].get('genes')} gene entries."
    if "model" in q or "trained" in q:
        return f"Available model artifacts: {', '.join(k for k,v in summary['models'].items() if v) or 'none'}. The Dravet model is not trained without curated labels."
    return "Local analysis is available for database counts, uploaded-file status, and model artifacts. Ask a specific question such as: How many variants are in the database?"


def ai_answer(question: str):
    local = local_answer(question)
    key = os.getenv("OPENROUTER_API_KEY")
    if not key:
        return {"answer": local, "provider": "local", "note": "Set OPENROUTER_API_KEY to enable optional free-router language answers."}
    context = json.dumps(db_summary(), indent=2)
    payload = json.dumps({"model": os.getenv("OPENROUTER_MODEL", "openrouter/free"), "messages": [{"role": "system", "content": "You are a cautious research assistant for DeepGene. Use only the supplied project context. Do not diagnose, estimate disease risk, or invent percentages. Distinguish functional-effect prediction from clinical interpretation."}, {"role": "user", "content": f"Project context:\n{context}\n\nQuestion: {question}"}], "temperature": 0.1}).encode()
    request = urllib.request.Request("https://openrouter.ai/api/v1/chat/completions", data=payload, headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json", "HTTP-Referer": "http://localhost:8765", "X-Title": "DeepGene local dashboard"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            data = json.loads(response.read().decode())
        return {"answer": data["choices"][0]["message"]["content"], "provider": "openrouter"}
    except Exception as exc:
        return {"answer": local, "provider": "local", "note": f"Optional AI provider unavailable: {type(exc).__name__}."}


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
    suffix = Path(filename).suffix.lower()
    if suffix not in {".csv", ".tsv", ".txt"}:
        raise ValueError("Only CSV or TSV files are accepted.")
    text = payload.decode("utf-8-sig", errors="strict")
    sep = "\t" if suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(io.StringIO(text), sep=sep, low_memory=False)


def detect_kind(df, requested):
    if requested in {"clinvar", "functional"}:
        return requested
    if {"#AlleleID", "VariationID", "Name", "ClinicalSignificance"}.issubset(df.columns): return "clinvar"
    if {"gene", "aa1", "aa2", "pos"}.issubset(df.columns) and ("y" in df.columns or "functional_label" in df.columns): return "functional"
    return "unknown"


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
                length = int(self.headers.get("Content-Length", 0)); body = json.loads(self.rfile.read(length)); return self.send_json(ai_answer(str(body.get("question", ""))))
            if path == "/api/upload":
                filename, requested, payload = read_upload(self); df = parse_table(filename, payload); kind = detect_kind(df, requested)
                if kind == "unknown": raise ValueError("Could not recognize the file. Choose ClinVar or functional data and use the documented columns.")
                destination = UPLOADS / filename; destination.write_bytes(payload)
                result = ingest_clinvar(df) if kind == "clinvar" else retrain_functional(df)
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
