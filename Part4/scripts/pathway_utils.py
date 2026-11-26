"""Shared inputs for Part 4. Run the analysis from the Part4 directory."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

LIBRARIES = {
    "MSigDB_Hallmark_2020": (50, "4275592957a1587652092bb398cf7721"),
    "Reactome_Pathways_2024": (2104, "1a7c43cb503d9db23494605169614a7b"),
    "GO_Biological_Process_2025": (5342, "0fc02458765cfce6c9348dce9f7a9397"),
}
URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"

def rank_genes(de):
    d = de.dropna(subset=["symbol", "stat"]).copy()
    d = d[np.isfinite(d.stat) & d.symbol.astype(str).str.strip().ne("")]
    d["abs"] = d.stat.abs()
    d = d.sort_values("abs", ascending=False).drop_duplicates("symbol")
    return d[["symbol", "stat"]].sort_values("stat", ascending=False).reset_index(drop=True)

def get_gmt(name):
    """Check the exact archived library by its documented SHA-256 prefix."""
    path = Path("genesets") / (name + ".gmt")
    expected_count, prefix = LIBRARIES[name]
    def validate(text):
        rows = [line.split("\t") for line in text.splitlines() if line.strip()]
        digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        if len(rows) != expected_count or any(len(row) < 3 for row in rows) or not digest.startswith(prefix):
            raise ValueError(f"{name}: library differs from the archived version (sets={len(rows)}, SHA-256={digest}). Use the original GMT; do not silently substitute another version.")
    if path.exists():
        text = path.read_text(encoding="utf-8")
        validate(text)
    else:
        import requests
        response = requests.get(URL.format(name), timeout=180)
        response.raise_for_status()
        text = response.content.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
        validate(text)
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(text.encode("utf-8"))
    return str(path)

def new_dir(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"{path} is not empty; choose a fresh output directory.")
    path.mkdir(parents=True, exist_ok=True)
    return path

def lead_genes(value):
    return [] if pd.isna(value) else [g for g in str(value).split(";") if g]
