"""Pathway / gene-set libraries in GMT form.

The Hallmark library ships in `backend/data/gmt/`. The others are downloaded once
from Enrichr and cached under `$SJVC_CACHE_DIR/gmt` (default `~/.cache/sjv/gmt`).
No `gseapy` dependency — a GMT line is just  `term \t description \t gene...`.

All libraries are human gene symbols; pathway mode is disabled for mouse.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import httpx

BACKEND_DIR = Path(__file__).resolve().parents[2]
BUNDLED_DIR = BACKEND_DIR / "data" / "gmt"
CACHE_DIR = Path(
    os.environ.get("SJVC_CACHE_DIR") or os.environ.get("SJV_CACHE_DIR") or Path.home() / ".cache" / "sjv"
) / "gmt"

_ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={name}"

# label shown in the UI -> Enrichr library name
LIBRARIES: Dict[str, str] = {
    "MSigDB Hallmark": "MSigDB_Hallmark_2020",
    "KEGG (human)": "KEGG_2021_Human",
    "Reactome": "Reactome_2022",
    "WikiPathways (human)": "WikiPathway_2023_Human",
    "GO Biological Process": "GO_Biological_Process_2023",
}


class PathwayError(ValueError):
    pass


@dataclass
class Library:
    label: str
    name: str
    terms: Dict[str, List[str]]   # term -> gene symbols


_loaded: Dict[str, Library] = {}


def list_libraries() -> List[str]:
    return list(LIBRARIES)


def _parse_gmt(text: str) -> Dict[str, List[str]]:
    terms: Dict[str, List[str]] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term = parts[0].strip()
        genes = [g.split(",")[0].strip() for g in parts[2:] if g.strip()]
        if term and genes:
            terms[term] = genes
    return terms


def _load(label: str) -> Library:
    if label in _loaded:
        return _loaded[label]
    if label not in LIBRARIES:
        raise PathwayError(f"unknown pathway library {label!r}")
    name = LIBRARIES[label]

    bundled = BUNDLED_DIR / f"{name}.gmt"
    cached = CACHE_DIR / f"{name}.gmt"
    text: Optional[str] = None
    if bundled.exists():
        text = bundled.read_text()
    elif cached.exists():
        text = cached.read_text()
    else:
        try:
            r = httpx.get(_ENRICHR.format(name=name), timeout=60, follow_redirects=True)
            r.raise_for_status()
            text = r.text
        except httpx.HTTPError as e:
            raise PathwayError(f"could not download library {label!r}: {e}")
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cached.write_text(text)

    lib = Library(label=label, name=name, terms=_parse_gmt(text))
    if not lib.terms:
        raise PathwayError(f"library {label!r} parsed to zero terms")
    _loaded[label] = lib
    return lib


def search(label: str, query: str, limit: int = 25) -> List[Dict[str, object]]:
    lib = _load(label)
    q = query.strip().lower()
    hits = [t for t in lib.terms if q in t.lower()] if q else list(lib.terms)[:limit]
    hits.sort(key=lambda t: (0 if t.lower().startswith(q) else 1, len(t), t))
    return [{"term": t, "n_genes": len(lib.terms[t])} for t in hits[:limit]]


def genes(label: str, term: str) -> List[str]:
    lib = _load(label)
    if term not in lib.terms:
        raise PathwayError(f"term {term!r} not in {label!r}")
    return list(lib.terms[term])
