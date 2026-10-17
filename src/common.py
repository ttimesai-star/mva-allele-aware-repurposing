"""Shared helpers: paths, cached HTTP (so every run is reproducible offline after the first fetch)."""
from __future__ import annotations

import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
OUTPUTS = ROOT / "outputs"
FIGURES = ROOT / "figures"
for _d in (RAW, PROCESSED, OUTPUTS, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

UA = "mva-rx-allele-aware/1.0 (hackathon research; contact via GitHub issues)"


def load_yaml(name: str):
    with open(CONFIG / name, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _cache_path(url: str, body: bytes | None) -> Path:
    key = hashlib.sha1(url.encode() + (body or b"")).hexdigest()[:20]
    return RAW / "http_cache" / f"{key}.json"


def fetch_json(url: str, body: dict | None = None, retries: int = 4, pause: float = 0.35, refresh: bool = False):
    """GET (or POST JSON) with an on-disk cache under data/raw/http_cache.

    The cache is committed to the repository, so the full pipeline reruns
    byte-identically without network access (`make offline`).
    """
    data = json.dumps(body).encode() if body is not None else None
    cp = _cache_path(url, data)
    if cp.exists() and not refresh:
        return json.loads(cp.read_text(encoding="utf-8"))["payload"]
    headers = {"User-Agent": UA, "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=60) as r:
                payload = json.loads(r.read().decode("utf-8"))
            cp.parent.mkdir(parents=True, exist_ok=True)
            cp.write_text(json.dumps({"url": url, "body": body, "fetched": time.strftime("%Y-%m-%d"), "payload": payload},
                                     ensure_ascii=False), encoding="utf-8")
            time.sleep(pause)
            return payload
        except Exception as exc:  # network hiccups: back off and retry
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"fetch failed: {url}: {last}")


def eutils(tool: str, **params) -> dict:
    params.setdefault("retmode", "json")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{tool}.fcgi?" + urllib.parse.urlencode(params)
    return fetch_json(url)
