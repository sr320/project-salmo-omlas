"""Shared ingestion helpers."""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from salmo_omlas.config import DATA_RAW, ensure_dirs


def now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def save_raw(name: str, content: str | bytes) -> Path:
    ensure_dirs()
    path = DATA_RAW / name
    if isinstance(content, str):
        path.write_text(content, encoding="utf-8")
    else:
        path.write_bytes(content)
    return path


def http_get(url: str, **kwargs: Any) -> requests.Response:
    r = requests.get(url, timeout=kwargs.pop("timeout", 120), **kwargs)
    r.raise_for_status()
    return r


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, na_filter=False)


def dumps_metadata(obj: dict[str, Any]) -> str:
    return json.dumps(obj, ensure_ascii=False)
