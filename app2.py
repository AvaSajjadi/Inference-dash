# app2.py — Inference Dash research-tool UI
# Finalized version with:
# - signature upload + column mapping + normalization preview
# - CIE / ORNOR execution
# - progress tracking
# - tabbed results
# - TF/regulator tables + plots
# - pathway enrichment table + plot
# - network graph shown directly in results
# - local Dash/Plotly asset serving
# - reproducible downloads
# - logs / files panel
# - readable regulator / target names in tables and graph
# - improved Dash visualization
# - TF names + IDs
# - better network readability for CIE and ORNOR
# - result display controls:
#   * Top N regulators
#   * minimum regulator score
#   * minimum absolute edge score
#   * max nodes
#   * max edges
#   * label top hubs
# - FIXED signature preprocessing:
#   * deduplicate by identifier before filtering
#   * prefer lowest p-value, then largest |fc|
#   * prevents duplicate Entrez/gene rows from inflating ORNOR input
# - FIXED stuck-running UI:
#   * if output files already exist and are non-empty,
#     Dash auto-promotes running -> done so results tabs always appear
# - FIXED target label mapping in network:
#   * explicitly builds target_name_map from cleaned edge table
#   * prefers trg_name, then entities lookup, so target labels are readable
# - FIXED ORNOR TF parsing:
#   * robust delimiter detection for result tables
#   * handles tab-separated ORNOR TF output even when filename ends with .csv
#   * uses X as main ORNOR display score and T as secondary confidence column
#   * keeps TF / symbol / X / T as separate columns in the UI
# - FIXED network selection for ORNOR:
#   * three_tissues.rels now correctly points to three_tissues.rels (was wrongly pointing to all_tissues_entrez.rels)
# - FIXED ORNOR sampler (2025-04):
#   * --min_targets 500  (reduces TF count to ~70, enabling convergence)
#   * --log2fc 2.0       (strong evidence only, better signal)
#   * --chains 3         (proper multi-chain GR convergence)
#   * --pval 1.0         (sw780-style files have pval=1; filter on FC only)

import base64
import csv
import json
import re
import os
import sys
import subprocess
import shutil
import threading
import time
import traceback
import uuid
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, dcc, html, dash_table, no_update

# =============================================================================
# Paths / config
# =============================================================================

ROOT = Path(__file__).resolve().parent
UPLOADS_DIR = ROOT / "uploads"
RESULTS_DIR = ROOT / "results"
JOBS_DIR = RESULTS_DIR / "_jobs"
NETWORKS_DIR = ROOT / "networks"

UPLOADS_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
JOBS_DIR.mkdir(parents=True, exist_ok=True)

CIE_RUNNER = ROOT / "backend" / "engines" / "cie" / "run_cie.R"
ORNOR_RUNNER = ROOT / "backend" / "engines" / "ornor" / "run_ornor_real.py"

PROGRESS_RE = re.compile(r"PROGRESS:\s*(\d{1,3})\b", re.IGNORECASE)


def first_existing(*paths: Path) -> str:
    for p in paths:
        if p.exists():
            return str(p)
    return str(paths[0]) if paths else ""


DEFAULT_NETWORKS: Dict[str, Dict[str, str]] = {
    # Professor's CIE default: three_tissue.rels + ChIPfilter.ents
    "tcChIP (three_tissue.rels) [default]": {
        "rels": first_existing(
            NETWORKS_DIR / "three_tissue.rels",
            NETWORKS_DIR / "three_tissues.rels",
        ),
        "entities": first_existing(
            NETWORKS_DIR / "ChIPfilter.ents",
            NETWORKS_DIR / "all_tissues.fixed.entities",
            NETWORKS_DIR / "all_tissues.entities",
        ),
    },

    "tcChIP (all_tissues.rels)": {
        "rels": first_existing(
            NETWORKS_DIR / "all_tissues.rels",
        ),
        "entities": first_existing(
            NETWORKS_DIR / "ChIPfilter.ents",
            NETWORKS_DIR / "all_tissues.fixed.entities",
            NETWORKS_DIR / "all_tissues.entities",
        ),
    },

    "ORNOR Entrez (three_tissues.entrez.rels)": {
        "rels": first_existing(
            NETWORKS_DIR / "three_tissues.entrez.rels",
        ),
        "entities": None,
    },

    "ORNOR Entrez (all_tissues_entrez.rels)": {
        "rels": first_existing(
            NETWORKS_DIR / "all_tissues_entrez.rels",
        ),
        "entities": None,
    },
}


def preferred_network_choice_for_engine(engine: Optional[str]) -> str:
    if engine == "ORNOR":
        # ORNOR defaults to all_tissues Entrez network (better coverage)
        return "ORNOR Entrez (all_tissues_entrez.rels)"
    else:
        # CIE defaults to professor's setup: three_tissue.rels + ChIPfilter
        return "tcChIP (three_tissue.rels) [default]"


# =============================================================================
# Theme
# =============================================================================

THEME = {
    "page_bg": "#EEF6FF",
    "panel_bg": "#FFFFFF",
    "panel_border": "#CFE3FF",
    "title": "#0B2A4A",
    "text": "#163B5E",
    "muted": "#5F7FA6",
    "blue": "#2F80ED",
    "blue2": "#56CCF2",
    "good": "#12B981",
    "warn": "#F59E0B",
    "bad": "#EF4444",
    "soft": "#F7FBFF",
    "pink": "#F28B82",
    "green2": "#57CC99",
    "purple": "#7C83FD",
}

BTN = {
    "width": "100%",
    "height": "44px",
    "borderRadius": "14px",
    "border": "none",
    "background": f"linear-gradient(90deg, {THEME['blue']} 0%, {THEME['blue2']} 100%)",
    "color": "white",
    "fontWeight": "900",
    "fontSize": "14px",
    "cursor": "pointer",
}
BTN_DISABLED = {**BTN, "opacity": 0.45, "cursor": "not-allowed"}

BTN_SECONDARY = {
    "width": "100%",
    "height": "42px",
    "borderRadius": "14px",
    "border": f"1px solid {THEME['panel_border']}",
    "background": "white",
    "color": THEME["title"],
    "fontWeight": "900",
    "cursor": "pointer",
}

UPLOAD_STYLE = {
    "width": "100%",
    "height": "72px",
    "lineHeight": "72px",
    "borderWidth": "1px",
    "borderStyle": "dashed",
    "borderColor": THEME["panel_border"],
    "borderRadius": "16px",
    "textAlign": "center",
    "cursor": "pointer",
    "background": THEME["soft"],
    "color": THEME["title"],
}

SMALL_UPLOAD_STYLE = {
    **UPLOAD_STYLE,
    "height": "58px",
    "lineHeight": "58px",
}

INPUT_STYLE = {
    "width": "100%",
    "height": "40px",
    "borderRadius": "10px",
    "border": f"1px solid {THEME['panel_border']}",
    "padding": "0 10px",
    "color": THEME["text"],
}


# =============================================================================
# Aliases / parsing rules
# =============================================================================

GENE_ALIASES = {
    "gene", "genes", "symbol", "gene_symbol", "genesymbol", "hgnc_symbol",
    "external_gene_name", "genename", "gene_name", "target"
}

ENTREZ_ALIASES = {
    "entrez", "entrezid", "entrez_id", "ncbi", "ncbi_gene_id", "geneid",
    "gene_id", "entrezgene", "entrez_gene_id", "id"
}

FC_ALIASES = {
    "fc", "foldchange", "fold_change", "logfc", "log2fc", "log_fold_change",
    "log2foldchange", "log2_fold_change", "avg_log2fc", "lfc"
}

PVAL_ALIASES = {
    "pval", "p.value", "p_value", "pvalue", "p", "padj", "adjp", "adj_p",
    "adj.p.val", "adj_p_val", "fdr", "qvalue", "q_value", "qval"
}

DIRECTION_ALIASES = {
    "direction", "dir", "sign", "regulation", "state", "trend"
}


# =============================================================================
# UI helpers
# =============================================================================

def badge(text, color):
    return html.Span(
        text,
        style={
            "display": "inline-block",
            "padding": "6px 10px",
            "borderRadius": "999px",
            "background": color,
            "color": "white",
            "fontWeight": "900",
            "fontSize": "12px",
        },
    )


def card(title, children):
    return html.Div(
        style={
            "background": THEME["panel_bg"],
            "border": f"1px solid {THEME['panel_border']}",
            "borderRadius": "16px",
            "padding": "14px",
            "boxShadow": "0 10px 25px rgba(15, 76, 129, 0.08)",
            "marginBottom": "14px",
        },
        children=[
            html.Div(title, style={"fontWeight": "900", "color": THEME["title"], "marginBottom": "10px"}),
            *children,
        ],
    )


def stat_card(label: str, value: str, color: str = None):
    return html.Div(
        style={
            "background": "#F8FBFF",
            "border": f"1px solid {THEME['panel_border']}",
            "borderRadius": "14px",
            "padding": "12px",
        },
        children=[
            html.Div(label, style={"fontSize": "12px", "fontWeight": "800", "color": THEME["muted"]}),
            html.Div(value, style={"fontSize": "20px", "fontWeight": "900", "color": color or THEME["title"]}),
        ],
    )


def info_row(label: str, value: str, color: Optional[str] = None):
    return html.Div(
        style={
            "display": "flex",
            "justifyContent": "space-between",
            "gap": "12px",
            "padding": "6px 0",
            "borderBottom": "1px solid #EDF4FF",
        },
        children=[
            html.Span(label, style={"color": THEME["muted"], "fontWeight": "700"}),
            html.Span(value, style={"color": color or THEME["text"], "fontWeight": "800", "textAlign": "right"}),
        ],
    )


def clamp_int(v, lo=0, hi=100) -> int:
    try:
        v = int(v)
    except Exception:
        v = lo
    return max(lo, min(hi, v))


def safe_int(v, default: int) -> int:
    try:
        if v is None or v == "":
            return default
        return int(v)
    except Exception:
        return default


def safe_float(v, default: float) -> float:
    try:
        if v is None or v == "":
            return default
        return float(v)
    except Exception:
        return default


def progress_fill_style(pct: int) -> Dict[str, str]:
    pct = clamp_int(pct)
    return {
        "height": "100%",
        "width": f"{pct}%",
        "borderRadius": "999px",
        "background": f"linear-gradient(90deg, {THEME['blue']} 0%, {THEME['blue2']} 100%)",
        "transition": "width 0.25s ease",
    }


def df_table(
    df: pd.DataFrame,
    page_size: int = 15,
    export_name: Optional[str] = None,
) -> dash_table.DataTable:
    table_kwargs = {
        "columns": [{"name": c, "id": c} for c in df.columns],
        "data": df.to_dict("records"),
        "page_size": page_size,
        "sort_action": "native",
        "filter_action": "native",
        "export_format": "csv",
        "style_table": {"overflowX": "auto"},
        "style_cell": {
            "fontFamily": "ui-monospace, Menlo, Monaco, Consolas, monospace",
            "fontSize": 12,
            "padding": "8px",
            "textAlign": "left",
            "whiteSpace": "normal",
            "height": "auto",
            "minWidth": "80px",
            "maxWidth": "260px",
        },
        "style_header": {
            "fontWeight": "800",
            "backgroundColor": "#EAF3FF",
        },
        "style_data_conditional": [
            {"if": {"row_index": "odd"}, "backgroundColor": "#FCFEFF"},
        ],
    }

    if export_name:
        table_kwargs["id"] = export_name

    return dash_table.DataTable(**table_kwargs)


# =============================================================================
# File helpers
# =============================================================================

def safe_name(filename: str) -> str:
    filename = (filename or "").strip().replace("\\", "/").split("/")[-1]
    filename = re.sub(r"[^A-Za-z0-9._-]+", "_", filename)
    return filename[:220] if filename else f"upload_{uuid.uuid4().hex}.dat"


def save_upload(contents_b64: str, filename: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / safe_name(filename)
    _, b64data = contents_b64.split(",", 1)
    path.write_bytes(base64.b64decode(b64data))
    return path


def zip_outputs(out_dir: Path, zip_path: Path, exclude_paths: Optional[Set[Path]] = None) -> None:
    zip_path.parent.mkdir(parents=True, exist_ok=True)

    exclude_resolved = set()
    for p in (exclude_paths or set()):
        try:
            exclude_resolved.add(p.resolve())
        except Exception:
            pass

    try:
        zip_self = zip_path.resolve()
    except Exception:
        zip_self = zip_path

    if zip_path.exists():
        zip_path.unlink()

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        if out_dir.exists():
            # Only include TF results files (*_tfs.tsv or *_tfs.csv)
            # This gives users just the summary results table, not all intermediate files
            for fp in out_dir.glob("*_tfs.tsv"):
                try:
                    resolved = fp.resolve()
                except Exception:
                    resolved = fp
                if resolved != zip_self and resolved not in exclude_resolved:
                    z.write(fp, arcname=fp.name)

            for fp in out_dir.glob("*_tfs.csv"):
                try:
                    resolved = fp.resolve()
                except Exception:
                    resolved = fp
                if resolved != zip_self and resolved not in exclude_resolved:
                    z.write(fp, arcname=fp.name)


# =============================================================================
# Parsing helpers
# =============================================================================

def strip_bom_text(s: str) -> str:
    return str(s).replace("\ufeff", "").strip()


def normalize_colname(name: str) -> str:
    s = strip_bom_text(name).lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def infer_separator(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="ignore")[:4000]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", "|"])
        return dialect.delimiter
    except Exception:
        if "\t" in sample:
            return "\t"
        if ";" in sample and sample.count(";") > sample.count(","):
            return ";"
        if "|" in sample and sample.count("|") > sample.count(","):
            return "|"
        return ","


def clean_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    original_cols = [str(c) for c in df.columns]
    cleaned_cols = [strip_bom_text(c) for c in original_cols]

    seen = {}
    unique_cols = []
    for c in cleaned_cols:
        if c not in seen:
            seen[c] = 0
            unique_cols.append(c)
        else:
            seen[c] += 1
            unique_cols.append(f"{c}__dup{seen[c]}")

    df.columns = unique_cols
    return df


def read_signature_file(path: Path) -> Tuple[pd.DataFrame, str]:
    sep = infer_separator(path)
    df = pd.read_csv(path, sep=sep, engine="python", encoding="utf-8-sig")
    df = df.loc[:, ~df.columns.astype(str).str.contains(r"^Unnamed", case=False, regex=True)]
    df = clean_dataframe_columns(df)
    return df, sep


def guess_role(columns: List[str], aliases: set) -> Optional[str]:
    exact_norm = {normalize_colname(c): c for c in columns}
    for a in aliases:
        if a in exact_norm:
            return exact_norm[a]

    for c in columns:
        nc = normalize_colname(c)
        if nc in aliases:
            return c

    for c in columns:
        nc = normalize_colname(c)
        if aliases is FC_ALIASES and (nc.startswith("log2foldchange") or nc.startswith("logfc")):
            return c
        if aliases is PVAL_ALIASES and (nc.startswith("padj") or nc.startswith("adj_p") or nc.startswith("pval")):
            return c
    return None


def parse_direction_value(x) -> int:
    if pd.isna(x):
        return 0
    s = strip_bom_text(x).lower()
    if s in {"1", "+1", "up", "pos", "positive", "increase", "increased"}:
        return 1
    if s in {"-1", "down", "neg", "negative", "decrease", "decreased"}:
        return -1
    try:
        f = float(s)
        if f > 0:
            return 1
        if f < 0:
            return -1
    except Exception:
        pass
    return 0


def sign_from_fc(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce")
    out = pd.Series(np.where(vals > 0, 1, np.where(vals < 0, -1, 0)), index=series.index)
    return out.astype(int)


def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def summarize_preview(df: pd.DataFrame, max_rows: int = 20) -> List[Dict]:
    out = df.head(max_rows).copy().replace({np.nan: ""})
    for c in out.columns:
        out[c] = out[c].astype(str)
    return out.to_dict("records")


def make_dropdown_options(columns: List[str]) -> List[Dict[str, str]]:
    opts = [{"label": "— None —", "value": "__NONE__"}]
    opts.extend([{"label": c, "value": c} for c in columns])
    return opts


def human_bool(v: bool) -> str:
    return "Yes" if bool(v) else "No"


def _make_dedup_key(entrez_series: pd.Series, gene_series: pd.Series) -> pd.Series:
    entrez_txt = entrez_series.copy()
    gene_txt = gene_series.copy()

    entrez_txt = entrez_txt.fillna("").astype(str).str.strip()
    gene_txt = gene_txt.fillna("").astype(str).str.strip()

    def _norm_one(x: str) -> str:
        if not x:
            return ""
        try:
            f = float(x)
            if np.isfinite(f) and float(f).is_integer():
                return str(int(f))
        except Exception:
            pass
        return x

    entrez_txt = entrez_txt.map(_norm_one)
    gene_txt = gene_txt.map(lambda x: x if x.lower() not in {"nan", "none"} else "")

    key = entrez_txt.copy()
    use_gene = key.eq("")
    key.loc[use_gene] = gene_txt.loc[use_gene]

    return key.fillna("").astype(str).str.strip()


def deduplicate_signature_rows(norm_df: pd.DataFrame) -> pd.DataFrame:
    if norm_df is None or norm_df.empty:
        return norm_df

    df = norm_df.copy()

    if "gene" not in df.columns:
        df["gene"] = ""
    if "entrez" not in df.columns:
        df["entrez"] = np.nan
    if "fc" not in df.columns:
        df["fc"] = np.nan
    if "pval" not in df.columns:
        df["pval"] = 1.0
    if "sign" not in df.columns:
        df["sign"] = 0

    df["_dedup_key"] = _make_dedup_key(df["entrez"], df["gene"])
    df = df[df["_dedup_key"].astype(str).str.len() > 0].copy()

    if df.empty:
        return df

    df["_pval_sort"] = pd.to_numeric(df["pval"], errors="coerce").fillna(np.inf)
    df["_abs_fc_sort"] = pd.to_numeric(df["fc"], errors="coerce").abs().fillna(-np.inf)
    df["_nonzero_sign_sort"] = (pd.to_numeric(df["sign"], errors="coerce").fillna(0).abs() > 0).astype(int)

    df = df.sort_values(
        by=["_dedup_key", "_pval_sort", "_abs_fc_sort", "_nonzero_sign_sort"],
        ascending=[True, True, False, False],
        kind="mergesort",
    )

    df = df.drop_duplicates(subset=["_dedup_key"], keep="first").copy()
    df = df.drop(columns=["_dedup_key", "_pval_sort", "_abs_fc_sort", "_nonzero_sign_sort"], errors="ignore")
    df = df.reset_index(drop=True)

    return df


# =============================================================================
# Name mapping helpers
# =============================================================================

def _best_existing_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    lower_map = {str(c).strip().lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def _smart_read_table(path: Path) -> pd.DataFrame:
    if not path or not path.exists():
        return pd.DataFrame()

    sep = None
    try:
        sample = path.read_text(encoding="utf-8", errors="ignore")[:5000]
        if "\t" in sample and sample.count("\t") >= sample.count(","):
            sep = "\t"
        else:
            sep = ","
    except Exception:
        sep = ","

    try:
        df = pd.read_csv(path, sep=sep, engine="python", encoding="utf-8-sig")
        df = clean_dataframe_columns(df)
        return df
    except Exception:
        try:
            df = pd.read_csv(path, sep=None, engine="python", encoding="utf-8-sig")
            df = clean_dataframe_columns(df)
            return df
        except Exception:
            return pd.DataFrame()


def _normalize_mapping_key(uid) -> str:
    s = str(uid).strip()
    if not s or s.lower() in {"nan", "none"}:
        return ""
    try:
        f = float(s)
        if f.is_integer():
            return str(int(f))
    except Exception:
        pass
    return s


def _extract_uid_name_map_from_df(df: pd.DataFrame) -> Dict[str, str]:
    if df is None or df.empty:
        return {}

    uid_col = _best_existing_col(df, ["uid", "id", "entrez", "entrezid", "gene_id", "ncbi_gene_id"])
    name_col = _best_existing_col(df, ["name", "symbol", "gene", "gene_symbol", "label", "gene_name"])

    if not uid_col or not name_col:
        return {}

    out = {}
    for _, row in df[[uid_col, name_col]].dropna().iterrows():
        uid = _normalize_mapping_key(row[uid_col])
        name = str(row[name_col]).strip()
        if not uid or not name or name.lower() in {"nan", "none"}:
            continue
        out[uid] = name
    return out


_entities_name_map_cache: Dict[str, Dict[str, str]] = {}
_entities_name_map_lock = __import__("threading").Lock()

def load_entities_name_map(entities_path: Optional[str]) -> Dict[str, str]:
    cache_key = str(entities_path or "")
    if cache_key in _entities_name_map_cache:
        return _entities_name_map_cache[cache_key]

    with _entities_name_map_lock:
        # double-checked locking: re-check inside lock
        if cache_key in _entities_name_map_cache:
            return _entities_name_map_cache[cache_key]

        # all expensive I/O happens inside the lock so concurrent callers
        # block here instead of each doing 10 s of duplicate work
        out: Dict[str, str] = {}

        if entities_path:
            p = Path(entities_path)
            if p.exists():
                df = _smart_read_table(p)
                out.update(_extract_uid_name_map_from_df(df))

        candidate_dirs = [ROOT / "uploads", ROOT / "networks", ROOT / "resources"]
        candidate_names = [
            "tcCHIP_uid_map.tsv",
            "entrez_to_symbol.tsv",
            "entrez_to_symbol.csv",
            "entrez_symbol.tsv",
            "entrez_symbol.csv",
            "gene_map.tsv",
            "gene_map.csv",
        ]

        for d in candidate_dirs:
            if not d.exists():
                continue

            paths = []
            for name in candidate_names:
                fp = d / name
                if fp.exists():
                    paths.append(fp)

            if not paths:
                paths.extend(sorted(d.glob("*entrez*symbol*.tsv")))
                paths.extend(sorted(d.glob("*entrez*symbol*.csv")))
                paths.extend(sorted(d.glob("*symbol*entrez*.tsv")))
                paths.extend(sorted(d.glob("*symbol*entrez*.csv")))

            for fp in paths:
                try:
                    df = _smart_read_table(fp)
                    extra = _extract_uid_name_map_from_df(df)
                    for k, v in extra.items():
                        out.setdefault(k, v)
                except Exception:
                    pass

        _entities_name_map_cache[cache_key] = out
        return out


def lookup_name(uid, name_map: Dict[str, str]) -> str:
    s = _normalize_mapping_key(uid)
    if not s:
        return ""
    if s in name_map:
        return str(name_map[s]).strip()
    return s


def annotate_tf_table_with_names(tf_df: pd.DataFrame, entities_name_map: Dict[str, str]) -> pd.DataFrame:
    if tf_df is None or tf_df.empty:
        return tf_df

    df = tf_df.copy()
    df.columns = [strip_bom_text(c) for c in df.columns]

    id_col = _best_existing_col(df, ["TF", "source", "id", "uid", "tf_id", "tf"])
    symbol_col = _best_existing_col(df, ["symbol", "display_name", "name", "tf"])

    def _is_missing_text(x: str) -> bool:
        s = str(x).strip()
        return s == "" or s.lower() in {"nan", "none"}

    def _looks_like_numeric_id(x: str) -> bool:
        s = str(x).strip()
        if _is_missing_text(s):
            return True
        try:
            return float(s).is_integer()
        except Exception:
            return False

    if id_col:
        raw_ids = df[id_col].astype(str).str.strip()
        mapped_ids = raw_ids.map(lambda x: lookup_name(x, entities_name_map))
    else:
        raw_ids = pd.Series([""] * len(df), index=df.index)
        mapped_ids = pd.Series([""] * len(df), index=df.index)

    if symbol_col:
        original_symbols = df[symbol_col].astype(str).str.strip()
    else:
        original_symbols = pd.Series([""] * len(df), index=df.index)

    display_names = []
    for idx in df.index:
        uid = raw_ids.loc[idx] if idx in raw_ids.index else ""
        mapped = mapped_ids.loc[idx] if idx in mapped_ids.index else ""
        original = original_symbols.loc[idx] if idx in original_symbols.index else ""

        if _is_missing_text(original):
            chosen = mapped if mapped else uid
        elif uid and original == uid and mapped and mapped != uid:
            chosen = mapped
        elif _looks_like_numeric_id(original) and mapped and mapped != uid:
            chosen = mapped
        else:
            chosen = original

        if _is_missing_text(chosen):
            chosen = mapped if mapped else uid

        display_names.append(str(chosen).strip())

    df["display_name"] = pd.Series(display_names, index=df.index).astype(str)

    if id_col:
        df["display_label"] = np.where(
            df["display_name"].astype(str).str.strip() == raw_ids.astype(str).str.strip(),
            raw_ids.astype(str),
            df["display_name"].astype(str) + " (" + raw_ids.astype(str) + ")",
        )
    else:
        df["display_label"] = df["display_name"].astype(str)

    if "X" in df.columns:
        df["display_score"] = pd.to_numeric(df["X"], errors="coerce").fillna(0.0)
    elif "T" in df.columns:
        df["display_score"] = pd.to_numeric(df["T"], errors="coerce").fillna(0.0)
    elif "display_score" not in df.columns:
        df["display_score"] = 0.0

    preferred = [c for c in ["display_name", "display_label", "display_score", "TF", "symbol", "X", "T"] if c in df.columns]
    remaining = [c for c in df.columns if c not in preferred]
    return df[preferred + remaining]


def annotate_edge_table_with_names(edge_df: pd.DataFrame, entities_name_map: Dict[str, str]) -> pd.DataFrame:
    if edge_df is None or edge_df.empty:
        return edge_df

    df = edge_df.copy()
    df.columns = [strip_bom_text(c) for c in df.columns]
    lower = {c.lower(): c for c in df.columns}

    src_col = lower.get("srcuid") or lower.get("source") or lower.get("src")
    trg_col = lower.get("trguid") or lower.get("target") or lower.get("trg")

    if src_col:
        df[src_col] = df[src_col].astype(str).str.strip()
        df = df[df[src_col].str.len() > 0].copy()
        src_lookup = {uid: lookup_name(uid, entities_name_map) for uid in df[src_col].unique()}
        df["src_name"] = df[src_col].map(src_lookup)
        df["src_label"] = np.where(
            df["src_name"].astype(str).str.strip() == df[src_col].astype(str).str.strip(),
            df[src_col].astype(str),
            df["src_name"].astype(str) + " (" + df[src_col].astype(str) + ")",
        )

    if trg_col:
        df[trg_col] = df[trg_col].astype(str).str.strip()
        df = df[df[trg_col].str.len() > 0].copy()
        trg_lookup = {uid: lookup_name(uid, entities_name_map) for uid in df[trg_col].unique()}
        df["trg_name"] = df[trg_col].map(trg_lookup)
        df["trg_label"] = np.where(
            df["trg_name"].astype(str).str.strip() == df[trg_col].astype(str).str.strip(),
            df[trg_col].astype(str),
            df["trg_name"].astype(str) + " (" + df[trg_col].astype(str) + ")",
        )

    ordered = [c for c in ["src_name", "trg_name", "src_label", "trg_label"] if c in df.columns]
    remaining = [c for c in df.columns if c not in ordered]
    return df[ordered + remaining]


def build_tf_name_map(tf_df: pd.DataFrame, entities_name_map: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    entities_name_map = entities_name_map or {}
    out = {}

    if tf_df is not None and not tf_df.empty:
        df = tf_df.copy()
        df.columns = [strip_bom_text(c) for c in df.columns]

        id_col = _best_existing_col(df, ["TF", "source", "id", "uid", "tf_id"])
        name_col = _best_existing_col(df, ["display_name", "symbol", "name", "tf"])

        if id_col and name_col:
            for _, row in df[[id_col, name_col]].dropna().iterrows():
                uid = str(row[id_col]).strip()
                name = str(row[name_col]).strip()
                if uid and name and uid.lower() not in {"nan", "none"} and name.lower() not in {"nan", "none"}:
                    out[uid] = name

    for uid, name in entities_name_map.items():
        out.setdefault(str(uid), str(name))

    return out


# =============================================================================
# Job status
# =============================================================================

@dataclass
class JobPaths:
    job_id: str
    job_dir: Path
    status_json: Path
    stdout_log: Path
    out_dir: Path
    zip_path: Path


def _job_paths(job_id: str) -> JobPaths:
    job_dir = JOBS_DIR / job_id
    out_dir = job_dir / "out"
    return JobPaths(
        job_id=job_id,
        job_dir=job_dir,
        status_json=job_dir / "status.json",
        stdout_log=job_dir / "stdout.log",
        out_dir=out_dir,
        zip_path=job_dir / "outputs.zip",
    )


def write_status(p: JobPaths, **kwargs):
    p.job_dir.mkdir(parents=True, exist_ok=True)
    p.out_dir.mkdir(parents=True, exist_ok=True)
    base = {}
    if p.status_json.exists():
        try:
            base = json.loads(p.status_json.read_text())
        except Exception:
            base = {}
    base.update(kwargs)
    base.setdefault("job_id", p.job_id)
    base.setdefault("state", "queued")
    base.setdefault("progress", 0)
    base.setdefault("message", "")
    base["progress"] = clamp_int(base.get("progress", 0))
    base["updated_at"] = time.time()
    p.status_json.write_text(json.dumps(base, indent=2))


def read_status(p: JobPaths) -> Dict:
    if not p.status_json.exists():
        return {"job_id": p.job_id, "state": "missing", "progress": 0, "message": ""}
    try:
        d = json.loads(p.status_json.read_text())
        d["progress"] = clamp_int(d.get("progress", 0))
        d.setdefault("state", "unknown")
        d.setdefault("message", "")
        return d
    except Exception:
        return {"job_id": p.job_id, "state": "corrupt", "progress": 0, "message": "Unreadable status.json"}


# =============================================================================
# Signature metadata / normalization
# =============================================================================

def build_signature_metadata(raw_path: Path) -> Dict:
    df, sep = read_signature_file(raw_path)
    cols = [str(c) for c in df.columns]

    meta = {
        "raw_path": str(raw_path),
        "separator": "\\t" if sep == "\t" else sep,
        "rows": int(len(df)),
        "cols": int(df.shape[1]),
        "columns": cols,
        "preview_records": summarize_preview(df, 20),
        "detected": {
            "gene_col": guess_role(cols, GENE_ALIASES) or "__NONE__",
            "entrez_col": guess_role(cols, ENTREZ_ALIASES) or "__NONE__",
            "fc_col": guess_role(cols, FC_ALIASES) or "__NONE__",
            "pval_col": guess_role(cols, PVAL_ALIASES) or "__NONE__",
            "direction_col": guess_role(cols, DIRECTION_ALIASES) or "__NONE__",
        },
    }
    return meta


def normalize_signature_for_engine(
    raw_path: Path,
    engine: str,
    out_dir: Path,
    gene_col: Optional[str],
    entrez_col: Optional[str],
    fc_col: Optional[str],
    pval_col: Optional[str],
    direction_col: Optional[str],
    pval_thresh: float,
    abs_fc_thresh: Optional[float],
    use_pval_filter: bool,
    use_abs_fc_filter: bool,
) -> Dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    df, sep = read_signature_file(raw_path)

    def noneify(x):
        return None if x in (None, "", "__NONE__") else x

    gene_col = noneify(gene_col)
    entrez_col = noneify(entrez_col)
    fc_col = noneify(fc_col)
    pval_col = noneify(pval_col)
    direction_col = noneify(direction_col)

    # Case-insensitive column resolver — find the actual column name in df
    def resolve_col(col):
        if col is None:
            return None
        if col in df.columns:
            return col
        col_lower = col.lower()
        for c in df.columns:
            if c.lower() == col_lower:
                return c
        return col  # return as-is, will fail later with a clear error

    errors = []
    warnings = []

    if engine == "CIE":
        if not entrez_col and not gene_col:
            errors.append("CIE requires an identifier column.")
        if not fc_col:
            errors.append("CIE requires a fold-change column.")
        if not pval_col:
            errors.append("CIE requires a p-value / adjusted p-value column.")
    else:
        if not entrez_col and not gene_col:
            errors.append("ORNOR requires an identifier column.")
        if not direction_col and not fc_col:
            errors.append("ORNOR requires either a direction column or a fold-change column.")

    if errors:
        return {"ok": False, "errors": errors, "warnings": warnings}

    if engine == "ORNOR" and use_pval_filter and not pval_col:
        errors.append(
            "ORNOR: p-value filtering is enabled but no p-value column is mapped. "
            "All genes would receive pval=1.0 and none would pass the p <= threshold filter. "
            "Either map a p-value column or uncheck 'Use p-value filtering'."
        )
        return {"ok": False, "errors": errors, "warnings": warnings}

    norm = pd.DataFrame(index=df.index)

    gene_col = resolve_col(gene_col)
    entrez_col = resolve_col(entrez_col)
    fc_col = resolve_col(fc_col)
    pval_col = resolve_col(pval_col)
    direction_col = resolve_col(direction_col)

    if gene_col:
        norm["gene"] = df[gene_col].astype(str).str.strip()
        norm.loc[norm["gene"].str.lower().isin({"nan", "none", ""}), "gene"] = ""
    else:
        norm["gene"] = ""

    if entrez_col:
        norm["entrez"] = pd.to_numeric(df[entrez_col], errors="coerce")
    elif gene_col:
        # Auto-convert gene symbols to Entrez IDs
        print(f"[ORNOR] No entrez column mapped — attempting symbol->Entrez conversion")
        sym_map = _symbol_to_entrez_map()
        norm["entrez"] = df[gene_col].astype(str).str.strip().str.upper().map(sym_map)
        n_mapped = norm["entrez"].notna().sum()
        print(f"[ORNOR] Converted {n_mapped}/{len(norm)} gene symbols to Entrez IDs")
    else:
        norm["entrez"] = np.nan

    if fc_col:
        norm["fc"] = coerce_numeric(df[fc_col])
    else:
        norm["fc"] = np.nan

    if pval_col:
        norm["pval"] = coerce_numeric(df[pval_col]).clip(lower=1e-300)
    else:
        norm["pval"] = 1.0
        if engine == "ORNOR":
            warnings.append(
                "No p-value column provided; all genes will receive pval=1.0. "
                "If p-value filtering is enabled, NO genes will pass the evidence threshold and ORNOR will fail. "
                "Either map a p-value column, or disable p-value filtering."
            )

    if direction_col:
        norm["direction"] = df[direction_col].apply(parse_direction_value).astype(int)
        sign_source = f"direction column '{direction_col}'"
    else:
        norm["direction"] = 0
        sign_source = "none"

    if direction_col:
        norm["sign"] = norm["direction"].astype(int)
    elif fc_col:
        norm["sign"] = sign_from_fc(norm["fc"])
        sign_source = f"inferred from fold-change column '{fc_col}'"
    else:
        norm["sign"] = 0

    before_dedup_rows = int(len(norm))
    norm = deduplicate_signature_rows(norm)
    after_dedup_rows = int(len(norm))

    include = pd.Series(True, index=norm.index)

    if use_pval_filter and pval_col and engine != "CIE":
        include = include & norm["pval"].notna() & (norm["pval"] <= float(pval_thresh))

    if use_abs_fc_filter and fc_col and abs_fc_thresh is not None and engine != "CIE":
        include = include & norm["fc"].notna() & (norm["fc"].abs() >= float(abs_fc_thresh))

    valid_id = pd.Series(False, index=norm.index)
    if entrez_col:
        valid_id = valid_id | norm["entrez"].notna()
    if gene_col:
        valid_id = valid_id | norm["gene"].astype(str).str.len().gt(0)
    include = include & valid_id

    norm["include"] = include
    filtered = norm[norm["include"]].copy()

    if filtered.empty:
        return {
            "ok": False,
            "errors": ["No rows remain after deduplication, identifier checks, and selected filters."],
            "warnings": warnings,
        }

    canonical_path = out_dir / "canonical_signature.tsv"
    norm.to_csv(canonical_path, sep="\t", index=False)

    cie_input_path = out_dir / "cie_input.tsv"
    if engine == "CIE":
        # Send full signature to CIE — let CIE apply pval/FC filtering internally
        cie_df = norm.copy()
        cie_df["entrez"] = pd.to_numeric(cie_df["entrez"], errors="coerce")
        cie_df = cie_df[cie_df["entrez"].notna() & cie_df["fc"].notna()].copy()
        cie_df["entrez"] = cie_df["entrez"].astype("Int64")
        print(f"[CIE] sending {len(cie_df)} genes to CIE (unfiltered, CIE handles thresholds internally)")
        cie_write = cie_df[["entrez", "fc", "pval"]].copy()
        cie_write.columns = ["entrez", "fc", "pval"]
        cie_write.to_csv(cie_input_path, sep="\t", index=False)

    ornor_input_path = out_dir / "ornor_input.tsv"
    if engine == "ORNOR":
        ornor_df = filtered.copy()

        ornor_df["entrez"] = pd.to_numeric(ornor_df.get("entrez"), errors="coerce")
        before_ornor_rows = int(len(ornor_df))
        ornor_df = ornor_df[ornor_df["entrez"].notna() & (ornor_df["entrez"] != 0)].copy()
        dropped_bad_entrez = before_ornor_rows - int(len(ornor_df))
        if dropped_bad_entrez > 0:
            warnings.append(
                f"Removed {dropped_bad_entrez} ORNOR rows with missing/zero Entrez IDs because the selected ORNOR network uses Entrez targets."
            )

        if not ornor_df.empty:
            ornor_df["entrez"] = ornor_df["entrez"].astype("Int64")
            ornor_df = ornor_df.rename(columns={"fc": "logfc"})

        ornor_write = ornor_df[["entrez", "logfc", "pval"]].copy()
        ornor_write = ornor_write.reset_index(drop=True)
        ornor_write.to_csv(ornor_input_path, sep="\t", index=False)

    stats = {
        "rows_loaded": int(len(df)),
        "rows_after_dedup": after_dedup_rows,
        "duplicates_removed": max(0, before_dedup_rows - after_dedup_rows),
        "rows_after_filter": int(len(filtered)),
        "separator": "\\t" if sep == "\t" else sep,
        "gene_col": gene_col or "",
        "entrez_col": entrez_col or "",
        "fc_col": fc_col or "",
        "pval_col": pval_col or "",
        "direction_col": direction_col or "",
        "sign_source": sign_source,
        "pval_filter_used": bool(use_pval_filter and pval_col),
        "abs_fc_filter_used": bool(use_abs_fc_filter and fc_col and abs_fc_thresh is not None),
        "positive_signs": int((filtered["sign"] > 0).sum()),
        "negative_signs": int((filtered["sign"] < 0).sum()),
        "zero_signs": int((filtered["sign"] == 0).sum()),
    }

    (out_dir / "normalization_summary.json").write_text(json.dumps(stats, indent=2))

    preview_norm = filtered.head(20).copy().replace({np.nan: ""})
    for c in preview_norm.columns:
        preview_norm[c] = preview_norm[c].astype(str)

    return {
        "ok": True,
        "errors": [],
        "warnings": warnings,
        "canonical_path": str(canonical_path),
        "engine_input_path": str(cie_input_path if engine == "CIE" else ornor_input_path),
        "stats": stats,
        "normalized_preview_records": preview_norm.to_dict("records"),
    }


# =============================================================================
# Runner helpers
# =============================================================================

def run_with_progress(cmd: list, job: JobPaths, cwd: Optional[Path] = None, extra_env: Optional[Dict[str, str]] = None) -> Tuple[int, str]:
    write_status(job, state="running", progress=0, message="Starting...", cmd=cmd)

    job.stdout_log.parent.mkdir(parents=True, exist_ok=True)
    last_msg = "Running..."

    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)

    with job.stdout_log.open("w", encoding="utf-8") as logf:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=str(cwd) if cwd else None,
            env=env,
        )
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                logf.write(line)
                logf.flush()
                s = line.strip()
                if not s:
                    continue
                m = PROGRESS_RE.search(s)
                if m:
                    pct = clamp_int(m.group(1))
                    write_status(job, progress=pct, message=f"Running... ({pct}%)")
                else:
                    last_msg = s[:240]
                    write_status(job, message=last_msg)
            rc = proc.wait()
            return rc, last_msg
        except Exception as e:
            try:
                proc.kill()
            except Exception:
                pass
            return 1, f"Runner crashed: {e}"




def _symbol_to_entrez_map() -> dict:
    """Load gene symbol -> Entrez ID mapping from resources."""
    try:
        import pandas as pd
        p = ROOT / "resources" / "entrez_to_symbol.tsv"
        df = pd.read_csv(p, sep="\t")
        df.columns = [c.strip().upper() for c in df.columns]
        return dict(zip(df["SYMBOL"].str.strip().str.upper(), df["ENTREZID"].astype(int)))
    except Exception as e:
        print(f"[symbol->entrez] failed to load map: {e}")
        return {}

def _auto_tune_min_targets(ornor_input_path: str, rels_path: str, target_tfs: int = 25) -> int:
    """Pick --min_targets so we get ~target_tfs TFs with enough evidence overlap."""
    try:
        import pandas as pd
        sig = pd.read_csv(ornor_input_path, sep="\t")
        net = pd.read_csv(rels_path, sep="\t")
        net.columns = [c.lower().strip() for c in net.columns]
        trg_col = next((c for c in net.columns if c in ("trguid", "trg", "target")), None)
        src_col = next((c for c in net.columns if c in ("srcuid", "src", "source")), None)
        if not trg_col or not src_col:
            return 25
        sig_ids = set(pd.to_numeric(sig["entrez"], errors="coerce").dropna().astype(int))
        counts = []
        for tf, grp in net.groupby(src_col):
            tf_targets = set(pd.to_numeric(grp[trg_col], errors="coerce").dropna().astype(int))
            counts.append(len(tf_targets & sig_ids))
        counts.sort(reverse=True)
        if not counts or counts[0] == 0:
            return 25
        # Find threshold that gives ~target_tfs TFs
        for threshold in range(1, max(counts) + 1):
            n_tfs = sum(1 for c in counts if c >= threshold)
            if n_tfs <= target_tfs:
                return max(1, threshold - 1)
        return max(1, max(counts))
    except Exception as e:
        print(f"[ORNOR] auto_tune_min_targets failed: {e}, using default 25")
        return 25

def job_thread(
    engine: str,
    normalized_input_path: str,
    network_choice: str,
    rels_uploaded: Optional[str],
    ents_uploaded: Optional[str],
    normalization_stats: Dict,
    job_id: str,
    max_edges: int = 100,
    pval_thresh: float = 0.05,
    log2fc_thresh: float = 0.5,
    skip_mcmc: bool = False,
    original_filename: str = "",
):
    job = _job_paths(job_id)
    job.job_dir.mkdir(parents=True, exist_ok=True)
    job.out_dir.mkdir(parents=True, exist_ok=True)

    try:
        try:
            max_edges = int(max_edges) if max_edges is not None else 100
        except Exception:
            max_edges = 100

        if network_choice == "__UPLOAD__":
            rels_path = rels_uploaded
            ents_path = ents_uploaded
        else:
            resolved_choice = network_choice

            if resolved_choice not in DEFAULT_NETWORKS:
                resolved_choice = next(iter(DEFAULT_NETWORKS.keys()), None)

            meta = DEFAULT_NETWORKS.get(resolved_choice, {})
            rels_path = meta.get("rels")
            ents_path = meta.get("entities")

        in_path = Path(normalized_input_path)

        if not in_path.exists():
            write_status(job, state="error", progress=0, message=f"Normalized input file missing on disk: {in_path}")
            return
        if not rels_path or not Path(rels_path).exists():
            write_status(job, state="error", progress=0, message=f"Network .rels missing on disk: {rels_path}")
            return
        if engine == "CIE" and (not ents_path or not Path(ents_path).exists()):
            write_status(job, state="error", progress=0, message=f"Entities file missing on disk: {ents_path}")
            return

        write_status(
            job,
            state="queued",
            progress=0,
            message=f"Queued ({engine})...",
            engine=engine,
            normalization_stats=normalization_stats,
            rels_path=rels_path,
            ents_path=ents_path,
        )

        run_request = {
            "engine": engine,
            "normalized_input_path": str(in_path),
            "network_choice": network_choice,
            "rels_path": str(rels_path) if rels_path is not None else None,
            "ents_path": str(ents_path) if ents_path is not None else None,
            "normalization_stats": normalization_stats,
            "max_edges": max_edges,
            "original_filename": original_filename,
        }
        (job.job_dir / "run_request.json").write_text(
            json.dumps(run_request, indent=2, default=str),
            encoding="utf-8",
        )

        if engine == "CIE":
            out_edges = job.out_dir / "cie_edges.csv"
            # Find Rscript in PATH, or check common container paths
            rscript_path = shutil.which("Rscript")
            if not rscript_path:
                # Fallback paths for Nixpacks and other container environments
                for path in ["/root/.nix-profile/bin/Rscript", "/usr/bin/Rscript", "/nix/var/nix/profiles/default/bin/Rscript"]:
                    if Path(path).exists():
                        rscript_path = path
                        break
            if not rscript_path:
                write_status(job, state="error", message="R (Rscript) not found in PATH or standard locations")
                return
            cmd = [
                rscript_path, str(CIE_RUNNER),
                "-s", str(in_path),
                "-o", str(out_edges),
                "--rels", str(rels_path),
                "--ents", str(ents_path),
                "--db", "tcChIP",
                "--tissue", "all",
                "-m", "Fisher",
                "-p", str(pval_thresh),
                "-f", str(log2fc_thresh),
                "-u", "1",
                "-c", "1",
            ]
            rc, last = run_with_progress(cmd, job)
            if rc != 0:
                write_status(job, state="error", message=f"CIE failed (rc={rc}). {last}")
                return
        else:
            out_tfs = job.out_dir / "ornor_tfs.tsv"
            out_edges = job.out_dir / "ornor_edges.csv"
            ornor_top_edges = 50000
            py_path = os.environ.get("PYTHONPATH", "")
            extra_env = {
                "PYTHONPATH": f"{ROOT}{os.pathsep}{py_path}" if py_path else str(ROOT)
            }
            cmd = [
                sys.executable,
                str(ORNOR_RUNNER),
                str(in_path),
                str(rels_path),
                str(out_edges),
                str(out_tfs),
                str(ornor_top_edges),
                # ── FIXED ORNOR sampler defaults ──────────────────────────────
                # min_targets=500 reduces TF count to ~70 so chains can converge
                # log2fc=2.0 uses only strongly DE genes as evidence
                # chains=3 enables proper multi-chain GR convergence checking
                # pval=1.0 handles files where all pvals are placeholder 1.0
                "--min_targets", "15",
                "--threshold_logic", "and",
                "--pval", str(pval_thresh) if normalization_stats.get("pval_filter_used") else "1.0",
                "--log2fc", str(log2fc_thresh) if normalization_stats.get("abs_fc_filter_used") else "0.0",
                "--chains", "1",
                "--max_samples", "800",
                "--gr_level", "1.1",
            ] + (["--skip_mcmc"] if skip_mcmc else []) \
              + (["--entities", str(ents_path)] if ents_path and Path(ents_path).exists() else [])
            rc, last = run_with_progress(cmd, job, cwd=ROOT, extra_env=extra_env)
            if rc != 0:
                write_status(job, state="error", message=f"ORNOR failed (rc={rc}). {last}")
                return

        zip_outputs(job.job_dir / "out", job.zip_path, exclude_paths={job.zip_path})
        write_status(job, state="done", progress=100, message="Finished.", zip_path=str(job.zip_path))
    except Exception as e:
        tb = traceback.format_exc()
        try:
            job.stdout_log.parent.mkdir(parents=True, exist_ok=True)
            with job.stdout_log.open("a", encoding="utf-8") as logf:
                logf.write("\n=== INTERNAL ERROR ===\n")
                logf.write(tb)
                logf.write("\n")
        except Exception:
            pass
        write_status(job, state="error", progress=0, message=f"{engine} launcher crashed before the runner started: {e}")


# =============================================================================
# Result parsers / display helpers
# =============================================================================

def _find_first_existing(paths: List[Path]) -> Optional[Path]:
    return next((p for p in paths if p.exists()), None)


def read_result_table(fp: Optional[Path], nrows: Optional[int] = None) -> pd.DataFrame:
    if fp is None or not fp.exists():
        return pd.DataFrame()

    attempts = []

    try:
        detected_sep = infer_separator(fp)
        attempts.append({"sep": detected_sep, "engine": "python", "encoding": "utf-8-sig"})
    except Exception:
        pass

    attempts.extend([
        {"sep": None, "engine": "python", "encoding": "utf-8-sig"},
        {"sep": "\t", "engine": "python", "encoding": "utf-8-sig"},
        {"sep": ",", "engine": "python", "encoding": "utf-8-sig"},
        {"sep": ";", "engine": "python", "encoding": "utf-8-sig"},
    ])

    seen = set()
    for kwargs in attempts:
        sig = tuple(sorted(kwargs.items()))
        if sig in seen:
            continue
        seen.add(sig)

        try:
            df = pd.read_csv(fp, nrows=nrows, **kwargs)
            df = df.loc[:, ~df.columns.astype(str).str.contains(r"^Unnamed", case=False, regex=True)]
            df = clean_dataframe_columns(df)

            if df.empty:
                continue

            if df.shape[1] == 1:
                col0 = str(df.columns[0])
                sample_cell = ""
                try:
                    sample_cell = str(df.iloc[0, 0])
                except Exception:
                    sample_cell = ""

                fused = ("\t" in col0) or ("," in col0) or ("\t" in sample_cell and kwargs.get("sep") != "\t")
                if fused:
                    continue

            return df
        except Exception:
            continue

    return pd.DataFrame()


def coerce_first_existing_col(df: pd.DataFrame, names: List[str]) -> Optional[str]:
    lower_to_real = {str(c).strip().lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lower_to_real:
            return lower_to_real[name.lower()]
    return None


def prepare_tf_display_df(
    tf_df: pd.DataFrame,
    engine: str,
    top_n: int,
    min_display_score: float,
    use_pval_filter: bool = True,
) -> Tuple[pd.DataFrame, str]:
    if tf_df is None or tf_df.empty:
        return pd.DataFrame(), "score"

    df = tf_df.copy()
    df.columns = [strip_bom_text(c) for c in df.columns]

    label_col = coerce_first_existing_col(df, ["display_label", "display_name", "symbol", "name", "TF", "tf", "uid", "id", "source"])
    raw_score_label = "score"

    if engine == "CIE":
        p_col = coerce_first_existing_col(df, ["pvalue", "pval"])
        fdr_col = coerce_first_existing_col(df, ["fdr"])
        score_col = p_col or fdr_col
        if score_col:
            vals = pd.to_numeric(df[score_col], errors="coerce").clip(lower=1e-300)
            df["display_score"] = -np.log10(vals)
            raw_score_label = f"-log10({score_col})"
        else:
            df["display_score"] = 0.0
            raw_score_label = "display_score"
    else:
        if "X" in df.columns:
            df["display_score"] = pd.to_numeric(df["X"], errors="coerce").fillna(0.0)
            raw_score_label = "X"
        elif "T" in df.columns:
            df["display_score"] = pd.to_numeric(df["T"], errors="coerce").fillna(0.0)
            raw_score_label = "T"
        else:
            score_col = coerce_first_existing_col(df, ["display_score", "score"])
            if score_col:
                df["display_score"] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
                raw_score_label = score_col
            else:
                df["display_score"] = 0.0
                raw_score_label = "display_score"

    if label_col and label_col != "display_label":
        df["display_label"] = df[label_col].astype(str)
    elif not label_col:
        df["display_label"] = np.arange(len(df)).astype(str)

    tf_id_col = coerce_first_existing_col(df, ["TF", "tf", "source", "id", "uid"])
    if "display_name" in df.columns and tf_id_col:
        same_as_id = df["display_name"].astype(str).str.strip() == df[tf_id_col].astype(str).str.strip()
        better_label = (~same_as_id) & df["display_name"].astype(str).str.strip().ne("")
        df.loc[better_label, "display_label"] = (
            df.loc[better_label, "display_name"].astype(str) + " (" + df.loc[better_label, tf_id_col].astype(str) + ")"
        )

    df = df.dropna(subset=["display_score"]).copy()
    # Only filter by p-value if enabled
    if use_pval_filter:
        df = df[df["display_score"] >= float(min_display_score)].copy()
    df = df.sort_values("display_score", ascending=False).copy()

    if top_n > 0:
        df = df.head(int(top_n)).copy()

    preferred = [c for c in ["display_name", "display_label", "display_score", "TF", "symbol", "X", "T"] if c in df.columns]
    remaining = [c for c in df.columns if c not in preferred]
    df = df[preferred + remaining]

    return df.reset_index(drop=True), raw_score_label


def prepare_edge_display_df(
    edge_df: pd.DataFrame,
    selected_tf_df: pd.DataFrame,
    min_abs_edge_score: float,
    max_edges: int,
) -> pd.DataFrame:
    if edge_df is None or edge_df.empty:
        return pd.DataFrame()

    df = edge_df.copy()
    df.columns = [strip_bom_text(c) for c in df.columns]

    src_col = coerce_first_existing_col(df, ["source", "srcuid", "src"])
    trg_col = coerce_first_existing_col(df, ["target", "trguid", "trg"])
    edge_score_col = coerce_first_existing_col(df, ["edge_score", "score"])
    tf_score_col = coerce_first_existing_col(df, ["X_tf", "tf_score"])

    print(f"[DEBUG-PREP] src_col={src_col!r} trg_col={trg_col!r} edge_df shape={df.shape} cols={list(df.columns[:8])}")

    if not src_col or not trg_col:
        return pd.DataFrame()

    df[src_col] = df[src_col].astype(str).str.strip()
    df[trg_col] = df[trg_col].astype(str).str.strip()
    df = df[(df[src_col].str.len() > 0) & (df[trg_col].str.len() > 0)].copy()

    if edge_score_col:
        df["_edge_score"] = pd.to_numeric(df[edge_score_col], errors="coerce").fillna(0.0)
    else:
        df["_edge_score"] = 0.0

    if tf_score_col:
        df["_tf_score"] = pd.to_numeric(df[tf_score_col], errors="coerce").fillna(0.0)
    else:
        df["_tf_score"] = 0.0

    df["_abs_edge_score"] = df["_edge_score"].abs()

    # Don't filter by selected TF here — let build_regulatory_network_figure
    # select regulators by significance from tf_score_map. This ensures top-significant
    # TFs always appear in the graph, even if they have few edges.
    # if selected_tf_df is not None and not selected_tf_df.empty:
    #     tf_id_col = coerce_first_existing_col(selected_tf_df, ["TF", "source", "uid", "id", "tf_id"])
    #     if tf_id_col:
    #         allowed = set(selected_tf_df[tf_id_col].astype(str))
    #         filtered = df[df[src_col].astype(str).isin(allowed)].copy()
    #         if not filtered.empty:
    #             df = filtered

    # FIRST: Reserve edges for top regulators (before filtering by score)
    df_top_regs = pd.DataFrame()
    if selected_tf_df is not None and not selected_tf_df.empty and max_edges > 0:
        tf_id_col = coerce_first_existing_col(selected_tf_df, ["TF", "source", "uid", "id", "tf_id"])
        print(f"[DEBUG-PREP] tf_id_col={tf_id_col!r} selected_tf_df has {len(selected_tf_df)} rows")
        if tf_id_col:
            # Get the top regulator IDs IN ORDER OF SIGNIFICANCE (from selected_tf_df)
            top_tf_ids_ordered = selected_tf_df[tf_id_col].astype(str).str.strip().tolist()
            print(f"[DEBUG-PREP] top regulators in order: {top_tf_ids_ordered[:5]}")
            # For each top regulator, grab their best edge (before filtering by score)
            edges_per_tf = []
            found_regs = []
            for tf_id in top_tf_ids_ordered:
                tf_edges = df[df[src_col].astype(str) == tf_id].sort_values("_abs_edge_score", ascending=False)
                if not tf_edges.empty:
                    edges_per_tf.append(tf_edges.iloc[0:1])
                    found_regs.append(tf_id)
            print(f"[DEBUG-PREP] searched for {top_tf_ids_ordered[:5]} in src_col={src_col!r}, found {found_regs}")

            if edges_per_tf:
                df_top_regs = pd.concat(edges_per_tf, ignore_index=True)
                print(f"[DEBUG-PREP] reserved {len(df_top_regs)} edges for {len(edges_per_tf)} top regulators")

    # THEN: Filter by score and select remaining edges
    df = df[df["_abs_edge_score"] >= float(min_abs_edge_score)].copy()
    df = df.sort_values(["_abs_edge_score", "_tf_score"], ascending=[False, False]).copy()

    if max_edges > 0:
        # Calculate how many remaining edges we can take
        remaining_budget = max(0, int(max_edges) - len(df_top_regs))
        df = df.head(remaining_budget).copy()
        # Combine: top regulator edges + remaining high-score edges
        if not df_top_regs.empty:
            df = pd.concat([df_top_regs, df], ignore_index=True)

    result_regs = df[src_col].unique() if src_col in df.columns else []
    print(f"[DEBUG-PREP-FINAL] returning {len(df)} edges with regulators: {sorted(result_regs)[:10]}")
    return df.reset_index(drop=True)


def regulator_bar_plot(
    tf_df: pd.DataFrame,
    score_label: str,
    title: str,
    top_n: int = 20,
) -> go.Figure:
    fig = go.Figure()
    if tf_df is None or tf_df.empty:
        fig.update_layout(title=title, height=420)
        return fig

    plot_df = tf_df.head(top_n).copy()
    fig.add_trace(
        go.Bar(
            x=plot_df["display_label"].astype(str),
            y=plot_df["display_score"].astype(float),
            text=np.round(plot_df["display_score"].astype(float), 4),
            textposition="outside",
        )
    )
    fig.update_layout(
        title=title,
        height=430,
        margin=dict(l=40, r=20, t=60, b=150),
        xaxis_title="Regulator",
        yaxis_title=score_label,
    )
    fig.update_xaxes(tickangle=45)
    return fig


def pathway_bar_plot(pathway_df: pd.DataFrame, top_n: int = 20) -> go.Figure:
    fig = go.Figure()
    if pathway_df is None or pathway_df.empty:
        fig.update_layout(title="Top pathways", height=420)
        return fig

    df = pathway_df.copy()
    label_col = coerce_first_existing_col(df, ["pathway", "name", "term", "id"])
    score_col = coerce_first_existing_col(df, ["pvalue", "pval", "fdr"])

    if not label_col or not score_col:
        fig.update_layout(title="Top pathways", height=420)
        return fig

    df[score_col] = pd.to_numeric(df[score_col], errors="coerce").clip(lower=1e-300)
    df["plot_score"] = -np.log10(df[score_col])
    df = df.sort_values("plot_score", ascending=False).head(top_n)

    fig.add_trace(
        go.Bar(
            x=df[label_col].astype(str),
            y=df["plot_score"].astype(float),
            text=np.round(df["plot_score"].astype(float), 4),
            textposition="outside",
        )
    )
    fig.update_layout(
        title="Top pathways",
        height=430,
        margin=dict(l=40, r=20, t=60, b=150),
        xaxis_title="Pathway",
        yaxis_title="-log10(p-value)",
    )
    fig.update_xaxes(tickangle=45)
    return fig


def build_regulatory_network_figure(
    edge_df: pd.DataFrame,
    tf_name_map: Optional[Dict[str, str]] = None,
    target_name_map: Optional[Dict[str, str]] = None,
    title: str = "Regulatory Network",
    max_nodes: int = 40,
    label_top_hubs: int = 12,
    tf_score_map: Optional[Dict[str, float]] = None,
) -> go.Figure:
    fig = go.Figure()
    tf_name_map = tf_name_map or {}
    tf_score_map = tf_score_map or {}
    target_name_map = target_name_map or {}

    if edge_df is None or edge_df.empty:
        fig.update_layout(
            title=title,
            height=980,
            annotations=[dict(text="No edge data available after filtering", x=0.5, y=0.5, showarrow=False)],
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    df = edge_df.copy()
    df.columns = [strip_bom_text(c) for c in df.columns]

    src_col = coerce_first_existing_col(df, ["source", "srcuid", "src"])
    trg_col = coerce_first_existing_col(df, ["target", "trguid", "trg"])
    edge_score_col = coerce_first_existing_col(df, ["edge_score", "score"])

    if not src_col or not trg_col:
        fig.update_layout(
            title=title,
            height=980,
            annotations=[dict(text="Missing source / target columns in edge table", x=0.5, y=0.5, showarrow=False)],
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    df["_src_uid"] = df[src_col].astype(str).str.strip()
    df["_trg_uid"] = df[trg_col].astype(str).str.strip()
    df = df[(df["_src_uid"].str.len() > 0) & (df["_trg_uid"].str.len() > 0)].copy()

    if edge_score_col:
        df["_edge_score"] = pd.to_numeric(df[edge_score_col], errors="coerce").fillna(0.0)
    else:
        df["_edge_score"] = 0.0

    reg_counts = df["_src_uid"].value_counts().to_dict()
    tgt_counts = df["_trg_uid"].value_counts().to_dict()

    max_regs = max(1, max_nodes // 2)
    max_tgts = max(1, max_nodes // 2)

    def _best_label(uid: str, mapping: Dict[str, str]) -> str:
        uid = str(uid).strip()
        mapped = str(mapping.get(uid, uid)).strip()
        return mapped if mapped else uid

    def _has_readable_name(uid: str, mapping: Dict[str, str]) -> bool:
        uid = str(uid).strip()
        mapped = str(mapping.get(uid, uid)).strip()
        return bool(mapped) and mapped != uid

    # FORCE top-significance regulators to appear by selecting from tf_score_map FIRST
    if tf_score_map:
        # Start with ALL significant regulators from tf_score_map, ranked by score
        sigs_ranked = sorted(tf_score_map.keys(), key=lambda x: -tf_score_map[x])
        regs = sigs_ranked[:max_regs]
        print(f"[DEBUG-NET-V2] using tf_score_map, regs_forced={regs[:5]}")
    else:
        # Fallback: get all unique regulators from edges
        all_regs = list(set(df["_src_uid"].astype(str).tolist()))
        all_regs_sorted = sorted(all_regs, key=lambda x: (-reg_counts.get(x, 0), _best_label(x, tf_name_map)))
        regs = all_regs_sorted[:max_regs]
        print(f"[DEBUG-NET-V2] no tf_score_map, using edge counts, regs={regs[:5]}")

    tgts = (
        df["_trg_uid"].value_counts()
        .sort_values(ascending=False)
        .head(max_tgts)
        .index.astype(str)
        .tolist()
    )

    # Filter edges, but if a selected regulator has no edges to selected targets,
    # include edges from that regulator to ANY target
    df_filtered = df[df["_src_uid"].isin(regs) & df["_trg_uid"].isin(tgts)].copy()
    if df_filtered.empty and regs:
        # Fallback: if no edges between top regs and top targets, show all edges from top regs
        df_filtered = df[df["_src_uid"].isin(regs)].copy()
        if not df_filtered.empty:
            tgts = list(set(df_filtered["_trg_uid"].astype(str).tolist()))

    df = df_filtered if not df_filtered.empty else df.copy()

    # Don't filter regs here - we want to keep force-selected top-significant regulators even if they have no edges
    # regs = [r for r in regs if r in df["_src_uid"].astype(str).values]

    # Only keep edges from selected regulators
    df = df[df["_src_uid"].isin(regs)].copy()

    tgts = list(dict.fromkeys(df["_trg_uid"].astype(str).tolist()))

    print(f"[DEBUG-NET] max_regs={max_regs} regs_selected={regs[:5]} score_map_keys={list(tf_score_map.keys())[:5]}")
    tgts = sorted(tgts, key=lambda x: (-tgt_counts.get(x, 0), _best_label(x, target_name_map)))

    pos = {}
    left_x = -2.8
    right_x = 2.8

    reg_y = np.linspace(2.8, -2.8, len(regs)) if regs else []
    tgt_y = np.linspace(2.8, -2.8, len(tgts)) if tgts else []

    for i, r in enumerate(regs):
        pos[r] = (left_x, float(reg_y[i]))
    for i, t in enumerate(tgts):
        pos[t] = (right_x, float(tgt_y[i]))

    pos_edge_x, pos_edge_y, pos_hover = [], [], []
    neg_edge_x, neg_edge_y, neg_hover = [], [], []
    zero_edge_x, zero_edge_y, zero_hover = [], [], []

    for _, row in df.iterrows():
        s = str(row["_src_uid"])
        t = str(row["_trg_uid"])
        score = float(row["_edge_score"])

        if s not in pos or t not in pos:
            continue

        x0, y0 = pos[s]
        x1, y1 = pos[t]
        hover = f"{_best_label(s, tf_name_map)} -> {_best_label(t, target_name_map)}<br>score={score:g}"

        if score > 0:
            pos_edge_x += [x0, x1, None]
            pos_edge_y += [y0, y1, None]
            pos_hover += [hover, hover, None]
        elif score < 0:
            neg_edge_x += [x0, x1, None]
            neg_edge_y += [y0, y1, None]
            neg_hover += [hover, hover, None]
        else:
            zero_edge_x += [x0, x1, None]
            zero_edge_y += [y0, y1, None]
            zero_hover += [hover, hover, None]

    if pos_edge_x:
        fig.add_trace(
            go.Scatter(
                x=pos_edge_x,
                y=pos_edge_y,
                mode="lines",
                line=dict(width=2.4, color="rgba(86, 193, 108, 0.42)"),
                hoverinfo="text",
                text=pos_hover,
                name="positive / agreement",
            )
        )

    if neg_edge_x:
        fig.add_trace(
            go.Scatter(
                x=neg_edge_x,
                y=neg_edge_y,
                mode="lines",
                line=dict(width=2.4, color="rgba(234, 108, 108, 0.40)"),
                hoverinfo="text",
                text=neg_hover,
                name="negative / disagreement",
            )
        )

    if zero_edge_x:
        fig.add_trace(
            go.Scatter(
                x=zero_edge_x,
                y=zero_edge_y,
                mode="lines",
                line=dict(width=1.2, color="rgba(160, 174, 192, 0.28)"),
                hoverinfo="text",
                text=zero_hover,
                name="neutral",
            )
        )

    reg_sizes = [30 + 3.0 * min(reg_counts.get(n, 1), 10) for n in regs]
    tgt_sizes = [22 + 2.4 * min(tgt_counts.get(n, 1), 10) for n in tgts]

    top_regs_for_labels = set(regs[: max(0, int(label_top_hubs))])
    top_tgts_for_labels = set(tgts[: max(0, int(label_top_hubs))])

    reg_text = [
        _best_label(n, tf_name_map) if (_has_readable_name(n, tf_name_map) or n in top_regs_for_labels) else ""
        for n in regs
    ]
    tgt_text = [
        _best_label(n, target_name_map) if (_has_readable_name(n, target_name_map) or n in top_tgts_for_labels) else ""
        for n in tgts
    ]

    fig.add_trace(
        go.Scatter(
            x=[pos[n][0] for n in regs],
            y=[pos[n][1] for n in regs],
            mode="markers+text",
            marker=dict(
                size=reg_sizes,
                color="#ff4d6d",
                line=dict(width=2.4, color="white"),
            ),
            text=reg_text,
            textposition="middle left",
            textfont=dict(size=20, color=THEME["title"]),
            hovertemplate=[
                f"regulator={_best_label(n, tf_name_map)}<br>uid={n}<br>degree={reg_counts.get(n, 0)}<extra></extra>"
                for n in regs
            ],
            name="regulators",
        )
    )

    fig.add_trace(
        go.Scatter(
            x=[pos[n][0] for n in tgts],
            y=[pos[n][1] for n in tgts],
            mode="markers+text",
            marker=dict(
                size=tgt_sizes,
                color="#2ec4b6",
                line=dict(width=2.2, color="white"),
            ),
            text=tgt_text,
            textposition="middle right",
            textfont=dict(size=18, color=THEME["title"]),
            hovertemplate=[
                f"target={_best_label(n, target_name_map)}<br>uid={n}<br>degree={tgt_counts.get(n, 0)}<extra></extra>"
                for n in tgts
            ],
            name="targets",
        )
    )

    fig.update_layout(
        title=title,
        height=980,
        margin=dict(l=110, r=110, t=80, b=50),
        xaxis=dict(visible=False, range=[-4.1, 4.1]),
        yaxis=dict(visible=False, range=[-3.2, 3.2]),
        showlegend=True,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(
            x=0.80,
            y=0.98,
            bgcolor="rgba(255,255,255,0.90)",
            bordercolor="#DCEBFF",
            borderwidth=1,
        ),
    )
    return fig


def make_summary_cards(summary: Dict[str, str]) -> html.Div:
    items = []
    for label, value, color in summary.get("cards", []):
        items.append(stat_card(label, value, color))
    return html.Div(
        style={"display": "grid", "gridTemplateColumns": "repeat(4, minmax(0,1fr))", "gap": "10px", "marginBottom": "14px"},
        children=items,
    )


def build_top_regulator_preview(tf_df: pd.DataFrame, n: int = 10) -> dash_table.DataTable:
    if tf_df is None or tf_df.empty:
        return df_table(pd.DataFrame({"message": ["No regulator table available"]}), page_size=5)

    preview_cols = [
        c for c in [
            "display_name",
            "display_label",
            "TF",
            "symbol",
            "display_score",
            "X",
            "T",
            "source",
            "name",
            "posterior_sd",
            "targets_found",
            "pvalue",
            "pval",
            "fdr",
            "proteinsfound",
        ] if c in tf_df.columns
    ]
    preview_df = tf_df[preview_cols].head(n).copy() if preview_cols else tf_df.head(n).copy()
    return df_table(preview_df, page_size=n)


# =============================================================================
# ORNOR quality helpers
# =============================================================================

ORNOR_GR_RE = re.compile(r"max GR=([0-9]*\.?[0-9]+)", re.IGNORECASE)
ORNOR_LAST_GR_RE = re.compile(r"last max GR=([0-9]*\.?[0-9]+)", re.IGNORECASE)

ORNOR_EVIDENCE_RE = re.compile(r"evidence:\s*(\d+)\s+genes", re.IGNORECASE)
ORNOR_EVIDENCE_LOOSE_RE = re.compile(r"overlap after loose evidence fallback:\s*(\d+)", re.IGNORECASE)

ORNOR_TF_CAND_RE = re.compile(r"active TFs[^:]*:\s*(\d+)", re.IGNORECASE)


def notice_box(title: str, text: str, accent: str, bg: str) -> html.Div:
    return html.Div(
        style={
            "background": bg,
            "borderLeft": f"5px solid {accent}",
            "borderRadius": "12px",
            "padding": "12px 14px",
            "marginBottom": "10px",
        },
        children=[
            html.Div(title, style={"fontWeight": "900", "color": THEME["title"], "marginBottom": "4px"}),
            html.Div(text, style={"color": THEME["text"], "fontSize": "13px", "lineHeight": "1.45"}),
        ],
    )


def _extract_last_float(pattern: re.Pattern, text: str) -> Optional[float]:
    matches = pattern.findall(text or "")
    if not matches:
        return None
    try:
        return float(matches[-1])
    except Exception:
        return None


def _extract_last_int(pattern: re.Pattern, text: str) -> Optional[int]:
    matches = pattern.findall(text or "")
    if not matches:
        return None
    try:
        return int(matches[-1])
    except Exception:
        return None


def _series_is_flat(values: pd.Series, tol: float = 1e-12) -> bool:
    vals = pd.to_numeric(values, errors="coerce").dropna()
    if len(vals) < 2:
        return False
    return bool(float(vals.max()) - float(vals.min()) <= tol)


def detect_ornor_quality(
    tf_df: pd.DataFrame,
    edge_df: pd.DataFrame,
    stdout_txt: str,
) -> Dict[str, object]:
    text = stdout_txt or ""
    text_lower = text.lower()

    last_gr = _extract_last_float(ORNOR_GR_RE, text)
    if last_gr is None:
        last_gr = _extract_last_float(ORNOR_LAST_GR_RE, text)

    evidence_overlap = _extract_last_int(ORNOR_EVIDENCE_RE, text)
    if evidence_overlap is None:
        evidence_overlap = _extract_last_int(ORNOR_EVIDENCE_LOOSE_RE, text)

    tf_candidates = _extract_last_int(ORNOR_TF_CAND_RE, text)
    if tf_candidates is None and tf_df is not None and not tf_df.empty:
        tf_candidates = int(len(tf_df))

    flat_x = False
    flat_t = False
    if tf_df is not None and not tf_df.empty:
        if "X" in tf_df.columns:
            flat_x = _series_is_flat(tf_df["X"])
        if "T" in tf_df.columns:
            flat_t = _series_is_flat(tf_df["T"])

    flat_detected = bool(flat_x and flat_t)
    flat_logged = (
        "flat/degenerate" in text_lower
        or "posterior is flat" in text_lower
        or "symmetric posterior" in text_lower
    )

    used_fallback = (
        "using fallback overlap ranking" in text_lower
        or "fallback overlap ranking" in text_lower
    )
    posterior_failed = "posterior sampling unavailable/failed" in text_lower

    if used_fallback:
        posterior_type = "Fallback ranking"
        posterior_color = THEME["warn"]
        posterior_desc = "Bayesian posterior was unavailable or unusable, so the app is showing overlap-based ranking instead."
    elif flat_detected or flat_logged:
        posterior_type = "Weak posterior"
        posterior_color = THEME["bad"]
        posterior_desc = "ORNOR produced a flat or nearly symmetric posterior. TF ranking is not strongly informative."
    else:
        posterior_type = "Bayesian posterior"
        posterior_color = THEME["good"]
        posterior_desc = "ORNOR produced a non-flat Bayesian posterior suitable for interpretation."

    if used_fallback:
        convergence_label = "N/A"
        convergence_color = THEME["warn"]
        convergence_desc = "Fallback mode does not represent a valid Bayesian convergence result."
    elif last_gr is None:
        convergence_label = "Unknown"
        convergence_color = THEME["warn"]
        convergence_desc = "No Gelman-Rubin value was recovered from the ORNOR run log."
    elif last_gr <= 1.10:
        convergence_label = "Good convergence"
        convergence_color = THEME["good"]
        convergence_desc = f"Final max GR = {last_gr:.4f}. This is consistent with a well-converged run."
    elif last_gr <= 1.25:
        convergence_label = "Usable but weak convergence"
        convergence_color = THEME["warn"]
        convergence_desc = f"Final max GR = {last_gr:.4f}. Posterior may be usable, but should be interpreted cautiously."
    else:
        convergence_label = "Poor convergence"
        convergence_color = THEME["bad"]
        convergence_desc = f"Final max GR = {last_gr:.4f}. Posterior likely did not mix well."

    if used_fallback:
        confidence_label = "Exploratory only"
        confidence_color = THEME["warn"]
    elif flat_detected or flat_logged:
        confidence_label = "Exploratory only"
        confidence_color = THEME["bad"]
    elif last_gr is not None and last_gr <= 1.10 and (evidence_overlap or 0) >= 10 and (tf_candidates or 0) >= 5:
        confidence_label = "Higher confidence"
        confidence_color = THEME["good"]
    elif last_gr is not None and last_gr <= 1.25 and (evidence_overlap or 0) >= 5:
        confidence_label = "Interpret cautiously"
        confidence_color = THEME["warn"]
    else:
        confidence_label = "Exploratory only"
        confidence_color = THEME["warn"]

    warnings = []

    if used_fallback:
        warnings.append({
            "title": "Fallback ranking used",
            "text": "These ORNOR results are not a strong Bayesian posterior. They are overlap-based fallback scores and should be treated as exploratory.",
            "accent": THEME["warn"],
            "bg": "#FFF7E8",
        })

    if (flat_detected or flat_logged) and not used_fallback:
        warnings.append({
            "title": "Flat / symmetric posterior detected",
            "text": "Most TFs received nearly identical posterior values. This usually means the model could not separate regulators meaningfully for this run.",
            "accent": THEME["bad"],
            "bg": "#FFF1F2",
        })

    if last_gr is not None and not used_fallback:
        if last_gr > 1.25:
            warnings.append({
                "title": "Poor convergence",
                "text": f"Gelman-Rubin remained high (max GR = {last_gr:.4f}). The posterior may be unstable.",
                "accent": THEME["bad"],
                "bg": "#FFF1F2",
            })
        elif last_gr > 1.10:
            warnings.append({
                "title": "Weak convergence",
                "text": f"Gelman-Rubin is acceptable but not ideal (max GR = {last_gr:.4f}). Interpret regulator ranking with caution.",
                "accent": THEME["warn"],
                "bg": "#FFF7E8",
            })

    if evidence_overlap is not None and evidence_overlap < 5:
        warnings.append({
            "title": "Low evidence overlap",
            "text": f"Only {evidence_overlap} evidence genes overlapped network targets. This limits ORNOR's ability to distinguish regulators.",
            "accent": THEME["warn"],
            "bg": "#FFF7E8",
        })

    return {
        "posterior_type": posterior_type,
        "posterior_color": posterior_color,
        "posterior_desc": posterior_desc,
        "convergence_label": convergence_label,
        "convergence_color": convergence_color,
        "convergence_desc": convergence_desc,
        "confidence_label": confidence_label,
        "confidence_color": confidence_color,
        "last_gr": last_gr,
        "evidence_overlap": evidence_overlap,
        "tf_candidates": tf_candidates,
        "flat_detected": bool(flat_detected or flat_logged),
        "used_fallback": used_fallback,
        "posterior_failed": posterior_failed,
        "warnings": warnings,
    }


def render_ornor_quality_strip(q: Dict[str, object]) -> html.Div:
    return html.Div(
        style={
            "display": "flex",
            "flexWrap": "wrap",
            "gap": "8px",
            "marginBottom": "10px",
        },
        children=[
            badge(str(q["posterior_type"]), str(q["posterior_color"])),
            badge(str(q["convergence_label"]), str(q["convergence_color"])),
            badge(str(q["confidence_label"]), str(q["confidence_color"])),
        ],
    )


def render_ornor_quality_panel(q: Dict[str, object]) -> html.Div:
    cards = [
        stat_card("Posterior type", str(q["posterior_type"]), str(q["posterior_color"])),
        stat_card("Convergence", str(q["convergence_label"]), str(q["convergence_color"])),
        stat_card("Evidence overlap", str(q["evidence_overlap"] if q["evidence_overlap"] is not None else "Unknown"), THEME["blue"]),
        stat_card("TF candidates", str(q["tf_candidates"] if q["tf_candidates"] is not None else "Unknown"), THEME["purple"]),
        stat_card("Interpretation", str(q["confidence_label"]), str(q["confidence_color"])),
    ]

    warning_boxes = [
        notice_box(w["title"], w["text"], w["accent"], w["bg"])
        for w in q.get("warnings", [])
    ]

    if not warning_boxes:
        warning_boxes = [
            notice_box(
                "No major ORNOR quality warning",
                "This run produced a non-flat posterior without an obvious fallback warning.",
                THEME["good"],
                "#ECFDF5",
            )
        ]

    return card("ORNOR quality summary", [
        render_ornor_quality_strip(q),
        html.Div(
            style={
                "display": "grid",
                "gridTemplateColumns": "repeat(5, minmax(0,1fr))",
                "gap": "10px",
                "marginBottom": "12px",
            },
            children=cards,
        ),
        info_row("Posterior assessment", str(q["posterior_desc"])),
        info_row("Convergence assessment", str(q["convergence_desc"])),
        html.Div(style={"height": "10px"}),
        *warning_boxes,
    ])


# =============================================================================
# Dash app
# =============================================================================

app = Dash(__name__, suppress_callback_exceptions=True)
app.scripts.config.serve_locally = True
app.css.config.serve_locally = True
server = app.server
server.config["MAX_CONTENT_LENGTH"] = 1024 * 1024 * 1024

app.layout = html.Div(
    style={"display": "flex", "background": THEME["page_bg"], "minHeight": "100vh"},
    children=[
        html.Div(
            style={
                "width": "430px",
                "padding": "16px",
                "background": "linear-gradient(180deg, #FFFFFF 0%, #EAF3FF 100%)",
                "borderRight": f"1px solid {THEME['panel_border']}",
                "overflowY": "auto",
            },
            children=[
                html.Div("Inference Dash", style={"fontSize": "24px", "fontWeight": "900", "color": THEME["title"]}),
                html.Div("Research workflow for regulatory inference with CIE and ORNOR.", style={"color": THEME["muted"], "marginBottom": "14px"}),

                dcc.Store(id="expr_path_store"),
                dcc.Store(id="expr_meta_store"),
                dcc.Store(id="job_id_store"),
                dcc.Store(id="validation_ok_store"),
                dcc.Store(id="rels_path_store"),
                dcc.Store(id="ents_path_store"),
                dcc.Store(id="results_tab_store", data="summary"),

                card("1) Engine", [
                    dcc.RadioItems(
                        id="engine_choice",
                        options=[{"label": " CIE", "value": "CIE"}, {"label": " ORNOR", "value": "ORNOR"}],
                        value="CIE",
                        labelStyle={"display": "block", "margin": "8px 0", "fontWeight": "700", "color": THEME["text"]},
                        inputStyle={"marginRight": "10px"},
                    )
                ]),

                card("2) Expression file (required)", [
                    dcc.Upload(
                        id="expr_upload",
                        children=html.Div(["Drag & drop or ", html.B("click to upload")]),
                        style=UPLOAD_STYLE,
                        multiple=False,
                    ),
                    html.Div(id="expr_upload_status", style={"marginTop": "10px", "fontSize": "12px", "color": THEME["muted"]}),
                ]),

                card("3) Signature column mapping", [
                    html.Div("Choose columns from the uploaded file. Auto-detection is applied first.", style={"fontSize": "12px", "color": THEME["muted"], "marginBottom": "10px"}),

                    html.Div("Gene column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="gene_col_dd", options=[{"label": "-- None --", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("Entrez column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="entrez_col_dd", options=[{"label": "-- None --", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("Fold-change column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="fc_col_dd", options=[{"label": "-- None --", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("P-value / adjusted p-value column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="pval_col_dd", options=[{"label": "-- None --", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("Direction column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="direction_col_dd", options=[{"label": "-- None --", "value": "__NONE__"}], value="__NONE__", clearable=False),
                ]),

                card("4) Filtering / normalization", [
                    dcc.Checklist(
                        id="use_pval_filter_ck",
                        options=[{"label": " Use p-value filtering", "value": "USE_PVAL"}],
                        value=["USE_PVAL"],
                    ),
                    html.Div("P-value threshold", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="pval_thresh_in", type="number", value=0.05, step=0.001, min=0, max=1, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    dcc.Checklist(
                        id="use_abs_fc_filter_ck",
                        options=[{"label": " Use absolute fold-change filtering", "value": "USE_ABS_FC"}],
                        value=[],
                    ),
                    html.Div("Absolute fold-change threshold", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="abs_fc_thresh_in", type="number", value=1.5, step=0.1, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    dcc.Checklist(
                        id="skip_mcmc_ck",
                        options=[{"label": " Skip MCMC (fast enrichment only, ORNOR)", "value": "SKIP_MCMC"}],
                        value=[],
                    ),

                    html.Div(id="validation_inline_hint", style={"marginTop": "10px", "fontSize": "12px", "color": THEME["muted"]}),
                ]),

                card("5) Network / database", [
                    dcc.Dropdown(
                        id="network_choice",
                        options=[{"label": k, "value": k} for k in DEFAULT_NETWORKS.keys()] + [{"label": "Upload custom network...", "value": "__UPLOAD__"}],
                        value=preferred_network_choice_for_engine("CIE") if DEFAULT_NETWORKS else "__UPLOAD__",
                        clearable=False,
                    ),
                    html.Div(
                        id="custom_network_box",
                        style={"marginTop": "12px", "display": "none"},
                        children=[
                            html.Div("Network .rels (required)", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                            dcc.Upload(id="rels_upload", children=html.Div(["Upload .rels"]), style=SMALL_UPLOAD_STYLE, multiple=False),
                            html.Div(id="rels_upload_status", style={"marginTop": "8px", "fontSize": "12px", "color": THEME["muted"]}),

                            html.Div(style={"height": "10px"}),
                            html.Div("Entities (required for CIE, optional for ORNOR)", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                            dcc.Upload(id="ents_upload", children=html.Div(["Upload .entities / .ents"]), style=SMALL_UPLOAD_STYLE, multiple=False),
                            html.Div(id="ents_upload_status", style={"marginTop": "8px", "fontSize": "12px", "color": THEME["muted"]}),
                        ],
                    ),
                ]),

                card("6) Run", [
                    html.Button("Run inference", id="run_btn", n_clicks=0, style=BTN_DISABLED),
                    html.Button("Clear results", id="clear_results_btn", n_clicks=0, style={
                        "width": "100%", "marginTop": "8px", "padding": "6px",
                        "background": "transparent", "color": THEME["muted"],
                        "border": f"1px solid {THEME['panel_border']}", "borderRadius": "6px",
                        "cursor": "pointer", "fontSize": "12px",
                    }),
                    html.Div(id="run_hint", style={"marginTop": "10px", "fontSize": "12px", "color": THEME["muted"]}),

                    html.Div(style={"height": "10px"}),
                    html.Div("Progress", style={"fontSize": "12px", "fontWeight": "800", "color": THEME["muted"]}),
                    html.Div(
                        id="progress_bar",
                        style={
                            "width": "100%",
                            "height": "16px",
                            "borderRadius": "999px",
                            "background": "#D9ECFF",
                            "overflow": "hidden",
                            "border": "1px solid #CFE3FF",
                        },
                        children=[html.Div(id="progress_fill", style=progress_fill_style(0))],
                    ),
                    html.Div(id="progress_text", style={"marginTop": "8px", "fontSize": "12px", "color": THEME["muted"]}),

                    html.Div(style={"height": "10px"}),
                    html.Button("Download outputs", id="download_btn", n_clicks=0, style=BTN_SECONDARY),
                    dcc.Download(id="download_outputs"),
                    html.Div(id="download_hint", style={"marginTop": "8px", "fontSize": "12px", "color": THEME["muted"]}),
                ]),

                card("7) Result display controls", [
                    html.Div("Top N regulators", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="top_n_regulators_in", type="number", value=20, min=1, step=1, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    html.Div("Filter by p-value", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Checklist(id="use_pval_filter_check", options=[{"label": " Enable p-value filtering", "value": "yes"}], value=[], style={"fontSize": "12px"}),

                    html.Div(style={"height": "10px"}),
                    html.Div("Minimum regulator score (if p-value filtering enabled)", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="min_reg_score_in", type="number", value=0.0, min=0, step=0.01, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    html.Div("Minimum absolute edge score", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="min_abs_edge_score_in", type="number", value=1.0, min=0, step=0.01, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    html.Div("Maximum visible nodes", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="max_nodes_in", type="number", value=40, min=4, step=1, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    html.Div("Maximum visible edges", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="max_edges_in", type="number", value=100, min=1, step=1, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    html.Div("Label top hubs", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Input(id="label_top_hubs_in", type="number", value=12, min=0, step=1, style=INPUT_STYLE),

                    html.Div(style={"height": "10px"}),
                    html.Button("Apply display filters", id="apply_display_btn", n_clicks=0, style={
                        "width": "100%", "padding": "8px", "background": THEME["blue"],
                        "color": "#fff", "border": "none", "borderRadius": "6px",
                        "fontWeight": "700", "cursor": "pointer", "fontSize": "13px",
                    }),
                    html.Div(
                        "These controls affect the regulator table, regulator plot, edge table, and network graph.",
                        style={"marginTop": "10px", "fontSize": "12px", "color": THEME["muted"]},
                    ),
                ]),
            ],
        ),

        html.Div(
            style={"flex": "1", "padding": "18px 22px", "overflowY": "auto"},
            children=[
                dcc.Interval(id="poll_interval", interval=1200, n_intervals=0, disabled=False),

                html.Div(
                    style={"display": "flex", "alignItems": "center", "justifyContent": "space-between", "gap": "10px", "marginBottom": "10px"},
                    children=[
                        html.Div(
                            [
                                html.Div("Inference Results", style={"fontSize": "22px", "fontWeight": "900", "color": THEME["title"]}),
                                html.Div(id="results_message", style={"color": THEME["muted"]}),
                            ]
                        ),
                        html.Div(id="results_badge"),
                    ],
                ),

                html.Div(id="signature_preview_block"),
                html.Div(id="validation_summary_block"),
                html.Div(id="normalized_preview_block"),
                html.Div(id="results_blocks"),
            ],
        ),
    ],
)


# =============================================================================
# Upload callbacks
# =============================================================================

@app.callback(Output("custom_network_box", "style"), Input("network_choice", "value"))
def toggle_custom(choice):
    return {"marginTop": "12px", "display": "block"} if choice == "__UPLOAD__" else {"marginTop": "12px", "display": "none"}


@app.callback(
    Output("network_choice", "value"),
    Input("engine_choice", "value"),
    State("network_choice", "value"),
    prevent_initial_call=True,
)
def sync_default_network_to_engine(engine, current_choice):
    if current_choice == "__UPLOAD__":
        return no_update

    current_choice = current_choice or ""
    preferred = preferred_network_choice_for_engine(engine)

    if engine == "ORNOR":
        if current_choice == "tcChIP (all_tissues.rels) [default]":
            return preferred
    elif engine == "CIE":
        if current_choice == "tcChIP (three_tissues.rels)":
            return preferred

    return no_update


@app.callback(
    Output("expr_upload_status", "children"),
    Output("expr_path_store", "data"),
    Output("expr_meta_store", "data"),
    Input("expr_upload", "contents"),
    State("expr_upload", "filename"),
    prevent_initial_call=True,
)
def on_expr_upload(contents, filename):
    if not contents:
        return badge("No file", THEME["warn"]), None, None
    try:
        p = save_upload(contents, filename, UPLOADS_DIR)
        meta = build_signature_metadata(p)
        msg = html.Span([
            badge("Uploaded", THEME["good"]),
            f"  {p.name} ({p.stat().st_size:,} bytes) -- {meta['rows']:,} rows, {meta['cols']} columns",
        ])
        return msg, str(p), meta
    except Exception as e:
        return badge(f"Upload failed: {e}", THEME["bad"]), None, None


@app.callback(
    Output("rels_upload_status", "children"),
    Output("rels_path_store", "data"),
    Input("rels_upload", "contents"),
    State("rels_upload", "filename"),
    prevent_initial_call=True,
)
def on_rels_upload(contents, filename):
    if not contents:
        return badge("No .rels", THEME["warn"]), None
    try:
        p = save_upload(contents, filename, UPLOADS_DIR)
        return html.Span([badge("Uploaded", THEME["good"]), f"  {p.name} ({p.stat().st_size:,} bytes)"]), str(p)
    except Exception as e:
        return badge(f"Upload failed: {e}", THEME["bad"]), None


@app.callback(
    Output("ents_upload_status", "children"),
    Output("ents_path_store", "data"),
    Input("ents_upload", "contents"),
    State("ents_upload", "filename"),
    prevent_initial_call=True,
)
def on_ents_upload(contents, filename):
    if not contents:
        return badge("No entities", THEME["warn"]), None
    try:
        p = save_upload(contents, filename, UPLOADS_DIR)
        return html.Span([badge("Uploaded", THEME["good"]), f"  {p.name} ({p.stat().st_size:,} bytes)"]), str(p)
    except Exception as e:
        return badge(f"Upload failed: {e}", THEME["bad"]), None


@app.callback(
    Output("gene_col_dd", "options"),
    Output("gene_col_dd", "value"),
    Output("entrez_col_dd", "options"),
    Output("entrez_col_dd", "value"),
    Output("fc_col_dd", "options"),
    Output("fc_col_dd", "value"),
    Output("pval_col_dd", "options"),
    Output("pval_col_dd", "value"),
    Output("direction_col_dd", "options"),
    Output("direction_col_dd", "value"),
    Input("expr_meta_store", "data"),
)
def populate_mapping_dropdowns(meta):
    opts = [{"label": "-- None --", "value": "__NONE__"}]
    defaults = ["__NONE__", "__NONE__", "__NONE__", "__NONE__", "__NONE__"]

    if meta and meta.get("columns"):
        columns = meta["columns"]
        opts = make_dropdown_options(columns)
        det = meta.get("detected", {})
        defaults = [
            det.get("gene_col", "__NONE__"),
            det.get("entrez_col", "__NONE__"),
            det.get("fc_col", "__NONE__"),
            det.get("pval_col", "__NONE__"),
            det.get("direction_col", "__NONE__"),
        ]

    return (
        opts, defaults[0],
        opts, defaults[1],
        opts, defaults[2],
        opts, defaults[3],
        opts, defaults[4],
    )


# =============================================================================
# Validation / preview
# =============================================================================

@app.callback(
    Output("signature_preview_block", "children"),
    Output("validation_summary_block", "children"),
    Output("normalized_preview_block", "children"),
    Output("validation_ok_store", "data"),
    Output("validation_inline_hint", "children"),
    Input("expr_meta_store", "data"),
    Input("engine_choice", "value"),
    Input("gene_col_dd", "value"),
    Input("entrez_col_dd", "value"),
    Input("fc_col_dd", "value"),
    Input("pval_col_dd", "value"),
    Input("direction_col_dd", "value"),
    Input("pval_thresh_in", "value"),
    Input("abs_fc_thresh_in", "value"),
    Input("use_pval_filter_ck", "value"),
    Input("use_abs_fc_filter_ck", "value"),
)
def validate_and_preview(
    meta,
    engine,
    gene_col,
    entrez_col,
    fc_col,
    pval_col,
    direction_col,
    pval_thresh,
    abs_fc_thresh,
    use_pval_filter_vals,
    use_abs_fc_filter_vals,
):
    empty_msg = html.Div("Upload a signature file to preview and validate it.", style={"color": THEME["muted"]})
    if not meta or not meta.get("raw_path"):
        return empty_msg, html.Div(), html.Div(), False, "Waiting for a signature file."

    raw_path = Path(meta["raw_path"])
    use_pval_filter = "USE_PVAL" in (use_pval_filter_vals or [])
    use_abs_fc_filter = "USE_ABS_FC" in (use_abs_fc_filter_vals or [])

    try:
        preview_df = pd.DataFrame(meta.get("preview_records", []))
        preview_block = card(
            "Uploaded signature preview",
            [
                html.Div(
                    style={"display": "grid", "gridTemplateColumns": "repeat(4, minmax(0,1fr))", "gap": "8px", "marginBottom": "10px"},
                    children=[
                        html.Div([badge(f"{meta.get('rows', 0):,} rows", THEME["blue"])]),
                        html.Div([badge(f"{meta.get('cols', 0)} columns", THEME["blue2"])]),
                        html.Div([badge(f"sep={meta.get('separator', ',')}", THEME["good"])]),
                        html.Div([badge(engine, THEME["warn"])]),
                    ],
                ),
                df_table(preview_df, page_size=10) if not preview_df.empty else html.Div("No preview rows available.", style={"color": THEME["muted"]}),
            ],
        )

        temp_out_dir = JOBS_DIR / "_preview_norm"
        temp_out_dir.mkdir(parents=True, exist_ok=True)

        res = normalize_signature_for_engine(
            raw_path=raw_path,
            engine=engine,
            out_dir=temp_out_dir,
            gene_col=gene_col,
            entrez_col=entrez_col,
            fc_col=fc_col,
            pval_col=pval_col,
            direction_col=direction_col,
            pval_thresh=float(pval_thresh if pval_thresh is not None else 0.05),
            abs_fc_thresh=float(abs_fc_thresh) if abs_fc_thresh not in (None, "") else None,
            use_pval_filter=use_pval_filter,
            use_abs_fc_filter=use_abs_fc_filter,
        )

        if not res["ok"]:
            children = [html.Div("Validation failed. Fix the mapping/filter settings below.", style={"fontWeight": "900", "color": THEME["bad"], "marginBottom": "10px"})]
            for e in res.get("errors", []):
                children.append(html.Div(f"* {e}", style={"color": THEME["bad"]}))
            return preview_block, card("Validation summary", children), html.Div(), False, "Validation failed. Review the selected columns and thresholds."

        stats = res["stats"]
        summary_block = card("Validation summary", [
            html.Div(
                style={"display": "grid", "gridTemplateColumns": "repeat(4, minmax(0,1fr))", "gap": "10px", "marginBottom": "12px"},
                children=[
                    stat_card("Rows loaded", f"{stats['rows_loaded']:,}", THEME["title"]),
                    stat_card("Rows after dedup", f"{stats['rows_after_dedup']:,}", THEME["blue"]),
                    stat_card("Rows after filters", f"{stats['rows_after_filter']:,}", THEME["good"]),
                    stat_card("Duplicates removed", f"{stats['duplicates_removed']:,}", THEME["warn"]),
                ],
            ),
            info_row("Engine", engine),
            info_row("Gene column", stats["gene_col"] or "None"),
            info_row("Entrez column", stats["entrez_col"] or "None"),
            info_row("Fold-change column", stats["fc_col"] or "None"),
            info_row("P-value column", stats["pval_col"] or "None"),
            info_row("Direction column", stats["direction_col"] or "None"),
            info_row("Sign source", stats["sign_source"]),
            info_row("P-value filter used", human_bool(stats["pval_filter_used"])),
            info_row("Absolute FC filter used", human_bool(stats["abs_fc_filter_used"])),
            info_row("Positive signs", f"{stats['positive_signs']:,}", THEME["good"]),
            info_row("Negative signs", f"{stats['negative_signs']:,}", THEME["bad"]),
            info_row("Zero signs", f"{stats['zero_signs']:,}"),
        ])

        norm_prev_df = pd.DataFrame(res.get("normalized_preview_records", []))
        normalized_preview_block = card(
            "Normalized signature preview",
            [
                html.Div("This is the deduplicated, filtered, normalized view that will be passed to the selected engine.", style={"fontSize": "12px", "color": THEME["muted"], "marginBottom": "10px"}),
                df_table(norm_prev_df, page_size=10) if not norm_prev_df.empty else html.Div("No normalized preview available.", style={"color": THEME["muted"]}),
            ],
        )

        hint = (
            f"Validation passed. {stats['rows_after_filter']:,} rows will be sent to {engine}. "
            f"Removed {stats['duplicates_removed']:,} duplicate rows before filtering."
        )
        return preview_block, summary_block, normalized_preview_block, True, hint

    except Exception as e:
        err = card("Validation summary", [html.Div(f"Validation crashed: {e}", style={"color": THEME["bad"], "fontWeight": "900"})])
        return empty_msg, err, html.Div(), False, f"Validation crashed: {e}"


# =============================================================================
# Run enablement
# =============================================================================

@app.callback(
    Output("run_btn", "style"),
    Output("run_hint", "children"),
    Input("validation_ok_store", "data"),
    Input("expr_path_store", "data"),
    Input("engine_choice", "value"),
    Input("network_choice", "value"),
    Input("rels_path_store", "data"),
    Input("ents_path_store", "data"),
)
def enable_run(validation_ok, expr_path, engine, network_choice, rels_uploaded, ents_uploaded):
    if not expr_path:
        return BTN_DISABLED, "Upload a signature file to enable Run."
    if not validation_ok:
        return BTN_DISABLED, "Fix validation issues before running."
    if network_choice == "__UPLOAD__":
        if not rels_uploaded:
            return BTN_DISABLED, "Upload a custom .rels network to enable Run."
        if engine == "CIE" and not ents_uploaded:
            return BTN_DISABLED, "Upload entities/ents for custom CIE runs."
    return BTN, "Ready. Click Run inference."


# =============================================================================
# Start run
# =============================================================================

@app.callback(
    Output("job_id_store", "data", allow_duplicate=True),
    Output("poll_interval", "disabled", allow_duplicate=True),
    Input("clear_results_btn", "n_clicks"),
    prevent_initial_call=True,
)
def clear_results(n):
    if not n:
        return no_update, no_update
    return None, True


@app.callback(
    Output("job_id_store", "data"),
    Output("poll_interval", "disabled"),
    Input("run_btn", "n_clicks"),
    State("engine_choice", "value"),
    State("expr_path_store", "data"),
    State("network_choice", "value"),
    State("rels_path_store", "data"),
    State("ents_path_store", "data"),
    State("gene_col_dd", "value"),
    State("entrez_col_dd", "value"),
    State("fc_col_dd", "value"),
    State("pval_col_dd", "value"),
    State("direction_col_dd", "value"),
    State("pval_thresh_in", "value"),
    State("abs_fc_thresh_in", "value"),
    State("use_pval_filter_ck", "value"),
    State("use_abs_fc_filter_ck", "value"),
    State("validation_ok_store", "data"),
    State("max_edges_in", "value"),
    State("skip_mcmc_ck", "value"),
    prevent_initial_call=True,
)
def start_run(
    n,
    engine,
    expr_path,
    network_choice,
    rels_uploaded,
    ents_uploaded,
    gene_col,
    entrez_col,
    fc_col,
    pval_col,
    direction_col,
    pval_thresh,
    abs_fc_thresh,
    use_pval_filter_vals,
    use_abs_fc_filter_vals,
    validation_ok,
    max_edges,
    skip_mcmc_vals,
):
    if not n or not expr_path or not validation_ok:
        return no_update, True

    if network_choice == "__UPLOAD__":
        if not rels_uploaded:
            return no_update, True
        if engine == "CIE" and not ents_uploaded:
            return no_update, True

    job_id = uuid.uuid4().hex[:10]
    job = _job_paths(job_id)
    job.job_dir.mkdir(parents=True, exist_ok=True)
    job.out_dir.mkdir(parents=True, exist_ok=True)

    use_pval_filter = "USE_PVAL" in (use_pval_filter_vals or [])
    use_abs_fc_filter = "USE_ABS_FC" in (use_abs_fc_filter_vals or [])
    skip_mcmc = "SKIP_MCMC" in (skip_mcmc_vals or [])

    norm_res = normalize_signature_for_engine(
        raw_path=Path(expr_path),
        engine=engine,
        out_dir=job.job_dir,
        gene_col=gene_col,
        entrez_col=entrez_col,
        fc_col=fc_col,
        pval_col=pval_col,
        direction_col=direction_col,
        pval_thresh=float(pval_thresh if pval_thresh is not None else 0.05),
        abs_fc_thresh=float(abs_fc_thresh) if abs_fc_thresh not in (None, "") else None,
        use_pval_filter=use_pval_filter,
        use_abs_fc_filter=use_abs_fc_filter,
    )

    if not norm_res["ok"]:
        write_status(job, state="error", progress=0, message="Normalization failed before run.", errors=norm_res.get("errors", []))
        return job_id, False

    threading.Thread(
        target=job_thread,
        kwargs=dict(
            engine=engine,
            normalized_input_path=norm_res["engine_input_path"],
            network_choice=network_choice,
            rels_uploaded=rels_uploaded,
            ents_uploaded=ents_uploaded,
            normalization_stats=norm_res["stats"],
            job_id=job_id,
            max_edges=int(max_edges) if max_edges is not None else 100,
            pval_thresh=float(pval_thresh if pval_thresh is not None else 0.05),
            log2fc_thresh=float(abs_fc_thresh if abs_fc_thresh not in (None, "") else 0.5),
            skip_mcmc=skip_mcmc,
            original_filename=Path(expr_path).stem if expr_path else "",
        ),
        daemon=True,
    ).start()

    return job_id, False


# =============================================================================
# Results tab memory
# =============================================================================

@app.callback(
    Output("results_tab_store", "data"),
    Input("results_tabs", "value"),
    prevent_initial_call=True,
)
def remember_selected_tab(tab_value):
    return tab_value or "summary"


# =============================================================================
# Poll / results
# =============================================================================

@app.callback(
    Output("progress_fill", "style"),
    Output("progress_text", "children"),
    Output("results_blocks", "children"),
    Output("results_badge", "children"),
    Output("results_message", "children"),
    Output("download_hint", "children"),
    Output("poll_interval", "disabled", allow_duplicate=True),
    Input("poll_interval", "n_intervals"),
    Input("job_id_store", "data"),
    Input("apply_display_btn", "n_clicks"),
    State("engine_choice", "value"),
    State("results_tab_store", "data"),
    State("top_n_regulators_in", "value"),
    State("use_pval_filter_check", "value"),
    State("min_reg_score_in", "value"),
    State("min_abs_edge_score_in", "value"),
    State("max_nodes_in", "value"),
    State("max_edges_in", "value"),
    State("label_top_hubs_in", "value"),
    prevent_initial_call="initial_duplicate",
)
def poll(
    n,
    job_id,
    apply_clicks,
    engine,
    current_tab,
    top_n_regulators,
    use_pval_filter,
    min_reg_score,
    min_abs_edge_score,
    max_nodes,
    max_edges,
    label_top_hubs,
):
    if not job_id:
        return (
            progress_fill_style(0),
            "",
            html.Div(),
            "",
            "Upload, map, validate, and then run.",
            "Run something first.",
            True,
        )

    top_n_regulators = max(1, safe_int(top_n_regulators, 20))
    min_reg_score = max(0.0, safe_float(min_reg_score, 0.0))
    min_abs_edge_score = max(0.0, safe_float(min_abs_edge_score, 0.0))
    max_nodes = max(4, safe_int(max_nodes, 40))
    max_edges = max(1, safe_int(max_edges, 100))
    label_top_hubs = max(0, safe_int(label_top_hubs, 12))

    job = _job_paths(job_id)
    st = read_status(job)

    # Always trust the engine recorded in status.json over the UI default
    engine = st.get("engine") or engine or "ORNOR"

    pct = clamp_int(st.get("progress", 0))
    state = st.get("state", "unknown")
    msg = st.get("message", "")
    fill = progress_fill_style(pct)

    print(f"[POLL] job_id={job_id} state={state} pct={pct} engine={engine}", flush=True)

    def _nonempty(fp: Path, min_bytes: int = 16) -> bool:
        try:
            return fp.exists() and fp.is_file() and fp.stat().st_size >= min_bytes
        except Exception:
            return False

    stdout_txt_now = ""
    if job.stdout_log.exists():
        try:
            stdout_txt_now = job.stdout_log.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            stdout_txt_now = ""

    cie_edges_csv = job.out_dir / "cie_edges.csv"
    cie_edges_tsv = job.out_dir / "cie_edges.tsv"

    ornor_tfs_csv = job.out_dir / "ornor_tfs.csv"
    ornor_tfs_tsv = job.out_dir / "ornor_tfs.tsv"
    ornor_edges_csv = job.out_dir / "ornor_edges.csv"
    ornor_edges_tsv = job.out_dir / "ornor_edges.tsv"

    cie_outputs_ready = _nonempty(cie_edges_csv) or _nonempty(cie_edges_tsv)
    ornor_outputs_ready = (
        (_nonempty(ornor_tfs_csv) or _nonempty(ornor_tfs_tsv))
        and
        (_nonempty(ornor_edges_csv) or _nonempty(ornor_edges_tsv))
    )

    outputs_ready = cie_outputs_ready or ornor_outputs_ready

    error_in_log = any(
        token in stdout_txt_now
        for token in [
            "Traceback (most recent call last):",
            "ModuleNotFoundError:",
            "ImportError:",
        ]
    )

    if state in ("queued", "running"):
        if outputs_ready:
            try:
                zip_outputs(job.job_dir / "out", job.zip_path, exclude_paths={job.zip_path})
            except Exception:
                pass

            write_status(
                job,
                state="done",
                progress=100,
                message="Finished.",
                zip_path=str(job.zip_path),
            )

            st = read_status(job)
            pct = clamp_int(st.get("progress", 100))
            state = st.get("state", "done")
            msg = st.get("message", "Finished.")
            fill = progress_fill_style(pct)

        elif error_in_log and not outputs_ready:
            write_status(
                job,
                state="error",
                progress=pct,
                message="Run failed. See stdout log.",
            )
            st = read_status(job)
            pct = clamp_int(st.get("progress", pct))
            state = st.get("state", "error")
            msg = st.get("message", "Run failed.")
            fill = progress_fill_style(pct)

    if state in ("queued", "running"):
        b = badge("Running", THEME["blue"])
        live_log = job.stdout_log.read_text(encoding="utf-8", errors="ignore")[-12000:] if job.stdout_log.exists() else ""
        blocks = card("Run in progress", [
            html.Div(msg, style={"color": THEME["muted"], "marginBottom": "8px"}),
            html.Details([
                html.Summary("Show live stdout"),
                html.Pre(live_log, style={"whiteSpace": "pre-wrap", "fontSize": "12px"})
            ]),
        ])
        return fill, f"{msg} ({pct}%)", blocks, b, f"Job {job_id}", "Download will appear when finished.", False

    if state == "error":
        b = badge("Error", THEME["bad"])
        log_txt = job.stdout_log.read_text(encoding="utf-8", errors="ignore") if job.stdout_log.exists() else "(no log)"
        errors = st.get("errors", [])
        blocks = card("Run failed", [
            html.Div(msg, style={"color": THEME["bad"], "fontWeight": "900", "marginBottom": "8px"}),
            *(html.Div(f"* {e}", style={"color": THEME["bad"]}) for e in errors),
            html.Details([
                html.Summary("Show stdout"),
                html.Pre(log_txt, style={"whiteSpace": "pre-wrap", "fontSize": "12px"})
            ]),
        ])
        return fill, msg, blocks, b, f"Job {job_id}", "No outputs to download.", True

    try:
     return _render_done(job_id, job, st, engine, current_tab, top_n_regulators, use_pval_filter, min_reg_score, min_abs_edge_score, max_nodes, max_edges, label_top_hubs)
    except Exception as _done_exc:
        import traceback as _tb
        print(f"[POLL DONE ERROR] {_done_exc}", flush=True)
        print(_tb.format_exc(), flush=True)
        raise

def _render_done(job_id, job, st, engine, current_tab, top_n_regulators, use_pval_filter, min_reg_score, min_abs_edge_score, max_nodes, max_edges, label_top_hubs):
    b = badge("Done", THEME["good"])
    out_dir = job.out_dir

    norm_json = job.job_dir / "normalization_summary.json"
    norm_stats = {}
    if norm_json.exists():
        try:
            norm_stats = json.loads(norm_json.read_text())
        except Exception:
            norm_stats = {}

    cie_tf_fp = _find_first_existing([
        out_dir / "cie_edges.csv_tfs.tsv",
        out_dir / "cie_edges.tsv_tfs.tsv",
        out_dir / "cie_tfs.tsv",
        out_dir / "cie_tfs.csv",
        out_dir / "cie_edges.csv_tfs.csv",
    ])
    cie_edge_fp = _find_first_existing([
        out_dir / "cie_edges.csv",
        out_dir / "cie_edges.tsv",
        out_dir / "edges.csv",
        out_dir / "edges.tsv",
    ])
    cie_pathway_fp = _find_first_existing([
        out_dir / "cie_edges.csv_pathwayEnrichment.tsv",
        out_dir / "cie_edges.tsv_pathwayEnrichment.tsv",
        out_dir / "pathwayEnrichment.tsv",
    ])

    ornor_tf_fp = _find_first_existing([out_dir / "ornor_tfs.csv", out_dir / "ornor_tfs.tsv"])
    ornor_edge_fp = _find_first_existing([out_dir / "ornor_edges.csv", out_dir / "ornor_edges.tsv"])

    tf_df_raw = read_result_table(cie_tf_fp if engine == "CIE" else ornor_tf_fp)
    edge_df_raw = read_result_table(cie_edge_fp if engine == "CIE" else ornor_edge_fp, nrows=5000)
    pathway_df = read_result_table(cie_pathway_fp if engine == "CIE" else None)

    run_meta_fp = _find_first_existing([out_dir / "cie_run_meta.json", job.job_dir / "run_request.json"])
    run_meta = {}
    if run_meta_fp and run_meta_fp.exists():
        try:
            run_meta = json.loads(run_meta_fp.read_text())
        except Exception:
            run_meta = {}

    ents_path = st.get("ents_path") or run_meta.get("ents_path")
    entities_name_map = load_entities_name_map(ents_path)

    tf_df = annotate_tf_table_with_names(tf_df_raw, entities_name_map)
    edge_df = annotate_edge_table_with_names(edge_df_raw, entities_name_map)

    stdout_txt = job.stdout_log.read_text(encoding="utf-8", errors="ignore") if job.stdout_log.exists() else ""
    ornor_quality = None
    if engine == "ORNOR":
        try:
            ornor_quality = detect_ornor_quality(
                tf_df=tf_df,
                edge_df=edge_df,
                stdout_txt=stdout_txt,
            )
        except Exception:
            ornor_quality = None

    if engine == "ORNOR" and ornor_quality is not None:
        b = html.Div(
            style={"display": "flex", "gap": "8px", "alignItems": "center", "flexWrap": "wrap"},
            children=[
                badge("Done", THEME["good"]),
                badge(str(ornor_quality["posterior_type"]), str(ornor_quality["posterior_color"])),
                badge(str(ornor_quality["convergence_label"]), str(ornor_quality["convergence_color"])),
            ],
        )

    tf_name_map = build_tf_name_map(tf_df, entities_name_map)

    # entities_name_map already covers all genes; if the edge table has a
    # dedicated trg_name column (e.g. custom networks), merge those in too.
    target_name_map = dict(entities_name_map)

    if edge_df is not None and not edge_df.empty:
        edge_cols = {str(c).strip().lower(): c for c in edge_df.columns}
        target_col = edge_cols.get("target") or edge_cols.get("trguid") or edge_cols.get("trg")
        trg_name_col = edge_cols.get("trg_name")

        if target_col and trg_name_col:
            # Build a tid→name dict from first occurrence per unique target (vectorized)
            sub = edge_df[[target_col, trg_name_col]].copy()
            sub[target_col] = sub[target_col].astype(str).str.strip()
            sub[trg_name_col] = sub[trg_name_col].astype(str).str.strip()
            sub = sub[sub[trg_name_col].str.lower().notna()]
            sub = sub.drop_duplicates(subset=[target_col])
            for tid, tname in zip(sub[target_col], sub[trg_name_col]):
                if tid and tname and tname.lower() not in {"nan", "none"} and tname != tid:
                    target_name_map[tid] = tname

    tf_display_df, tf_score_label = prepare_tf_display_df(
        tf_df=tf_df,
        engine=engine,
        top_n=top_n_regulators,
        min_display_score=min_reg_score,
        use_pval_filter="yes" in (use_pval_filter or []),
    )

    edge_display_df = prepare_edge_display_df(
        edge_df=edge_df,
        selected_tf_df=tf_display_df,
        min_abs_edge_score=min_abs_edge_score,
        max_edges=max_edges,
    )
    print(f"[DEBUG-EDGES] edge_df rows={len(edge_df)} edge_display_df rows={len(edge_display_df)} tf_display_df rows={len(tf_display_df)} tf_cols={list(tf_display_df.columns[:10])}")

    tf_plot = regulator_bar_plot(
        tf_df=tf_display_df,
        score_label=tf_score_label,
        title=("Top ORNOR regulators" if engine == "ORNOR" else "Top CIE regulators"),
        top_n=max(1, int(top_n_regulators or 20)),
    )
    _tf_id_col = coerce_first_existing_col(tf_display_df, ["TF", "tf", "source", "uid", "id"])
    _tf_score_map: Dict[str, float] = {}
    if _tf_id_col and not tf_display_df.empty and "display_score" in tf_display_df.columns:
        _tf_score_map = dict(zip(
            tf_display_df[_tf_id_col].astype(str).str.strip(),
            pd.to_numeric(tf_display_df["display_score"], errors="coerce").fillna(0.0),
        ))
    print(f"[DEBUG] _tf_id_col={_tf_id_col!r} score_map_keys={list(_tf_score_map.keys())[:5]} cols={list(tf_display_df.columns[:8])}")

    network_fig = build_regulatory_network_figure(
        edge_df=edge_display_df,
        tf_name_map=tf_name_map,
        target_name_map=target_name_map,
        title=("ORNOR Regulatory Network" if engine == "ORNOR" else "CIE Regulatory Network"),
        max_nodes=max(1, int(max_nodes or 40)),
        label_top_hubs=max(1, int(label_top_hubs or 12)),
        tf_score_map=_tf_score_map,
    )
    pathway_plot = pathway_bar_plot(pathway_df if engine == "CIE" else pd.DataFrame(), top_n=20)
    try:
        files_list = [html.Div(str(p.relative_to(job.job_dir))) for p in sorted(job.job_dir.rglob("*")) if p.is_file()]
    except Exception:
        files_list = []

    card_items = [
        ("Rows loaded", f"{norm_stats.get('rows_loaded', 0):,}", THEME["title"]),
        ("Rows after dedup", f"{norm_stats.get('rows_after_dedup', 0):,}", THEME["blue"]),
        ("Rows after filters", f"{norm_stats.get('rows_after_filter', 0):,}", THEME["good"]),
        ("Duplicates removed", f"{norm_stats.get('duplicates_removed', 0):,}", THEME["warn"]),
    ]
    if engine == "CIE":
        card_items.extend([
            ("Regulators", f"{len(tf_df):,}", THEME["purple"]),
            ("Edges", f"{len(edge_df):,}", THEME["blue"]),
            ("Pathways", f"{len(pathway_df):,}", THEME["warn"]),
            ("Positive signs", f"{norm_stats.get('positive_signs', 0):,}", THEME["good"]),
        ])
    else:
        card_items.extend([
            ("Regulators", f"{len(tf_df):,}", THEME["purple"]),
            ("Edges", f"{len(edge_df):,}", THEME["blue"]),
            ("Positive signs", f"{norm_stats.get('positive_signs', 0):,}", THEME["good"]),
            ("Negative signs", f"{norm_stats.get('negative_signs', 0):,}", THEME["bad"]),
        ])

    summary_cards = make_summary_cards({"cards": card_items})

    summary_children = []

    if engine == "ORNOR" and ornor_quality is not None:
        summary_children.append(render_ornor_quality_panel(ornor_quality))

    summary_children.append(
        card("Run summary", [
            summary_cards,
            info_row("Engine", engine),
            info_row("Network", str(st.get("rels_path", "")) if st.get("rels_path") else str(run_meta.get("rels_path", ""))),
            info_row("Rows loaded", f"{norm_stats.get('rows_loaded', 0):,}"),
            info_row("Rows after dedup", f"{norm_stats.get('rows_after_dedup', 0):,}", THEME["blue"]),
            info_row("Duplicates removed", f"{norm_stats.get('duplicates_removed', 0):,}", THEME["warn"]),
            info_row("Rows after filters", f"{norm_stats.get('rows_after_filter', 0):,}", THEME["good"]),
            info_row("P-value filter used", human_bool(norm_stats.get("pval_filter_used", False))),
            info_row("Absolute FC filter used", human_bool(norm_stats.get("abs_fc_filter_used", False))),
            info_row("Sign source", str(norm_stats.get("sign_source", ""))),
            html.Div(style={"height": "10px"}),
            html.Div("Display controls currently applied", style={"fontWeight": "900", "color": THEME["title"], "marginBottom": "8px"}),
            info_row("Top N regulators", str(top_n_regulators)),
            info_row("Minimum regulator score", f"{min_reg_score:g}"),
            info_row("Minimum |edge score|", f"{min_abs_edge_score:g}"),
            info_row("Maximum visible nodes", str(max_nodes)),
            info_row("Maximum visible edges", str(max_edges)),
            info_row("Label top hubs", str(label_top_hubs)),
            html.Div(style={"height": "10px"}),
            html.Div("Top regulator preview", style={"fontWeight": "900", "color": THEME["title"], "marginBottom": "8px"}),
            build_top_regulator_preview(tf_display_df, n=min(10, len(tf_display_df) if not tf_display_df.empty else 10)),
        ])
    )

    summary_details = html.Div(summary_children)

    tabs = dcc.Tabs(
        id="results_tabs",
        value=current_tab or "summary",
        colors={
            "border": THEME["panel_border"],
            "primary": THEME["blue"],
            "background": "#F8FBFF",
        },
        children=[
            dcc.Tab(
                label="Summary",
                value="summary",
                children=[
                    html.Div(style={"paddingTop": "14px"}, children=[summary_details]),
                ],
            ),
            dcc.Tab(
                label="TF / Regulator Results",
                value="tfs",
                children=[
                    html.Div(
                        style={"paddingTop": "14px"},
                        children=[
                            render_ornor_quality_strip(ornor_quality) if engine == "ORNOR" and ornor_quality is not None else html.Div(),
                            card(
                                "Regulator table",
                                [
                                    html.Div(
                                        f"Showing {len(tf_display_df):,} regulators after display filters. Labels include names + IDs when available.",
                                        style={"fontSize": "12px", "color": THEME["muted"], "marginBottom": "10px"},
                                    ),
                                    df_table(
                                        tf_display_df if not tf_display_df.empty else pd.DataFrame({"message": ["No regulator results available after filtering"]}),
                                        page_size=15,
                                    ),
                                ],
                            ),
                            card("Regulator plot", [dcc.Graph(figure=tf_plot, config={"displaylogo": False})]),
                        ],
                    )
                ],
            ),
            dcc.Tab(
                label="Network",
                value="network",
                children=[
                    html.Div(
                        style={"paddingTop": "14px"},
                        children=[
                            render_ornor_quality_strip(ornor_quality) if engine == "ORNOR" and ornor_quality is not None else html.Div(),
                            card("Network controls in effect", [
                                info_row("Minimum |edge score|", f"{min_abs_edge_score:g}"),
                                info_row("Maximum visible edges", str(max_edges)),
                                info_row("Maximum visible nodes", str(max_nodes)),
                                info_row("Label top hubs", str(label_top_hubs)),
                                html.Div(
                                    "The graph uses a bipartite layout with regulators on the left and targets on the right.",
                                    style={"color": THEME["muted"], "marginTop": "8px"},
                                ),
                            ]),
                            card("Interactive network graph", [dcc.Graph(figure=network_fig, config={"displaylogo": False})]),
                            card("Edge table", [
                                html.Div(
                                    f"Showing {len(edge_display_df):,} edges after display filters.",
                                    style={"fontSize": "12px", "color": THEME["muted"], "marginBottom": "10px"},
                                ),
                                df_table(
                                    edge_display_df.head(1000) if not edge_display_df.empty else pd.DataFrame({"message": ["No edge results available after filtering"]}),
                                    page_size=15,
                                )
                            ]),
                        ],
                    )
                ],
            ),
            dcc.Tab(
                label="Pathway Enrichment",
                value="pathways",
                children=[
                    html.Div(
                        style={"paddingTop": "14px"},
                        children=[
                            card("Pathway table", [df_table(pathway_df if not pathway_df.empty else pd.DataFrame({"message": ["No pathway enrichment output available"]}), page_size=15)]),
                            card("Pathway plot", [dcc.Graph(figure=pathway_plot, config={"displaylogo": False})]),
                        ],
                    )
                ],
            ) if engine == "CIE" else dcc.Tab(
                label="Pathway Enrichment",
                value="pathways",
                children=[
                    html.Div(
                        style={"paddingTop": "14px"},
                        children=[card("Pathway enrichment", [html.Div("Pathway enrichment is currently available for CIE outputs only.", style={"color": THEME["muted"]})])],
                    )
                ],
            ),
            dcc.Tab(
                label="Logs / Files",
                value="logs",
                children=[
                    html.Div(
                        style={"paddingTop": "14px"},
                        children=[
                            card("Logs / files", [
                                html.Details([html.Summary("Show stdout log"), html.Pre(stdout_txt[-15000:], style={"whiteSpace": "pre-wrap", "fontSize": "12px"})]),
                                html.Div(style={"height": "10px"}),
                                html.Details([html.Summary("Show job files"), html.Div(files_list)]),
                            ])
                        ],
                    )
                ],
            ),
        ],
    )

    blocks = html.Div([tabs])

    print(f"[RENDER_DONE] returning done results for job={job_id} engine={engine} tf_rows={len(tf_df) if tf_df is not None else 'None'}", flush=True)
    return progress_fill_style(100), "Done (100%)", blocks, b, f"Job {job_id}", "Outputs ready. Click Download outputs.", True


# =============================================================================
# Download
# =============================================================================

@app.callback(
    Output("download_outputs", "data"),
    Input("download_btn", "n_clicks"),
    State("job_id_store", "data"),
    prevent_initial_call=True,
)
def download_outputs_callback(n_clicks, job_id):
    if not job_id:
        return no_update

    job = _job_paths(job_id)
    st = read_status(job)
    if st.get("state") != "done":
        return no_update

    try:
        zip_outputs(job.job_dir / "out", job.zip_path, exclude_paths={job.zip_path})
        if not job.zip_path.exists():
            return no_update
        run_meta_fp = job.job_dir / "run_request.json"
        original_filename = ""
        if run_meta_fp.exists():
            try:
                original_filename = json.loads(run_meta_fp.read_text()).get("original_filename", "")
            except Exception:
                pass
        engine = st.get("engine", "inference")
        stem = original_filename or "results"
        download_name = f"{stem}_{engine.lower()}_outputs.zip"
        return dcc.send_file(str(job.zip_path), filename=download_name)
    except Exception:
        return no_update


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    # Pre-warm the entities cache for every known network so _render_done
    # never does expensive I/O on the hot polling path.
    _warmed: set = set()
    for _net in DEFAULT_NETWORKS.values():
        _ents = _net.get("entities")
        if _ents and str(_ents) not in _warmed and Path(_ents).exists():
            print(f"[STARTUP] Pre-warming entities cache: {_ents}", flush=True)
            load_entities_name_map(str(_ents))
            _warmed.add(str(_ents))
    print(f"[STARTUP] Entities cache ready ({len(_warmed)} file(s) loaded).", flush=True)

    # Use PORT env var if available (Railway, Heroku), otherwise default to 8050
    port = int(os.environ.get("PORT", 8050))
    app.run(host="0.0.0.0", port=port, debug=False, threaded=True)
