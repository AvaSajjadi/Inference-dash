# app1.py — Inference Dash research-tool UI
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
# - CLEARER NETWORK:
#   * strong bipartite layout
#   * fewer default edges/nodes
#   * bigger labels and nodes
#   * lighter edges
#   * much more visible regulator/target separation

import base64
import csv
import json
import re
import subprocess
import threading
import time
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
    "tcChIP (all_tissues.rels) [default]": {
        "rels": first_existing(
            NETWORKS_DIR / "all_tissues.rels",
        ),
        "entities": first_existing(
            NETWORKS_DIR / "all_tissues.fixed.entities",
            NETWORKS_DIR / "all_tissues.entities",
        ),
    },
    "tcChIP (three_tissues.rels)": {
        "rels": first_existing(
            NETWORKS_DIR / "three_tissues.rels",
            NETWORKS_DIR / "three_tissue.rels",
        ),
        "entities": first_existing(
            NETWORKS_DIR / "all_tissues.fixed.entities",
            NETWORKS_DIR / "all_tissues.entities",
        ),
    },
}

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
            for fp in out_dir.rglob("*"):
                if not fp.is_file():
                    continue
                try:
                    resolved = fp.resolve()
                except Exception:
                    resolved = fp

                if resolved == zip_self:
                    continue
                if resolved in exclude_resolved:
                    continue

                z.write(fp, arcname=str(fp.relative_to(out_dir)))

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


def load_entities_name_map(entities_path: Optional[str]) -> Dict[str, str]:
    if not entities_path:
        return {}

    p = Path(entities_path)
    if not p.exists():
        return {}

    df = _smart_read_table(p)
    if df.empty:
        return {}

    uid_col = _best_existing_col(df, ["uid", "id", "entrez", "entrezid", "gene_id"])
    name_col = _best_existing_col(df, ["name", "symbol", "gene", "gene_symbol", "label"])

    if not uid_col or not name_col:
        return {}

    out = {}
    for _, row in df[[uid_col, name_col]].dropna().iterrows():
        uid = str(row[uid_col]).strip()
        name = str(row[name_col]).strip()
        if uid and name and uid.lower() != "nan" and name.lower() != "nan":
            out[uid] = name
    return out


def lookup_name(uid, name_map: Dict[str, str]) -> str:
    s = str(uid).strip()
    if s in name_map:
        return name_map[s]
    try:
        f = float(s)
        if f.is_integer():
            s2 = str(int(f))
            if s2 in name_map:
                return name_map[s2]
    except Exception:
        pass
    return s


def annotate_tf_table_with_names(tf_df: pd.DataFrame, entities_name_map: Dict[str, str]) -> pd.DataFrame:
    if tf_df is None or tf_df.empty:
        return tf_df

    df = tf_df.copy()
    df.columns = [strip_bom_text(c) for c in df.columns]

    id_col = _best_existing_col(df, ["id", "uid", "tf_id"])
    name_col = _best_existing_col(df, ["name", "tf", "symbol"])

    if name_col:
        df["display_name"] = df[name_col].astype(str)
    elif id_col:
        df["display_name"] = df[id_col].astype(str).map(lambda x: lookup_name(x, entities_name_map))
    else:
        df["display_name"] = ""

    if id_col:
        df["display_label"] = df["display_name"].astype(str) + " (" + df[id_col].astype(str) + ")"
    else:
        df["display_label"] = df["display_name"].astype(str)

    preferred = []
    for c in ["display_name", "display_label"]:
        if c in df.columns:
            preferred.append(c)

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
        df["src_name"] = df[src_col].astype(str).map(lambda x: lookup_name(x, entities_name_map))
        df["src_label"] = df["src_name"].astype(str) + " (" + df[src_col].astype(str) + ")"

    if trg_col:
        df["trg_name"] = df[trg_col].astype(str).map(lambda x: lookup_name(x, entities_name_map))
        df["trg_label"] = df["trg_name"].astype(str) + " (" + df[trg_col].astype(str) + ")"

    ordered = []
    for c in ["src_name", "trg_name", "src_label", "trg_label"]:
        if c in df.columns:
            ordered.append(c)
    remaining = [c for c in df.columns if c not in ordered]
    return df[ordered + remaining]


def build_tf_name_map(tf_df: pd.DataFrame, entities_name_map: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    entities_name_map = entities_name_map or {}
    out = {}

    if tf_df is not None and not tf_df.empty:
        df = tf_df.copy()
        df.columns = [strip_bom_text(c) for c in df.columns]

        id_col = _best_existing_col(df, ["id", "uid", "tf_id"])
        name_col = _best_existing_col(df, ["display_name", "name", "tf", "symbol", "display_label"])

        if id_col and name_col:
            for _, row in df[[id_col, name_col]].dropna().iterrows():
                uid = str(row[id_col]).strip()
                name = str(row[name_col]).strip()
                if uid and name and uid.lower() != "nan" and name.lower() != "nan":
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

    norm = pd.DataFrame(index=df.index)

    if gene_col:
        norm["gene"] = df[gene_col].astype(str).str.strip()
        norm.loc[norm["gene"].str.lower().isin({"nan", "none", ""}), "gene"] = ""
    else:
        norm["gene"] = ""

    if entrez_col:
        norm["entrez"] = pd.to_numeric(df[entrez_col], errors="coerce")
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
            warnings.append("No p-value column provided; ORNOR input will use pval=1.0.")

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

    include = pd.Series(True, index=norm.index)

    if use_pval_filter and pval_col:
        include = include & norm["pval"].notna() & (norm["pval"] <= float(pval_thresh))

    if use_abs_fc_filter and fc_col and abs_fc_thresh is not None:
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
            "errors": ["No rows remain after applying identifier checks and selected filters."],
            "warnings": warnings,
        }

    canonical_path = out_dir / "canonical_signature.tsv"
    norm.to_csv(canonical_path, sep="\t", index=False)

    cie_input_path = out_dir / "cie_input.tsv"
    if engine == "CIE":
        cie_df = filtered.copy()
        cie_df["entrez"] = cie_df["entrez"].astype("Int64")
        cie_write = cie_df[["entrez", "fc", "pval"]].copy()
        cie_write.to_csv(cie_input_path, sep="\t", index=False)

    ornor_input_path = out_dir / "ornor_input.tsv"
    if engine == "ORNOR":
        ornor_df = filtered.copy()
        if entrez_col:
            ornor_df["entrez"] = ornor_df["entrez"].astype("Int64")
        ornor_write = ornor_df[["gene", "entrez", "direction", "fc", "pval", "sign"]].copy()
        ornor_write.to_csv(ornor_input_path, sep="\t", index=False)

    stats = {
        "rows_loaded": int(len(df)),
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

def run_with_progress(cmd: list, job: JobPaths) -> Tuple[int, str]:
    write_status(job, state="running", progress=0, message="Starting…", cmd=cmd)

    job.stdout_log.parent.mkdir(parents=True, exist_ok=True)
    last_msg = "Running…"

    with job.stdout_log.open("w", encoding="utf-8") as logf:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
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
                    write_status(job, progress=pct, message=f"Running… ({pct}%)")
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


def job_thread(
    engine: str,
    normalized_input_path: str,
    network_choice: str,
    rels_uploaded: Optional[str],
    ents_uploaded: Optional[str],
    normalization_stats: Dict,
    job_id: str,
):
    job = _job_paths(job_id)
    job.job_dir.mkdir(parents=True, exist_ok=True)
    job.out_dir.mkdir(parents=True, exist_ok=True)

    if network_choice == "__UPLOAD__":
        rels_path = rels_uploaded
        ents_path = ents_uploaded
    else:
        meta = DEFAULT_NETWORKS.get(network_choice, {})
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
        message=f"Queued ({engine})…",
        engine=engine,
        normalization_stats=normalization_stats,
        rels_path=rels_path,
        ents_path=ents_path,
    )

    (job.job_dir / "run_request.json").write_text(
        json.dumps(
            {
                "engine": engine,
                "normalized_input_path": str(in_path),
                "network_choice": network_choice,
                "rels_path": rels_path,
                "ents_path": ents_path,
                "normalization_stats": normalization_stats,
            },
            indent=2,
        )
    )

    if engine == "CIE":
        out_edges = job.out_dir / "cie_edges.csv"
        cmd = [
            "Rscript", str(CIE_RUNNER),
            "-s", str(in_path),
            "-o", str(out_edges),
            "--rels", str(rels_path),
            "--ents", str(ents_path),
            "--db", "tcChIP",
            "--tissue", "all",
            "-m", "Fisher",
            "-p", "0.05",
            "-f", "1.5",
            "-u", "1",
            "-c", "1",
        ]
        rc, last = run_with_progress(cmd, job)
        if rc != 0:
            write_status(job, state="error", message=f"CIE failed (rc={rc}). {last}")
            return
    else:
        out_tfs = job.out_dir / "ornor_tfs.csv"
        out_edges = job.out_dir / "ornor_edges.csv"
        cmd = [
            "python3", str(ORNOR_RUNNER),
            "--signature", str(in_path),
            "--network", str(rels_path),
            "--out_tfs", str(out_tfs),
            "--out_edges", str(out_edges),
        ]
        rc, last = run_with_progress(cmd, job)
        if rc != 0:
            write_status(job, state="error", message=f"ORNOR failed (rc={rc}). {last}")
            return

    zip_outputs(job.job_dir, job.zip_path, exclude_paths={job.zip_path})
    write_status(job, state="done", progress=100, message="Finished.", zip_path=str(job.zip_path))

# =============================================================================
# Result parsers / plots
# =============================================================================

def _find_first_existing(paths: List[Path]) -> Optional[Path]:
    return next((p for p in paths if p.exists()), None)


def read_result_table(fp: Optional[Path]) -> pd.DataFrame:
    if fp is None or not fp.exists():
        return pd.DataFrame()
    try:
        if fp.suffix.lower() == ".tsv":
            df = pd.read_csv(fp, sep="\t")
        else:
            df = pd.read_csv(fp)
        df.columns = [strip_bom_text(c) for c in df.columns]
        return df
    except Exception:
        return pd.DataFrame()


def coerce_first_existing_col(df: pd.DataFrame, names: List[str]) -> Optional[str]:
    lower_to_real = {str(c).strip().lower(): c for c in df.columns}
    for name in names:
        if name.lower() in lower_to_real:
            return lower_to_real[name.lower()]
    return None


def top_n_bar_plot(
    df: pd.DataFrame,
    label_col_candidates: List[str],
    score_col_candidates: List[str],
    title: str,
    transform_neglog10: bool = False,
    top_n: int = 20,
) -> go.Figure:
    fig = go.Figure()
    if df is None or df.empty:
        fig.update_layout(title=title, height=420)
        return fig

    label_col = coerce_first_existing_col(df, label_col_candidates)
    score_col = coerce_first_existing_col(df, score_col_candidates)

    if not label_col or not score_col:
        fig.update_layout(title=title, height=420)
        return fig

    plot_df = df[[label_col, score_col]].copy()
    plot_df[label_col] = plot_df[label_col].astype(str)
    plot_df[score_col] = pd.to_numeric(plot_df[score_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[score_col])

    if plot_df.empty:
        fig.update_layout(title=title, height=420)
        return fig

    y_title = score_col
    if transform_neglog10:
        plot_df[score_col] = plot_df[score_col].clip(lower=1e-300)
        plot_df["_plot_val"] = -np.log10(plot_df[score_col])
        y_col = "_plot_val"
        y_title = f"-log10({score_col})"
    else:
        y_col = score_col

    plot_df = plot_df.sort_values(y_col, ascending=False).head(top_n)

    fig.add_trace(
        go.Bar(
            x=plot_df[label_col].tolist(),
            y=plot_df[y_col].astype(float).tolist(),
        )
    )
    fig.update_layout(
        title=title,
        height=420,
        margin=dict(l=40, r=20, t=60, b=140),
        xaxis_title=label_col,
        yaxis_title=y_title,
    )
    fig.update_xaxes(tickangle=45)
    return fig


def build_cie_network_figure(
    edge_df: pd.DataFrame,
    tf_name_map: Optional[Dict[str, str]] = None,
    target_name_map: Optional[Dict[str, str]] = None,
    max_edges: int = 60,
    max_nodes: int = 30,
) -> go.Figure:
    """
    Very clear bipartite layout:
    - regulators on the left
    - targets on the right
    - stronger edges only
    - large labels
    - light edges
    """
    fig = go.Figure()
    tf_name_map = tf_name_map or {}
    target_name_map = target_name_map or {}

    if edge_df is None or edge_df.empty:
        fig.update_layout(
            title="CIE Regulatory Network",
            height=850,
            annotations=[dict(text="No edge data available", x=0.5, y=0.5, showarrow=False)],
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    df = edge_df.copy()
    df.columns = [strip_bom_text(c) for c in df.columns]
    lower = {c.lower(): c for c in df.columns}

    src_col = lower.get("srcuid") or lower.get("source")
    trg_col = lower.get("trguid") or lower.get("target")
    score_col = lower.get("score")

    if not src_col or not trg_col:
        fig.update_layout(
            title="CIE Regulatory Network",
            height=850,
            annotations=[dict(text="Edge file is missing srcuid/trguid columns", x=0.5, y=0.5, showarrow=False)],
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    df["_src_uid"] = df[src_col].astype(str)
    df["_trg_uid"] = df[trg_col].astype(str)
    if score_col:
        df["_score"] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
    else:
        df["_score"] = 0.0

    # Prefer strongest edges
    df["_abs_score"] = df["_score"].abs()
    df = df.sort_values("_abs_score", ascending=False).head(max_edges).copy()

    if df.empty:
        fig.update_layout(
            title="CIE Regulatory Network",
            height=850,
            annotations=[dict(text="No edges available after filtering", x=0.5, y=0.5, showarrow=False)],
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    # Node counts
    reg_counts = df["_src_uid"].value_counts().to_dict()
    tgt_counts = df["_trg_uid"].value_counts().to_dict()

    # Keep a manageable number of visible nodes
    max_regs = max(1, max_nodes // 2)
    max_tgts = max(1, max_nodes // 2)

    regs = (
        df["_src_uid"].value_counts()
        .sort_values(ascending=False)
        .head(max_regs)
        .index.astype(str)
        .tolist()
    )

    tgts = (
        df["_trg_uid"].value_counts()
        .sort_values(ascending=False)
        .head(max_tgts)
        .index.astype(str)
        .tolist()
    )

    df = df[df["_src_uid"].isin(regs) & df["_trg_uid"].isin(tgts)].copy()

    regs = list(dict.fromkeys(df["_src_uid"].astype(str).tolist()))
    tgts = list(dict.fromkeys(df["_trg_uid"].astype(str).tolist()))

    if not regs or not tgts:
        fig.update_layout(
            title="CIE Regulatory Network",
            height=850,
            annotations=[dict(text="No nodes available after filtering", x=0.5, y=0.5, showarrow=False)],
            plot_bgcolor="white",
            paper_bgcolor="white",
        )
        return fig

    # Sort by degree so the important nodes appear near the center/top
    regs = sorted(regs, key=lambda x: (-reg_counts.get(x, 0), str(tf_name_map.get(x, x))))
    tgts = sorted(tgts, key=lambda x: (-tgt_counts.get(x, 0), str(target_name_map.get(x, x))))

    def reg_label(uid: str) -> str:
        return str(tf_name_map.get(str(uid), str(uid)))

    def tgt_label(uid: str) -> str:
        return str(target_name_map.get(str(uid), str(uid)))

    # -------------------------------------------------------------------------
    # Strong, easy-to-read bipartite layout
    # -------------------------------------------------------------------------
    pos = {}

    left_x = -2.0
    right_x = 2.0

    reg_y = np.linspace(2.0, -2.0, len(regs))
    tgt_y = np.linspace(2.0, -2.0, len(tgts))

    for i, r in enumerate(regs):
        pos[r] = (left_x, float(reg_y[i]))

    for i, t in enumerate(tgts):
        pos[t] = (right_x, float(tgt_y[i]))

    # -------------------------------------------------------------------------
    # Edges by sign
    # -------------------------------------------------------------------------
    pos_edge_x, pos_edge_y, pos_hover = [], [], []
    neg_edge_x, neg_edge_y, neg_hover = [], [], []
    zero_edge_x, zero_edge_y, zero_hover = [], [], []

    for _, row in df.iterrows():
        s = str(row["_src_uid"])
        t = str(row["_trg_uid"])
        score = float(row["_score"])

        if s not in pos or t not in pos:
            continue

        x0, y0 = pos[s]
        x1, y1 = pos[t]

        hover = f"{reg_label(s)} → {tgt_label(t)}<br>score={score:g}"

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
                line=dict(width=2.2, color="rgba(86, 193, 108, 0.42)"),
                hoverinfo="text",
                text=pos_hover,
                name="agreement (+1)",
            )
        )

    if neg_edge_x:
        fig.add_trace(
            go.Scatter(
                x=neg_edge_x,
                y=neg_edge_y,
                mode="lines",
                line=dict(width=2.0, color="rgba(234, 108, 108, 0.36)"),
                hoverinfo="text",
                text=neg_hover,
                name="disagreement (-1)",
            )
        )

    if zero_edge_x:
        fig.add_trace(
            go.Scatter(
                x=zero_edge_x,
                y=zero_edge_y,
                mode="lines",
                line=dict(width=1.1, color="rgba(160, 174, 192, 0.20)"),
                hoverinfo="text",
                text=zero_hover,
                name="neutral (0)",
            )
        )

    # -------------------------------------------------------------------------
    # Nodes
    # -------------------------------------------------------------------------
    reg_sizes = [28 + 3 * min(reg_counts.get(n, 1), 8) for n in regs]
    tgt_sizes = [18 + 2 * min(tgt_counts.get(n, 1), 8) for n in tgts]

    fig.add_trace(
        go.Scatter(
            x=[pos[n][0] for n in regs],
            y=[pos[n][1] for n in regs],
            mode="markers+text",
            marker=dict(
                size=reg_sizes,
                color="#ff4d6d",
                line=dict(width=2.2, color="white"),
            ),
            text=[reg_label(n) for n in regs],
            textposition="middle left",
            textfont=dict(size=18, color=THEME["title"]),
            hovertemplate=[
                f"regulator={reg_label(n)}<br>uid={n}<br>degree={reg_counts.get(n, 0)}<extra></extra>"
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
                line=dict(width=2.0, color="white"),
            ),
            text=[tgt_label(n) for n in tgts],
            textposition="middle right",
            textfont=dict(size=16, color=THEME["title"]),
            hovertemplate=[
                f"target={tgt_label(n)}<br>uid={n}<br>degree={tgt_counts.get(n, 0)}<extra></extra>"
                for n in tgts
            ],
            name="targets",
        )
    )

    fig.update_layout(
        title=f"CIE Regulatory Network (top_edges={max_edges}, max_nodes={max_nodes})",
        height=900,
        margin=dict(l=80, r=80, t=70, b=40),
        xaxis=dict(visible=False, range=[-3.2, 3.2]),
        yaxis=dict(visible=False, range=[-2.5, 2.5]),
        showlegend=True,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(
            x=0.86,
            y=0.98,
            bgcolor="rgba(255,255,255,0.88)",
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


def build_top_regulator_preview(tf_df: pd.DataFrame, engine: str, n: int = 10) -> dash_table.DataTable:
    if tf_df is None or tf_df.empty:
        return df_table(pd.DataFrame({"message": ["No regulator table available"]}), page_size=5)

    df = tf_df.copy()
    label_col = coerce_first_existing_col(df, ["display_name", "display_label", "name", "tf", "uid", "id"])
    if engine == "CIE":
        score_col = coerce_first_existing_col(df, ["pvalue", "pval", "fdr"])
        if score_col:
            df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
            df = df.sort_values(score_col, ascending=True)
    else:
        score_col = coerce_first_existing_col(df, ["posterior", "posterior_mean", "mean", "score"])
        if score_col:
            df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
            df = df.sort_values(score_col, ascending=False)

    preview_cols = [
        c for c in [
            label_col,
            score_col,
            coerce_first_existing_col(df, ["fdr"]),
            coerce_first_existing_col(df, ["proteinsfound"]),
        ] if c
    ]
    preview_df = df[preview_cols].head(n).copy() if preview_cols else df.head(n).copy()
    return df_table(preview_df, page_size=n)

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
                "width": "410px",
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
                    dcc.Dropdown(id="gene_col_dd", options=[{"label": "— None —", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("Entrez column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="entrez_col_dd", options=[{"label": "— None —", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("Fold-change column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="fc_col_dd", options=[{"label": "— None —", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("P-value / adjusted p-value column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="pval_col_dd", options=[{"label": "— None —", "value": "__NONE__"}], value="__NONE__", clearable=False),

                    html.Div(style={"height": "10px"}),
                    html.Div("Direction column", style={"fontWeight": "800", "fontSize": "12px", "color": THEME["muted"]}),
                    dcc.Dropdown(id="direction_col_dd", options=[{"label": "— None —", "value": "__NONE__"}], value="__NONE__", clearable=False),
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

                    html.Div(id="validation_inline_hint", style={"marginTop": "10px", "fontSize": "12px", "color": THEME["muted"]}),
                ]),

                card("5) Network / database", [
                    dcc.Dropdown(
                        id="network_choice",
                        options=[{"label": k, "value": k} for k in DEFAULT_NETWORKS.keys()] + [{"label": "Upload custom network…", "value": "__UPLOAD__"}],
                        value="tcChIP (all_tissues.rels) [default]" if DEFAULT_NETWORKS else "__UPLOAD__",
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
            ],
        ),

        html.Div(
            style={"flex": "1", "padding": "18px 22px", "overflowY": "auto"},
            children=[
                dcc.Interval(id="poll_interval", interval=1200, n_intervals=0, disabled=True),

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
            f"  {p.name} ({p.stat().st_size:,} bytes) — {meta['rows']:,} rows, {meta['cols']} columns",
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
    opts = [{"label": "— None —", "value": "__NONE__"}]
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
                children.append(html.Div(f"• {e}", style={"color": THEME["bad"]}))
            return preview_block, card("Validation summary", children), html.Div(), False, "Validation failed. Review the selected columns and thresholds."

        stats = res["stats"]
        summary_block = card("Validation summary", [
            html.Div(
                style={"display": "grid", "gridTemplateColumns": "repeat(4, minmax(0,1fr))", "gap": "10px", "marginBottom": "12px"},
                children=[
                    stat_card("Rows loaded", f"{stats['rows_loaded']:,}", THEME["title"]),
                    stat_card("Rows after filters", f"{stats['rows_after_filter']:,}", THEME["good"]),
                    stat_card("Positive signs", f"{stats['positive_signs']:,}", THEME["good"]),
                    stat_card("Negative signs", f"{stats['negative_signs']:,}", THEME["bad"]),
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
            info_row("Zero signs", f"{stats['zero_signs']:,}"),
        ])

        norm_prev_df = pd.DataFrame(res.get("normalized_preview_records", []))
        normalized_preview_block = card(
            "Normalized signature preview",
            [
                html.Div("This is the filtered / normalized view that will be passed to the selected engine.", style={"fontSize": "12px", "color": THEME["muted"], "marginBottom": "10px"}),
                df_table(norm_prev_df, page_size=10) if not norm_prev_df.empty else html.Div("No normalized preview available.", style={"color": THEME["muted"]}),
            ],
        )

        hint = f"Validation passed. {stats['rows_after_filter']:,} rows will be sent to {engine}."
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
    State("job_id_store", "data"),
    State("engine_choice", "value"),
    State("results_tab_store", "data"),
    prevent_initial_call=True,
)
def poll(n, job_id, engine, current_tab):
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

    job = _job_paths(job_id)
    st = read_status(job)

    pct = clamp_int(st.get("progress", 0))
    state = st.get("state", "unknown")
    msg = st.get("message", "")
    fill = progress_fill_style(pct)

    if state in ("queued", "running"):
        b = badge("Running", THEME["blue"])
        live_log = job.stdout_log.read_text(encoding="utf-8", errors="ignore")[-12000:] if job.stdout_log.exists() else ""
        blocks = card("Run in progress", [
            html.Div(msg, style={"color": THEME["muted"], "marginBottom": "8px"}),
            html.Details([html.Summary("Show live stdout"), html.Pre(live_log, style={"whiteSpace": "pre-wrap", "fontSize": "12px"})]),
        ])
        return fill, f"{msg} ({pct}%)", blocks, b, f"Job {job_id}", "Download will appear when finished.", False

    if state == "error":
        b = badge("Error", THEME["bad"])
        log_txt = job.stdout_log.read_text(encoding="utf-8", errors="ignore") if job.stdout_log.exists() else "(no log)"
        errors = st.get("errors", [])
        blocks = card("Run failed", [
            html.Div(msg, style={"color": THEME["bad"], "fontWeight": "900", "marginBottom": "8px"}),
            *(html.Div(f"• {e}", style={"color": THEME["bad"]}) for e in errors),
            html.Details([html.Summary("Show stdout"), html.Pre(log_txt, style={"whiteSpace": "pre-wrap", "fontSize": "12px"})]),
        ])
        return fill, msg, blocks, b, f"Job {job_id}", "No outputs to download.", True

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
    edge_df_raw = read_result_table(cie_edge_fp if engine == "CIE" else ornor_edge_fp)
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

    tf_name_map = build_tf_name_map(tf_df, entities_name_map)
    target_name_map = entities_name_map

    card_items = [
        ("Rows loaded", f"{norm_stats.get('rows_loaded', 0):,}", THEME["title"]),
        ("Rows after filters", f"{norm_stats.get('rows_after_filter', 0):,}", THEME["good"]),
        ("Regulators", f"{len(tf_df):,}", THEME["purple"]),
        ("Edges", f"{len(edge_df):,}", THEME["blue"]),
    ]
    if engine == "CIE":
        card_items.extend([
            ("Pathways", f"{len(pathway_df):,}", THEME["warn"]),
            ("Positive signs", f"{norm_stats.get('positive_signs', 0):,}", THEME["good"]),
            ("Negative signs", f"{norm_stats.get('negative_signs', 0):,}", THEME["bad"]),
            ("Zero signs", f"{norm_stats.get('zero_signs', 0):,}", THEME["muted"]),
        ])
    else:
        card_items.extend([
            ("Positive signs", f"{norm_stats.get('positive_signs', 0):,}", THEME["good"]),
            ("Negative signs", f"{norm_stats.get('negative_signs', 0):,}", THEME["bad"]),
            ("Zero signs", f"{norm_stats.get('zero_signs', 0):,}", THEME["muted"]),
            ("Engine", engine, THEME["title"]),
        ])

    summary_cards = make_summary_cards({"cards": card_items})

    summary_details = card("Run summary", [
        summary_cards,
        info_row("Engine", engine),
        info_row("Network", str(st.get("rels_path", "")) if st.get("rels_path") else str(run_meta.get("rels_path", ""))),
        info_row("Rows loaded", f"{norm_stats.get('rows_loaded', 0):,}"),
        info_row("Rows after filters", f"{norm_stats.get('rows_after_filter', 0):,}", THEME["good"]),
        info_row("P-value filter used", human_bool(norm_stats.get("pval_filter_used", False))),
        info_row("Absolute FC filter used", human_bool(norm_stats.get("abs_fc_filter_used", False))),
        info_row("Sign source", str(norm_stats.get("sign_source", ""))),
        html.Div(style={"height": "10px"}),
        html.Div("Top regulator preview", style={"fontWeight": "900", "color": THEME["title"], "marginBottom": "8px"}),
        build_top_regulator_preview(tf_df, engine, n=10),
    ])

    tf_plot = top_n_bar_plot(
        tf_df,
        label_col_candidates=["display_name", "display_label", "name", "tf", "uid", "id"],
        score_col_candidates=["pvalue", "pval", "posterior", "posterior_mean", "mean", "score"],
        title="Top regulators",
        transform_neglog10=(engine == "CIE"),
        top_n=20,
    )

    pathway_plot = top_n_bar_plot(
        pathway_df,
        label_col_candidates=["pathway", "name", "term", "id"],
        score_col_candidates=["pvalue", "pval", "fdr"],
        title="Top pathways",
        transform_neglog10=True,
        top_n=20,
    )

    stdout_txt = job.stdout_log.read_text(encoding="utf-8", errors="ignore") if job.stdout_log.exists() else "(no stdout log)"
    files_list = [
        html.Div(str(fp.relative_to(job.job_dir)), style={"fontFamily": "ui-monospace, monospace", "fontSize": "12px"})
        for fp in sorted(job.job_dir.rglob("*"))
        if fp.is_file()
    ]

    if engine == "CIE":
        network_fig = build_cie_network_figure(
            edge_df=edge_df_raw,
            tf_name_map=tf_name_map,
            target_name_map=target_name_map,
            max_edges=60,
            max_nodes=30,
        )
        network_graph_card = card(
            "Interactive network graph",
            [dcc.Graph(figure=network_fig, config={"displaylogo": False})],
        )
        network_controls = card("Network controls", [
            html.Div(
                "High-visibility defaults applied: top_edges=60, max_nodes=30, large labels, lighter edges, strong left/right separation.",
                style={"color": THEME["muted"]},
            )
        ])
    else:
        network_controls = html.Div()
        network_graph_card = card(
            "Network graph",
            [html.Div("Network graph is currently available for CIE outputs only.", style={"color": THEME["muted"]})],
        )

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
                            card("Regulator table", [df_table(tf_df if not tf_df.empty else pd.DataFrame({"message": ["No regulator results available"]}), page_size=15)]),
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
                            network_controls,
                            network_graph_card,
                            card("Edge table", [
                                df_table(
                                    edge_df.head(500) if not edge_df.empty else pd.DataFrame({"message": ["No edge results available"]}),
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

    return progress_fill_style(100), "Done (100%)", blocks, b, f"Job {job_id}", "✅ Outputs ready. Click Download outputs.", True

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
        zip_outputs(job.job_dir, job.zip_path, exclude_paths={job.zip_path})
        if not job.zip_path.exists():
            return no_update
        return dcc.send_file(str(job.zip_path))
    except Exception:
        return no_update

# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8050, debug=False, threaded=True)
