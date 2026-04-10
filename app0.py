# app.py
# Inference Dash — single-file Dash app
#
# Improvements:
# - Professional readable network plot:
#   * prune aggressively
#   * keep largest connected component
#   * increase spacing in spring layout
#   * hover for all labels, only label top hubs
#   * reduce edges shown

import os
import re
import uuid
import time
import base64
import subprocess
import math
from pathlib import Path

import numpy as np
import pandas as pd

from dash import Dash, dcc, html, dash_table, Input, Output, State, no_update
import plotly.graph_objects as go

try:
    import networkx as nx
except Exception:
    nx = None

# -----------------------------
# Paths / folders
# -----------------------------
APP_ROOT = Path(__file__).resolve().parent
UPLOADS_DIR = APP_ROOT / "uploads"
RESULTS_DIR = APP_ROOT / "results"
NETWORKS_DIR = APP_ROOT / "networks"
BACKEND_DIR = APP_ROOT / "backend"

UPLOADS_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
NETWORKS_DIR.mkdir(parents=True, exist_ok=True)

DEFAULT_NETWORK_EXT = ".rels"

CIE_R_SCRIPT = BACKEND_DIR / "engines" / "cie" / "run_cie.R"
ORNOR_PY_SCRIPT = BACKEND_DIR / "engines" / "ornor" / "run_ornor_real.py"

ENTREZ_TO_SYMBOL_PATH = UPLOADS_DIR / "entrez_to_symbol.tsv"


# -----------------------------
# Helpers
# -----------------------------
def _safe_name(name: str) -> str:
    name = os.path.basename(name or "uploaded.tsv")
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name[:180]


def save_upload(contents: str, filename: str) -> Path:
    if not contents:
        raise ValueError("No contents to save.")
    _, content_string = contents.split(",", 1)
    raw = base64.b64decode(content_string.encode("utf-8"))
    suffix = Path(_safe_name(filename)).suffix or ".tsv"
    stem = Path(_safe_name(filename)).stem or "uploaded"
    out = UPLOADS_DIR / f"uploaded_{uuid.uuid4().hex[:12]}_{stem}{suffix}"
    out.write_bytes(raw)
    return out


def list_default_networks():
    return sorted([p.name for p in NETWORKS_DIR.glob(f"*{DEFAULT_NETWORK_EXT}")])


def load_entrez_symbol_map() -> dict:
    if not ENTREZ_TO_SYMBOL_PATH.exists():
        return {}
    try:
        m = pd.read_csv(ENTREZ_TO_SYMBOL_PATH, sep="\t")
        if "ENTREZID" not in m.columns or "SYMBOL" not in m.columns:
            return {}
        m = m.dropna(subset=["ENTREZID", "SYMBOL"]).copy()
        m["ENTREZID"] = pd.to_numeric(m["ENTREZID"], errors="coerce").astype("Int64")
        m = m.dropna(subset=["ENTREZID"])
        return dict(zip(m["ENTREZID"].astype(int).tolist(), m["SYMBOL"].astype(str).tolist()))
    except Exception:
        return {}


def maybe_add_symbols(df: pd.DataFrame, col: str, mapping: dict) -> pd.Series:
    if col not in df.columns:
        return pd.Series([""] * len(df))
    vals = df[col].copy()
    ids = pd.to_numeric(vals, errors="coerce")
    labels = []
    for v, i in zip(vals.tolist(), ids.tolist()):
        if i is not None and not (isinstance(i, float) and np.isnan(i)):
            ii = int(i)
            sym = mapping.get(ii)
            labels.append(f"{sym} ({ii})" if sym else str(ii))
        else:
            labels.append(str(v))
    return pd.Series(labels)


def read_edges_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def build_score_hist(df: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(title="Edge score distribution", height=320, margin=dict(l=30, r=20, t=50, b=30))
    if df is None or df.empty or "score" not in df.columns:
        return fig
    scores = pd.to_numeric(df["score"], errors="coerce").dropna().values
    fig.add_trace(go.Histogram(x=scores, nbinsx=30))
    fig.update_layout(xaxis_title="score", yaxis_title="count")
    return fig


def build_bipartite_network(df: pd.DataFrame, mapping: dict, top_edges: int = 120) -> go.Figure:
    """Fallback if networkx missing."""
    fig = go.Figure()
    fig.update_layout(
        title=f"Network view (bipartite fallback, top {top_edges} edges)",
        height=520,
        margin=dict(l=30, r=20, t=50, b=30),
        showlegend=True,
    )
    if df is None or df.empty or "source" not in df.columns or "target" not in df.columns:
        return fig

    work = df.copy()
    if "score" in work.columns:
        work["_abs_score"] = pd.to_numeric(work["score"], errors="coerce").abs()
    else:
        work["_abs_score"] = 0.0

    if "padj" in work.columns:
        work["_p"] = pd.to_numeric(work["padj"], errors="coerce")
    elif "pvalue" in work.columns:
        work["_p"] = pd.to_numeric(work["pvalue"], errors="coerce")
    else:
        work["_p"] = np.nan

    work = work.sort_values(by=["_abs_score", "_p"], ascending=[False, True]).head(int(top_edges)).copy()
    src_lbl = maybe_add_symbols(work, "source", mapping)
    tgt_lbl = maybe_add_symbols(work, "target", mapping)

    src_nodes = sorted(set(src_lbl.tolist()))
    tgt_nodes = sorted(set(tgt_lbl.tolist()))

    x_src, x_tgt = 0.0, 1.0
    y_src = np.linspace(0, 1, len(src_nodes)) if src_nodes else np.array([])
    y_tgt = np.linspace(0, 1, len(tgt_nodes)) if tgt_nodes else np.array([])

    src_pos = {n: (x_src, float(y)) for n, y in zip(src_nodes, y_src)}
    tgt_pos = {n: (x_tgt, float(y)) for n, y in zip(tgt_nodes, y_tgt)}

    edge_x, edge_y = [], []
    for s, t in zip(src_lbl.tolist(), tgt_lbl.tolist()):
        xs, ys = src_pos[s]
        xt, yt = tgt_pos[t]
        edge_x += [xs, xt, None]
        edge_y += [ys, yt, None]

    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines", name="edges", hoverinfo="none", line=dict(width=1), opacity=0.20))
    fig.add_trace(go.Scatter(
        x=[src_pos[n][0] for n in src_nodes],
        y=[src_pos[n][1] for n in src_nodes],
        mode="markers",
        name="sources",
        marker=dict(size=8),
        hovertext=src_nodes,
        hoverinfo="text",
    ))
    fig.add_trace(go.Scatter(
        x=[tgt_pos[n][0] for n in tgt_nodes],
        y=[tgt_pos[n][1] for n in tgt_nodes],
        mode="markers",
        name="targets",
        marker=dict(size=8),
        hovertext=tgt_nodes,
        hoverinfo="text",
    ))
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(plot_bgcolor="white")
    return fig


def build_clean_network(
    df: pd.DataFrame,
    mapping: dict,
    top_edges: int = 160,
    max_nodes: int = 45,
    label_top_k: int = 10,
    min_abs_score_quantile: float = 0.70,
) -> go.Figure:
    """
    Actually readable network:
    - Take top_edges by abs(score), then apply an abs(score) quantile cutoff (removes weak edges).
    - Build graph, prune nodes to max_nodes by degree.
    - Keep ONLY the largest connected component.
    - Spring layout with more spacing.
    - No labels except top hubs; hover shows everything.
    """
    fig = go.Figure()
    fig.update_layout(
        title=f"Network view (readable, top_edges={top_edges}, max_nodes={max_nodes})",
        height=620,
        margin=dict(l=20, r=20, t=50, b=20),
        showlegend=True,
        hovermode="closest",
    )

    if df is None or df.empty or "source" not in df.columns or "target" not in df.columns:
        return fig

    if nx is None:
        return build_bipartite_network(df, mapping=mapping, top_edges=min(120, int(top_edges)))

    work = df.copy()

    # abs(score)
    if "score" in work.columns:
        work["_abs_score"] = pd.to_numeric(work["score"], errors="coerce").abs()
    else:
        work["_abs_score"] = 0.0

    # p preference
    if "padj" in work.columns:
        work["_p"] = pd.to_numeric(work["padj"], errors="coerce")
    elif "pvalue" in work.columns:
        work["_p"] = pd.to_numeric(work["pvalue"], errors="coerce")
    else:
        work["_p"] = np.nan

    # pick strongest edges
    work = work.sort_values(by=["_abs_score", "_p"], ascending=[False, True]).head(int(top_edges)).copy()

    # remove weak edges (quantile threshold)
    abs_vals = work["_abs_score"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(abs_vals) > 5:
        thr = float(abs_vals.quantile(min_abs_score_quantile))
        work = work[work["_abs_score"] >= thr].copy()

    if work.empty:
        return fig

    # Labels
    src_lbl = maybe_add_symbols(work, "source", mapping)
    tgt_lbl = maybe_add_symbols(work, "target", mapping)

    # Build graph
    G = nx.Graph()
    for i, (s, t) in enumerate(zip(src_lbl.tolist(), tgt_lbl.tolist())):
        score = work["score"].iloc[i] if "score" in work.columns else None
        pval = work["pvalue"].iloc[i] if "pvalue" in work.columns else None
        padj = work["padj"].iloc[i] if "padj" in work.columns else None
        abs_score = float(work["_abs_score"].iloc[i]) if "_abs_score" in work.columns else 0.0

        G.add_node(s, side="source")
        G.add_node(t, side="target")

        # keep max abs edge if duplicates
        if G.has_edge(s, t):
            if abs_score > float(G[s][t].get("abs_score", -1)):
                G[s][t].update(score=score, abs_score=abs_score, pvalue=pval, padj=padj)
        else:
            G.add_edge(s, t, score=score, abs_score=abs_score, pvalue=pval, padj=padj)

    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return fig

    # prune by degree
    if G.number_of_nodes() > int(max_nodes):
        deg = dict(G.degree())
        keep = sorted(deg, key=lambda n: deg[n], reverse=True)[: int(max_nodes)]
        G = G.subgraph(keep).copy()

    # keep largest connected component only (critical for readability)
    if G.number_of_nodes() > 0:
        comps = list(nx.connected_components(G))
        if comps:
            largest = max(comps, key=len)
            G = G.subgraph(largest).copy()

    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return fig

    # Layout: more spacing
    n = G.number_of_nodes()
    k = 2.2 / math.sqrt(max(n, 1))  # bigger => more spread
    pos = nx.spring_layout(G, seed=7, k=k, iterations=350)

    # Node sizes by degree
    deg = dict(G.degree())
    max_deg = max(deg.values()) if deg else 1

    nodes = list(G.nodes())
    x_nodes = [pos[n][0] for n in nodes]
    y_nodes = [pos[n][1] for n in nodes]

    sizes = [10 + 22 * (deg[n] / max_deg) for n in nodes]
    colors = ["#ff6b6b" if G.nodes[n].get("side") == "source" else "#22c55e" for n in nodes]

    # Only label top hubs
    top_label_nodes = set(sorted(deg, key=lambda n: deg[n], reverse=True)[: int(label_top_k)])
    texts = [n if n in top_label_nodes else "" for n in nodes]

    node_hover = [f"{n}<br>degree={deg.get(n,0)}<br>side={G.nodes[n].get('side','?')}" for n in nodes]

    # Edge trace: SINGLE trace (much cleaner + faster). Hover is on nodes, not edges.
    edge_x, edge_y = [], []
    abs_scores = []
    for (u, v, data) in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
        abs_scores.append(float(data.get("abs_score", 0.0) or 0.0))

    max_abs = max(abs_scores) if abs_scores else 1.0
    # A nice moderate width; not per-edge width (keeps clean)
    edge_width = 1.2

    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            name="edges",
            hoverinfo="skip",
            line=dict(width=edge_width),
            opacity=0.18,
            showlegend=True,
        )
    )

    # Nodes on top
    fig.add_trace(
        go.Scatter(
            x=x_nodes,
            y=y_nodes,
            mode="markers+text",
            name="nodes",
            text=texts,
            textposition="top center",
            textfont=dict(size=10),
            marker=dict(size=sizes, color=colors, line=dict(width=0.7, color="white")),
            hoverinfo="text",
            hovertext=node_hover,
            showlegend=True,
        )
    )

    # Legend entries for groups
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", marker=dict(size=10, color="#ff6b6b"), name="sources"))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", marker=dict(size=10, color="#22c55e"), name="targets"))

    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(plot_bgcolor="white")
    return fig


def run_subprocess(cmd, cwd=None) -> tuple[int, str, str]:
    p = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    out, err = p.communicate()
    return p.returncode, out, err


def choose_network_path(uploaded_network_path: str | None, selected_default: str | None) -> Path:
    if uploaded_network_path:
        return Path(uploaded_network_path)
    if selected_default:
        return NETWORKS_DIR / selected_default
    nets = list_default_networks()
    if nets:
        return NETWORKS_DIR / nets[0]
    raise FileNotFoundError("No network provided and no default .rels exists in networks/.")


# -----------------------------
# Dash app
# -----------------------------
app = Dash(__name__)
app.title = "Inference Dash"

DEFAULT_NETS = list_default_networks()

app.layout = html.Div(
    style={"maxWidth": "980px", "margin": "28px auto 60px auto", "fontFamily": "system-ui, -apple-system, Segoe UI, Roboto, Arial"},
    children=[
        html.H1("Inference Dash", style={"marginBottom": "0.25rem"}),
        html.Div("Upload expression + (optional) network. Choose ONE engine. Run.", style={"color": "#5b6777", "marginBottom": "16px"}),
        html.Hr(),

        dcc.Store(id="store_expr_path"),
        dcc.Store(id="store_net_path"),
        dcc.Store(id="store_run_state"),
        dcc.Store(id="store_progress"),

        html.Div(
            style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "16px"},
            children=[
                html.Div(
                    style={"border": "1px solid #e6e8ee", "borderRadius": "12px", "padding": "14px", "background": "white"},
                    children=[
                        html.H4("1) Upload expression", style={"marginTop": 0}),
                        dcc.Upload(
                            id="upload_expr",
                            children=html.Div(["Drag & drop or select expression/signature file (csv/tsv)"]),
                            style={
                                "width": "100%",
                                "height": "64px",
                                "lineHeight": "64px",
                                "borderWidth": "2px",
                                "borderStyle": "dashed",
                                "borderRadius": "10px",
                                "textAlign": "center",
                                "borderColor": "#b7c2ff",
                                "background": "#fbfcff",
                                "cursor": "pointer",
                            },
                            multiple=False,
                        ),
                        html.Div(id="expr_selected", style={"marginTop": "10px", "color": "#304056"}),
                    ],
                ),
                html.Div(
                    style={"border": "1px solid #e6e8ee", "borderRadius": "12px", "padding": "14px", "background": "white"},
                    children=[
                        html.H4("2) Upload network (optional)", style={"marginTop": 0}),
                        html.Div("If you don’t upload, we use the default network.", style={"color": "#5b6777", "marginBottom": "8px"}),
                        dcc.Upload(
                            id="upload_net",
                            children=html.Div(["Drag & drop or select network (.rels)"]),
                            style={
                                "width": "100%",
                                "height": "64px",
                                "lineHeight": "64px",
                                "borderWidth": "2px",
                                "borderStyle": "dashed",
                                "borderRadius": "10px",
                                "textAlign": "center",
                                "borderColor": "#b7c2ff",
                                "background": "#fbfcff",
                                "cursor": "pointer",
                            },
                            multiple=False,
                        ),
                        html.Div(id="net_selected", style={"marginTop": "10px", "color": "#304056"}),

                        html.Div(style={"height": "10px"}),
                        html.Div("Or choose default network:", style={"fontWeight": 600}),
                        dcc.Dropdown(
                            id="default_network_dd",
                            options=[{"label": n, "value": n} for n in DEFAULT_NETS],
                            value=DEFAULT_NETS[0] if DEFAULT_NETS else None,
                            clearable=False,
                        ),
                        html.Div(id="default_net_path", style={"fontSize": "12px", "color": "#6a778a", "marginTop": "6px"}),
                    ],
                ),
            ],
        ),

        html.Div(style={"height": "14px"}),

        html.Div(
            style={"border": "1px solid #e6e8ee", "borderRadius": "12px", "padding": "14px", "background": "white"},
            children=[
                html.H4("3) Choose engine (only one)", style={"marginTop": 0}),
                dcc.RadioItems(
                    id="engine_radio",
                    options=[{"label": " CIE", "value": "cie"}, {"label": " ORNOR", "value": "ornor"}],
                    value="cie",
                    inline=True,
                    style={"fontWeight": 600},
                ),
                html.Div(style={"height": "10px"}),

                html.Div(
                    style={"display": "flex", "gap": "10px", "alignItems": "center", "flexWrap": "wrap"},
                    children=[
                        html.Div("ORNOR top_edges", style={"fontWeight": 600}),
                        dcc.Input(id="ornor_top_edges", type="number", value=5000, min=1, step=1, style={"width": "120px"}),
                        html.Div("(Only used for ORNOR)", style={"color": "#5b6777"}),
                    ],
                ),

                html.Div(style={"height": "12px"}),

                html.Button(
                    "Run",
                    id="run_btn",
                    n_clicks=0,
                    style={"background": "#0b1320", "color": "white", "border": "none", "padding": "9px 18px", "borderRadius": "10px", "cursor": "pointer", "fontWeight": 700},
                ),
                html.Div(id="network_hint", style={"marginTop": "10px", "color": "#5b6777"}),
            ],
        ),

        html.Div(style={"height": "16px"}),

        html.Div(
            style={"border": "1px solid #e6e8ee", "borderRadius": "12px", "padding": "14px", "background": "white"},
            children=[
                html.H4("Run status", style={"marginTop": 0}),
                html.Progress(id="progress_bar", value=0, max=100, style={"width": "100%", "height": "14px"}),
                html.Div(id="progress_text", style={"marginTop": "6px", "color": "#5b6777"}),

                html.Div(
                    id="status_box",
                    style={
                        "marginTop": "12px",
                        "background": "#0b1320",
                        "color": "#eaf0ff",
                        "borderRadius": "12px",
                        "padding": "12px",
                        "whiteSpace": "pre-wrap",
                        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace",
                        "fontSize": "12px",
                        "maxHeight": "360px",
                        "overflowY": "auto",
                    },
                    children="Waiting…",
                ),
            ],
        ),

        html.Div(style={"height": "16px"}),

        html.Div(
            style={"border": "1px solid #e6e8ee", "borderRadius": "12px", "padding": "14px", "background": "white"},
            children=[
                html.H4("Results", style={"marginTop": 0}),
                html.Div("Results will appear here after the run finishes.", style={"color": "#5b6777"}),

                html.Div(style={"height": "10px"}),

                html.Div(
                    style={"display": "flex", "gap": "10px", "flexWrap": "wrap"},
                    children=[
                        html.Button("Download edges CSV", id="btn_dl_edges"),
                        html.Button("Download TFs CSV", id="btn_dl_tfs"),
                        html.Button("Download normalized signature TSV", id="btn_dl_sig"),
                        dcc.Download(id="dl_edges"),
                        dcc.Download(id="dl_tfs"),
                        dcc.Download(id="dl_sig"),
                    ],
                ),

                html.Div(style={"height": "10px"}),
                html.Div("Top rows (preview):", style={"fontWeight": 700, "marginBottom": "8px"}),

                dash_table.DataTable(
                    id="table_preview",
                    columns=[],
                    data=[],
                    page_size=10,
                    style_table={"overflowX": "auto"},
                    style_header={"fontWeight": 700, "backgroundColor": "#f3f5fb"},
                    style_cell={"padding": "8px", "fontFamily": "system-ui"},
                ),

                html.Div(style={"height": "16px"}),

                dcc.Graph(id="fig_score_hist", figure=go.Figure(), config={"displaylogo": False}),
                dcc.Graph(id="fig_network", figure=go.Figure(), config={"displaylogo": False}),
            ],
        ),
    ],
)

# -----------------------------
# Callbacks: file selection
# -----------------------------
@app.callback(
    Output("store_expr_path", "data"),
    Output("expr_selected", "children"),
    Input("upload_expr", "contents"),
    State("upload_expr", "filename"),
    prevent_initial_call=True,
)
def on_upload_expr(contents, filename):
    if not contents:
        return no_update, no_update
    p = save_upload(contents, filename)
    return str(p), f"Selected: {Path(filename).name}  (saved as {p.name})"


@app.callback(
    Output("store_net_path", "data"),
    Output("net_selected", "children"),
    Input("upload_net", "contents"),
    State("upload_net", "filename"),
    prevent_initial_call=True,
)
def on_upload_net(contents, filename):
    if not contents:
        return no_update, no_update
    p = save_upload(contents, filename)
    return str(p), f"Uploaded network: {Path(filename).name}  (saved as {p.name})"


@app.callback(
    Output("default_net_path", "children"),
    Output("network_hint", "children"),
    Input("default_network_dd", "value"),
)
def show_default_path(netname):
    if not netname:
        return "Default network path: (none found)", "Pick a default network (dropdown) or upload a .rels file."
    p = NETWORKS_DIR / netname
    return f"Default network path: {p}", "Network: using selected default .rels."


# -----------------------------
# Run callback
# -----------------------------
@app.callback(
    Output("store_run_state", "data"),
    Output("store_progress", "data"),
    Output("status_box", "children"),
    Input("run_btn", "n_clicks"),
    State("store_expr_path", "data"),
    State("store_net_path", "data"),
    State("default_network_dd", "value"),
    State("engine_radio", "value"),
    State("ornor_top_edges", "value"),
    prevent_initial_call=True,
)
def run_engine(n, expr_path, net_path, default_net, engine, ornor_top_edges):
    if not n:
        return no_update, no_update, no_update

    if not expr_path:
        return {}, 0, "ERROR: please upload an expression/signature file first."

    try:
        expr_p = Path(expr_path)
        if not expr_p.exists():
            return {}, 0, f"ERROR: signature file not found: {expr_p}"
        network_p = choose_network_path(net_path, default_net)
        if not network_p.exists():
            return {}, 0, f"ERROR: network not found: {network_p}"
    except Exception as e:
        return {}, 0, f"ERROR selecting network: {e}"

    edges_out = RESULTS_DIR / f"{engine}_edges.csv"
    tfs_out = RESULTS_DIR / f"{engine}_edges.csv_tfs.csv"
    sig_out = RESULTS_DIR / f"{engine}_normalized_signature.tsv"

    status_lines = []
    status_lines.append(f"Engine: {engine.upper()}")
    status_lines.append(f"Signature: {expr_p}")
    status_lines.append(f"Network:   {network_p}")
    status_lines.append(f"Output:    {edges_out}")
    status_lines.append("")
    status_lines.append("Running:")

    progress = 10

    if engine == "cie":
        cmd = ["Rscript", str(CIE_R_SCRIPT), str(expr_p), str(network_p), str(edges_out)]
        status_lines.append(f"Rscript {' '.join(cmd[1:])}")
        code, out, err = run_subprocess(cmd, cwd=APP_ROOT)
    else:
        topn = int(ornor_top_edges or 5000)
        cmd = ["python3", str(ORNOR_PY_SCRIPT), str(expr_p), str(network_p), str(edges_out), str(topn)]
        status_lines.append(" ".join(cmd))
        code, out, err = run_subprocess(cmd, cwd=APP_ROOT)

    progress = 100

    status_lines.append("")
    status_lines.append("STDOUT:")
    status_lines.append(out.strip() if out else "(empty)")
    status_lines.append("")
    status_lines.append("STDERR:")
    status_lines.append(err.strip() if err else "(empty)")
    status_lines.append("")
    status_lines.append(f"Exit code: {code}")

    state = {
        "engine": engine,
        "expr_path": str(expr_p),
        "net_path": str(network_p),
        "edges_out": str(edges_out),
        "tfs_out": str(tfs_out),
        "sig_out": str(sig_out),
        "exit_code": int(code),
        "stdout": out,
        "stderr": err,
        "ts": time.time(),
    }

    return state, progress, "\n".join(status_lines)


# -----------------------------
# Progress UI
# -----------------------------
@app.callback(
    Output("progress_bar", "value"),
    Output("progress_text", "children"),
    Input("store_progress", "data"),
)
def render_progress(p):
    if p is None:
        return 0, ""
    try:
        p = int(p)
    except Exception:
        p = 0
    return p, f"{p}%"


# -----------------------------
# Render results (table + plots)
# -----------------------------
@app.callback(
    Output("table_preview", "columns"),
    Output("table_preview", "data"),
    Output("fig_score_hist", "figure"),
    Output("fig_network", "figure"),
    Input("store_run_state", "data"),
)
def render_results(state):
    if not state or "edges_out" not in state:
        empty_fig = go.Figure()
        return [], [], empty_fig, empty_fig

    edges_path = Path(state["edges_out"])
    if not edges_path.exists() or state.get("exit_code", 1) != 0:
        empty_fig = go.Figure()
        return [], [], empty_fig, empty_fig

    df = read_edges_csv(edges_path)

    mapping = load_entrez_symbol_map()
    if mapping and "source" in df.columns and "target" in df.columns:
        df = df.copy()
        df["source_name"] = maybe_add_symbols(df, "source", mapping)
        df["target_name"] = maybe_add_symbols(df, "target", mapping)
        front = ["source_name", "target_name"]
        keep = [c for c in df.columns if c not in front]
        df = df[front + keep]

    preview = df.head(50).copy()
    columns = [{"name": c, "id": c} for c in preview.columns]
    data = preview.to_dict("records")

    fig_hist = build_score_hist(df)

    # ✅ readable network defaults
    fig_net = build_clean_network(
        df,
        mapping=mapping,
        top_edges=160,        # fewer edges
        max_nodes=45,         # fewer nodes
        label_top_k=10,       # only label top hubs
        min_abs_score_quantile=0.70,  # drop weak edges automatically
    )

    return columns, data, fig_hist, fig_net


# -----------------------------
# Downloads
# -----------------------------
@app.callback(
    Output("dl_edges", "data"),
    Input("btn_dl_edges", "n_clicks"),
    State("store_run_state", "data"),
    prevent_initial_call=True,
)
def download_edges(n, state):
    if not n or not state:
        return no_update
    p = Path(state.get("edges_out", ""))
    if not p.exists():
        return no_update
    return dcc.send_file(str(p))


@app.callback(
    Output("dl_tfs", "data"),
    Input("btn_dl_tfs", "n_clicks"),
    State("store_run_state", "data"),
    prevent_initial_call=True,
)
def download_tfs(n, state):
    if not n or not state:
        return no_update
    p = Path(state.get("tfs_out", ""))
    if not p.exists():
        return no_update
    return dcc.send_file(str(p))


@app.callback(
    Output("dl_sig", "data"),
    Input("btn_dl_sig", "n_clicks"),
    State("store_run_state", "data"),
    prevent_initial_call=True,
)
def download_sig(n, state):
    if not n or not state:
        return no_update
    p = Path(state.get("sig_out", ""))
    if not p.exists():
        return no_update
    return dcc.send_file(str(p))


# -----------------------------
# Main
# -----------------------------
if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8050, debug=False)
