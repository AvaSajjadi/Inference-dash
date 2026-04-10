# app.py
# FINAL (Label mapping removed completely)
# ---------------------------------------
# - Removes Entrez→Symbol mapping upload + all related code/stores
# - Everything displays Entrez IDs (numbers) only
# - Keeps: Three/Five/Other network selection + "Other" network upload
# - Keeps: edge-hover highlighting, sparse network, ranking, plots, table

from __future__ import annotations

import base64
import io
import traceback
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

import dash
from dash import Input, Output, State, dcc, html, dash_table


# ----------------------------
# Paths / config
# ----------------------------
BASE_DIR = Path(__file__).resolve().parent
UPLOAD_DIR = BASE_DIR / "uploads"
NETWORK_DIR = BASE_DIR / "networks"
CUSTOM_NETWORK_DIR = UPLOAD_DIR / "custom_networks"

UPLOAD_DIR.mkdir(exist_ok=True)
NETWORK_DIR.mkdir(exist_ok=True)
CUSTOM_NETWORK_DIR.mkdir(exist_ok=True)

NETWORK_FILES = {
    "three": NETWORK_DIR / "three_tissues.rels",
    "five": NETWORK_DIR / "five_tissues.rels",
    # "other" comes from upload
}

DEFAULT_TOP_N = 25
DEFAULT_LOGFC_THRESH = 1.0
DEFAULT_TOP_K_PER_REG = 5


# ----------------------------
# Utilities
# ----------------------------
def _b64_to_text(contents: str) -> str:
    _, b64 = contents.split(",", 1)
    raw = base64.b64decode(b64)
    return raw.decode("utf-8", errors="replace")


def _save_upload(contents: str, out_path: Path) -> None:
    _, b64 = contents.split(",", 1)
    out_path.write_bytes(base64.b64decode(b64))


def _read_table_autodelim(text: str) -> pd.DataFrame:
    return pd.read_csv(io.StringIO(text), sep=None, engine="python")


def _normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def parse_signature(contents: str, filename: str) -> Tuple[pd.DataFrame, Path]:
    """
    Required columns (synonyms supported):
      - entrez (or gene_id, entrezid, etc.)
      - pval (or pvalue)
      - fc (or logfc)
    """
    text = _b64_to_text(contents)
    df = _read_table_autodelim(text)
    df = _normalize_cols(df)

    rename = {}
    for c in df.columns:
        if c in {"entrez", "entrezid", "geneid", "gene_id", "gene"}:
            rename[c] = "entrez"
        elif c in {"pval", "pvalue", "p_value"}:
            rename[c] = "pval"
        elif c in {"fc", "logfc", "log_fc", "foldchange"}:
            rename[c] = "fc"
    df = df.rename(columns=rename)

    required = {"entrez", "pval", "fc"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Signature missing columns {sorted(missing)}. Found: {list(df.columns)}")

    df["entrez"] = pd.to_numeric(df["entrez"], errors="coerce")
    df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    df["fc"] = pd.to_numeric(df["fc"], errors="coerce")
    df = df.dropna(subset=["entrez", "pval", "fc"]).copy()
    df["entrez"] = df["entrez"].astype(int)

    # pval safety
    df.loc[df["pval"] <= 0, "pval"] = np.nextafter(0, 1)

    out_path = UPLOAD_DIR / filename
    df.to_csv(out_path, index=False)
    return df, out_path


def load_rels(path: Path) -> pd.DataFrame:
    """
    Reads *.rels like: uid srcuid trguid type pmids nls  (whitespace separated)
    Returns: source(int), target(int), type(str)
    """
    if not path.exists():
        raise FileNotFoundError(f"Network file not found: {path}")

    df = pd.read_csv(path, sep=r"\s+")
    df = _normalize_cols(df)

    required = {"srcuid", "trguid", "type"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} missing columns {sorted(missing)}. Found: {list(df.columns)}")

    out = pd.DataFrame(
        {
            "source": pd.to_numeric(df["srcuid"], errors="coerce"),
            "target": pd.to_numeric(df["trguid"], errors="coerce"),
            "type": df["type"].astype(str).str.lower().str.strip(),
        }
    ).dropna(subset=["source", "target"])
    out["source"] = out["source"].astype(int)
    out["target"] = out["target"].astype(int)
    out.loc[~out["type"].isin(["increase", "decrease"]), "type"] = "unknown"
    return out


def baseline_rank_regulators(sig: pd.DataFrame, rels: pd.DataFrame, top_n: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    score(r) = sum_{targets connected to r and present in sig} |fc(target)|
    """
    sig_targets = set(sig["entrez"].tolist())
    edges_hit = rels[rels["target"].isin(sig_targets)].copy()

    if edges_hit.empty:
        regs = pd.DataFrame(columns=["regulator", "score", "n_targets", "mean_abs_fc"])
        return regs, edges_hit

    sig_map = sig.set_index("entrez")[["fc", "pval"]]
    edges_hit = edges_hit.join(sig_map, on="target")
    edges_hit["abs_fc"] = edges_hit["fc"].abs()

    regs = (
        edges_hit.groupby("source", as_index=False)
        .agg(
            n_targets=("target", "nunique"),
            score=("abs_fc", "sum"),
            mean_abs_fc=("abs_fc", "mean"),
        )
        .sort_values(["score", "n_targets"], ascending=[False, False])
        .head(max(1, int(top_n)))
        .rename(columns={"source": "regulator"})
    )

    return regs, edges_hit


def _spread(n: int) -> List[float]:
    if n <= 1:
        return [0.5]
    return list(np.linspace(0, 1, n))


def build_sparse_subnetwork(
    edges_hit: pd.DataFrame,
    regs_df: pd.DataFrame,
    sig: pd.DataFrame,
    logfc_thresh: float,
    top_k_per_reg: int,
) -> pd.DataFrame:
    if regs_df.empty or edges_hit.empty:
        return pd.DataFrame()

    sig_f = sig[sig["fc"].abs() >= float(logfc_thresh)].copy()
    keep_targets = set(sig_f["entrez"].astype(int).tolist())
    if not keep_targets:
        return pd.DataFrame()

    k = int(top_k_per_reg) if top_k_per_reg is not None else DEFAULT_TOP_K_PER_REG
    k = max(1, min(50, k))

    sub_list = []
    for r in regs_df["regulator"].astype(int).tolist():
        r_edges = edges_hit[(edges_hit["source"] == r) & (edges_hit["target"].isin(keep_targets))].copy()
        if r_edges.empty:
            continue
        r_edges["abs_fc"] = r_edges["fc"].abs()
        r_edges = r_edges.sort_values("abs_fc", ascending=False).head(k)
        sub_list.append(r_edges)

    if not sub_list:
        return pd.DataFrame()

    return pd.concat(sub_list, ignore_index=True)


def make_network_figure_from_sub(
    sub: pd.DataFrame,
    sig: pd.DataFrame,
    logfc_thresh: float,
    top_k_per_reg: int,
    highlight_edge_id: Optional[int],
) -> go.Figure:
    fig = go.Figure()

    if sub.empty:
        fig.update_layout(
            title="Network view (no edges after filtering / no matches)",
            margin=dict(l=20, r=20, t=60, b=20),
        )
        return fig

    regs = sorted(sub["source"].unique().tolist())
    genes = sorted(sub["target"].unique().tolist())

    reg_pos = {r: (0.0, _spread(len(regs))[i]) for i, r in enumerate(regs)}
    gene_pos = {g: (1.0, _spread(len(genes))[i]) for i, g in enumerate(genes)}

    fc_map = sig.set_index("entrez")["fc"].to_dict()

    base_color = "#b0b0b0"
    base_width = 1.4
    base_opacity = 0.55

    hi_color = "#111111"
    hi_width = 5.0
    hi_opacity = 1.0

    # edges: one trace per edge so hover can highlight exactly that edge
    for idx, row in sub.reset_index(drop=True).iterrows():
        r = int(row["source"])
        g = int(row["target"])
        x0, y0 = reg_pos[r]
        x1, y1 = gene_pos[g]

        is_hi = highlight_edge_id is not None and idx == int(highlight_edge_id)

        fig.add_trace(
            go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(width=hi_width if is_hi else base_width, color=hi_color if is_hi else base_color),
                opacity=hi_opacity if is_hi else base_opacity,
                hoverinfo="text",
                hovertext=(
                    f"<b>{r}</b> → <b>{g}</b><br>"
                    f"type: {row.get('type','')}<br>"
                    f"logFC(target): {float(fc_map.get(g, np.nan)):.3f}"
                ),
                customdata=[idx],
                showlegend=False,
            )
        )

    # regulator nodes
    fig.add_trace(
        go.Scatter(
            x=[reg_pos[r][0] for r in regs],
            y=[reg_pos[r][1] for r in regs],
            mode="markers+text",
            text=[str(r) for r in regs],
            textposition="middle left",
            marker=dict(size=14, symbol="circle", color="#7b2cbf", line=dict(width=1)),
            hovertext=[f"Regulator: <b>{r}</b>" for r in regs],
            hoverinfo="text",
            name="Regulators",
        )
    )

    # gene nodes (colored by logFC)
    gene_fc = [float(fc_map.get(g, 0.0)) for g in genes]
    fig.add_trace(
        go.Scatter(
            x=[gene_pos[g][0] for g in genes],
            y=[gene_pos[g][1] for g in genes],
            mode="markers+text",
            text=[str(g) for g in genes],
            textposition="middle right",
            marker=dict(
                size=10,
                color=gene_fc,
                colorscale="RdBu",
                reversescale=True,
                colorbar=dict(title="logFC"),
                line=dict(width=0.7),
            ),
            hovertext=[f"Gene: <b>{g}</b><br>logFC: {float(fc_map.get(g, np.nan)):.3f}" for g in genes],
            hoverinfo="text",
            name="Genes",
        )
    )

    fig.update_layout(
        title=f"Network view (bipartite, sparse): hover an EDGE to highlight it • top-{top_k_per_reg} targets/reg, |logFC| ≥ {logfc_thresh}",
        margin=dict(l=20, r=20, t=70, b=20),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.15, 1.15]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.05, 1.05]),
        hovermode="closest",
    )
    return fig


# ----------------------------
# Dash UI
# ----------------------------
app = dash.Dash(__name__)
app.title = "Inference Dash"

app.layout = html.Div(
    style={"maxWidth": "1200px", "margin": "0 auto", "padding": "18px", "fontFamily": "Arial"},
    children=[
        html.H2("Inference Dash"),
        html.Div("Signature + Prior network → regulator ranking + sparse bipartite network view"),
        html.Hr(),

        html.Div(
            children=[
                html.Div(html.B("1) Upload signature (entrez, pval, fc)")),
                dcc.Upload(
                    id="upload_sig",
                    children=html.Div(["Drag & drop or ", html.B("Select signature file")]),
                    style={
                        "width": "100%", "height": "70px", "lineHeight": "70px",
                        "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "6px",
                        "textAlign": "center",
                    },
                    multiple=False,
                ),
                html.Div(id="sig_status", style={"marginTop": "8px"}),
            ]
        ),

        html.Hr(),

        html.Div(
            style={"display": "flex", "gap": "18px", "alignItems": "center", "flexWrap": "wrap"},
            children=[
                html.Div(
                    children=[
                        html.Div(html.B("Network")),
                        dcc.RadioItems(
                            id="network_choice",
                            options=[
                                {"label": "Three tissues", "value": "three"},
                                {"label": "Five tissues", "value": "five"},
                                {"label": "Other", "value": "other"},
                            ],
                            value="three",
                            inline=True,
                        ),
                    ]
                ),
                html.Div(
                    children=[
                        html.Div(html.B("Top N regulators (ranking)")),
                        dcc.Input(
                            id="top_n",
                            type="number",
                            min=5,
                            max=200,
                            step=1,
                            value=DEFAULT_TOP_N,
                            style={"width": "110px"},
                        ),
                    ]
                ),
                html.Div(
                    children=[
                        html.Div(html.B("Network filter: |logFC| ≥")),
                        dcc.Input(
                            id="logfc_thresh",
                            type="number",
                            min=0,
                            max=10,
                            step=0.1,
                            value=DEFAULT_LOGFC_THRESH,
                            style={"width": "110px"},
                        ),
                    ]
                ),
                html.Div(
                    children=[
                        html.Div(html.B("Network sparsity: top-K targets/reg")),
                        dcc.Input(
                            id="top_k_per_reg",
                            type="number",
                            min=1,
                            max=50,
                            step=1,
                            value=DEFAULT_TOP_K_PER_REG,
                            style={"width": "110px"},
                        ),
                    ]
                ),
                html.Button("Run", id="run_btn", n_clicks=0, style={"height": "36px"}),
            ],
        ),

        # Network upload (only when "Other" selected)
        html.Div(
            id="other_network_box",
            style={"marginTop": "12px", "display": "none"},
            children=[
                html.Div(html.B("Upload your network (.rels)")),
                dcc.Upload(
                    id="upload_net",
                    children=html.Div(["Drag & drop or ", html.B("Select network file")]),
                    style={
                        "width": "100%",
                        "height": "70px",
                        "lineHeight": "70px",
                        "borderWidth": "1px",
                        "borderStyle": "dashed",
                        "borderRadius": "6px",
                        "textAlign": "center",
                        "marginTop": "6px",
                    },
                    multiple=False,
                ),
                html.Div(id="net_status", style={"marginTop": "8px"}),
            ],
        ),

        html.Div(id="run_status", style={"marginTop": "12px", "whiteSpace": "pre-wrap"}),

        html.Hr(),

        html.H3("Summary plot"),
        dcc.Loading(type="default", children=dcc.Graph(id="summary_plot", figure=go.Figure())),

        html.H3("Top regulators (table)"),
        dash_table.DataTable(
            id="reg_table",
            columns=[
                {"name": "Regulator (Entrez)", "id": "regulator"},
                {"name": "n_targets", "id": "n_targets"},
                {"name": "score", "id": "score"},
                {"name": "mean_abs_fc", "id": "mean_abs_fc"},
            ],
            data=[],
            sort_action="native",
            filter_action="native",
            page_size=12,
            style_table={"overflowX": "auto"},
            style_cell={"padding": "6px", "fontSize": "14px"},
        ),

        html.H3("Network view (hover an EDGE to highlight it)"),
        dcc.Loading(
            type="default",
            children=dcc.Graph(
                id="network_plot",
                figure=go.Figure(),
                style={"height": "720px"},
                clear_on_unhover=True,
            ),
        ),

        dcc.Store(id="sig_store"),
        dcc.Store(id="custom_net_store"),
        dcc.Store(id="network_store"),
    ],
)


# ----------------------------
# UI: show/hide other-network upload
# ----------------------------
@app.callback(
    Output("other_network_box", "style"),
    Input("network_choice", "value"),
)
def toggle_other_network_box(choice: str):
    if choice == "other":
        return {"marginTop": "12px", "display": "block"}
    return {"marginTop": "12px", "display": "none"}


# ----------------------------
# Upload callbacks
# ----------------------------
@app.callback(
    Output("sig_status", "children"),
    Output("sig_store", "data"),
    Input("upload_sig", "contents"),
    State("upload_sig", "filename"),
    prevent_initial_call=True,
)
def on_upload_sig(contents, filename):
    try:
        if not contents or not filename:
            return "No signature uploaded.", None

        _save_upload(contents, UPLOAD_DIR / filename)
        sig, _ = parse_signature(contents, filename)

        msg = html.Div(
            [
                html.Div("✅ Signature uploaded", style={"color": "green", "fontWeight": "bold"}),
                html.Div(f"Saved: uploads/{filename}"),
                html.Div(f"Rows: {sig.shape[0]}"),
                html.Div(f"Columns: {list(sig.columns)}"),
            ]
        )
        return msg, sig.to_json(date_format="iso", orient="split")

    except Exception:
        return html.Pre("❌ Signature upload failed:\n\n" + traceback.format_exc()), None


@app.callback(
    Output("net_status", "children"),
    Output("custom_net_store", "data"),
    Input("upload_net", "contents"),
    State("upload_net", "filename"),
    prevent_initial_call=True,
)
def on_upload_net(contents, filename):
    try:
        if not contents or not filename:
            return "No network uploaded.", None

        out_path = CUSTOM_NETWORK_DIR / filename
        _save_upload(contents, out_path)

        # validate format
        _ = load_rels(out_path)

        msg = html.Div(
            [
                html.Div("✅ Network uploaded", style={"color": "green", "fontWeight": "bold"}),
                html.Div(f"Saved: uploads/custom_networks/{filename}"),
            ]
        )
        return msg, str(out_path)

    except Exception:
        return html.Pre("❌ Network upload failed:\n\n" + traceback.format_exc()), None


# ----------------------------
# Run callback
# ----------------------------
@app.callback(
    Output("run_status", "children"),
    Output("summary_plot", "figure"),
    Output("reg_table", "data"),
    Output("network_store", "data"),
    Input("run_btn", "n_clicks"),
    State("sig_store", "data"),
    State("network_choice", "value"),
    State("custom_net_store", "data"),
    State("top_n", "value"),
    State("logfc_thresh", "value"),
    State("top_k_per_reg", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, sig_json, network_choice, custom_net_path, top_n, logfc_thresh, top_k_per_reg):
    empty = go.Figure()
    try:
        if not sig_json:
            return "❌ Upload a signature file first.", empty, [], None

        sig = pd.read_json(sig_json, orient="split")

        if network_choice == "other":
            if not custom_net_path:
                return "❌ You selected 'Other' but did not upload a network .rels file.", empty, [], None
            net_path = Path(custom_net_path)
        else:
            net_path = NETWORK_FILES.get(network_choice, NETWORK_FILES["three"])

        rels = load_rels(net_path)

        top_n = int(top_n) if top_n else DEFAULT_TOP_N
        top_n = max(5, min(200, top_n))

        regs_df, edges_hit = baseline_rank_regulators(sig, rels, top_n=top_n)

        # Summary plot
        if regs_df.empty:
            summary_fig = go.Figure()
            summary_fig.update_layout(title="Summary plot (no matching edges for this signature + network)")
        else:
            plot_df = regs_df.copy()
            plot_df["label"] = plot_df["regulator"].astype(str)
            summary_fig = px.bar(
                plot_df,
                x="label",
                y="score",
                hover_data=["n_targets", "mean_abs_fc", "regulator"],
                title=f"Top {len(plot_df)} regulators (score = sum(|logFC|) over connected targets)",
                labels={"label": "Regulator (Entrez)", "score": "Score"},
            )
            summary_fig.update_layout(xaxis_tickangle=45)

        # Table data
        table_data = []
        if not regs_df.empty:
            out = regs_df.copy()
            out["score"] = out["score"].astype(float).round(4)
            out["mean_abs_fc"] = out["mean_abs_fc"].astype(float).round(4)
            table_data = out[["regulator", "n_targets", "score", "mean_abs_fc"]].to_dict("records")

        logfc_thresh = float(logfc_thresh) if logfc_thresh is not None else DEFAULT_LOGFC_THRESH
        top_k_per_reg = int(top_k_per_reg) if top_k_per_reg is not None else DEFAULT_TOP_K_PER_REG
        top_k_per_reg = max(1, min(50, top_k_per_reg))

        sub = build_sparse_subnetwork(edges_hit, regs_df, sig, logfc_thresh, top_k_per_reg)

        network_store = {
            "sub": sub.to_dict("records"),
            "sig": sig.to_dict("records"),
            "logfc_thresh": logfc_thresh,
            "top_k_per_reg": top_k_per_reg,
        }

        status_lines = [
            "✅ Run complete.",
            f"Network: {net_path.name}",
            f"Network edges loaded: {rels.shape[0]}",
            f"Signature genes: {sig.shape[0]}",
            f"Edges matching signature targets (ranking): {edges_hit.shape[0]}",
            f"Top regulators shown: {0 if regs_df.empty else regs_df.shape[0]}",
            f"Network visualization filter: |logFC| ≥ {logfc_thresh}",
            f"Network visualization sparsity: top-K targets/regulator = {top_k_per_reg}",
        ]

        return "\n".join(status_lines), summary_fig, table_data, network_store

    except Exception:
        return "❌ Run failed:\n\n" + traceback.format_exc(), empty, [], None


# ----------------------------
# Network hover callback: highlight hovered EDGE
# ----------------------------
@app.callback(
    Output("network_plot", "figure"),
    Input("network_store", "data"),
    Input("network_plot", "hoverData"),
)
def update_network_figure(network_store: Optional[dict], hoverData: Optional[dict]):
    if not network_store:
        fig = go.Figure()
        fig.update_layout(title="Upload signature + click Run to build network")
        return fig

    try:
        sub = pd.DataFrame(network_store["sub"])
        sig = pd.DataFrame(network_store["sig"])
        logfc_thresh = float(network_store.get("logfc_thresh", DEFAULT_LOGFC_THRESH))
        top_k_per_reg = int(network_store.get("top_k_per_reg", DEFAULT_TOP_K_PER_REG))

        highlight_edge_id = None
        if hoverData and "points" in hoverData and hoverData["points"]:
            pt = hoverData["points"][0]
            if "customdata" in pt and pt["customdata"] is not None:
                try:
                    highlight_edge_id = int(pt["customdata"])
                except Exception:
                    highlight_edge_id = None

        return make_network_figure_from_sub(
            sub=sub,
            sig=sig,
            logfc_thresh=logfc_thresh,
            top_k_per_reg=top_k_per_reg,
            highlight_edge_id=highlight_edge_id,
        )

    except Exception:
        fig = go.Figure()
        fig.update_layout(title="❌ Network figure failed")
        return fig


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8050, debug=True)
