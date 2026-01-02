import base64
import subprocess
from pathlib import Path

from dash import Dash, html, dcc, Input, Output, State
import pandas as pd

# -------------------------
# Paths
# -------------------------
BASE_DIR = Path(__file__).parent
UPLOAD_DIR = BASE_DIR / "uploads"
RESULTS_DIR = BASE_DIR / "results"

UPLOAD_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

RUN_CIE_SCRIPT = BASE_DIR / "run_cie.R"

# -------------------------
# App
# -------------------------
app = Dash(__name__)
server = app.server

# -------------------------
# Layout
# -------------------------
app.layout = html.Div(
    style={"width": "70%", "margin": "40px auto", "fontFamily": "serif"},
    children=[
        html.H1("Inference Dash"),
        html.P(
            "Unified interface for CIE (statistical, R) and ORNOR (Bayesian) inference"
        ),
        html.Hr(),

        html.H3("Inference method"),
        dcc.RadioItems(
            id="method",
            options=[
                {"label": "CIE (Statistical, R)", "value": "cie"},
                {"label": "ORNOR (Bayesian)", "value": "ornor"},
            ],
            value="cie",
            inline=True,
        ),

        html.Br(),
        html.H3("Upload expression matrix (CSV)"),
        dcc.Upload(
            id="upload",
            children=html.Div(
                ["Drag and Drop or ", html.A("Select File")]
            ),
            style={
                "width": "100%",
                "height": "80px",
                "lineHeight": "80px",
                "borderWidth": "1px",
                "borderStyle": "dashed",
                "borderRadius": "5px",
                "textAlign": "center",
            },
            multiple=False,
        ),

        html.Br(),
        html.Button("Run inference", id="run", n_clicks=0),

        html.Br(),
        html.Br(),
        html.Div(id="status"),
        html.Hr(),
        html.Div(id="output-preview"),
    ],
)

# -------------------------
# Callback
# -------------------------
@app.callback(
    Output("status", "children"),
    Output("output-preview", "children"),
    Input("run", "n_clicks"),
    State("method", "value"),
    State("upload", "contents"),
    State("upload", "filename"),
    prevent_initial_call=True,
)
def run_inference(n_clicks, method, contents, filename):
    if contents is None:
        return html.Span("❌ No file uploaded", style={"color": "red"}), None

    try:
        # -------------------------
        # Decode uploaded file
        # -------------------------
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)

        expr_path = UPLOAD_DIR / filename
        with open(expr_path, "wb") as f:
            f.write(decoded)

        # -------------------------
        # Run CIE
        # -------------------------
        if method == "cie":
            network_csv = BASE_DIR / "networks" / "human_network.csv"
            out_csv = RESULTS_DIR / "cie_output.csv"

            cmd = [
                "Rscript",
                str(RUN_CIE_SCRIPT),
                str(expr_path),
                str(network_csv),
                str(out_csv),
            ]

            subprocess.run(cmd, check=True)

            df = pd.read_csv(out_csv)
            preview = df.head(20)

            return (
                html.Span("✅ CIE inference completed", style={"color": "green"}),
                html.Div(
                    [
                        html.H4("Output preview"),
                        dcc.Markdown(f"Rows: {df.shape[0]} | Columns: {df.shape[1]}"),
                        dcc.Markdown(preview.to_markdown(index=False)),
                    ]
                ),
            )

        # -------------------------
        # ORNOR placeholder
        # -------------------------
        else:
            return (
                html.Span(
                    "⚠️ ORNOR backend not wired yet",
                    style={"color": "orange"},
                ),
                None,
            )

    except subprocess.CalledProcessError as e:
        return (
            html.Div(
                [
                    html.Span("❌ Inference failed:", style={"color": "red"}),
                    html.Pre(str(e)),
                ]
            ),
            None,
        )

    except Exception as e:
        return (
            html.Div(
                [
                    html.Span("❌ Error:", style={"color": "red"}),
                    html.Pre(str(e)),
                ]
            ),
            None,
        )


# -------------------------
# Run
# -------------------------
if __name__ == "__main__":
    app.run(debug=True)
