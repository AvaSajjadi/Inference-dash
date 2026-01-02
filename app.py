from pathlib import Path
import base64
import subprocess
import pandas as pd

from dash import Dash, html, dcc, Input, Output, State


# =====================
# Paths
# =====================
BASE_DIR = Path(__file__).resolve().parent
UPLOADS_DIR = BASE_DIR / "uploads"
RESULTS_DIR = BASE_DIR / "results"

UPLOADS_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)


# =====================
# Dash app
# =====================
app = Dash(__name__)

app.layout = html.Div(
    style={"maxWidth": "900px", "margin": "40px"},
    children=[
        html.H1("Inference Dash"),
        html.P("Unified interface for CIE (statistical, R) and ORNOR (Bayesian) inference"),
        html.Hr(),

        html.H3("Inference method"),
        dcc.RadioItems(
            id="method",
            options=[
                {"label": "CIE (Statistical, R)", "value": "cie"},
                {"label": "ORNOR (Bayesian)", "value": "ornor"},
            ],
            value="cie",
        ),

        html.Br(),
        html.H3("Upload expression matrix (CSV)"),
        dcc.Upload(
            id="upload",
            children=html.Div(["Drag and Drop or Select File"]),
            style={
                "width": "100%",
                "height": "60px",
                "lineHeight": "60px",
                "borderWidth": "1px",
                "borderStyle": "dashed",
                "borderRadius": "5px",
                "textAlign": "center",
            },
        ),

        html.Br(),
        html.Button("Run inference", id="run-btn"),
        html.Br(),
        html.Br(),

        html.Div(id="status"),
        html.Div(id="preview"),
    ],
)


# =====================
# Helpers
# =====================
def render_success(df: pd.DataFrame, message: str):
    return html.Div(
        children=[
            html.P(f"✅ {message}", style={"color": "green"}),
            html.Hr(),
            html.H4("Output preview"),
            html.P(f"Rows: {df.shape[0]} | Columns: {df.shape[1]}"),
            html.Pre(df.head(10).to_string(index=False)),
        ]
    )


def render_error(msg: str):
    return html.Div(
        children=[
            html.P("❌ Error:", style={"color": "red", "fontWeight": "bold"}),
            html.Pre(msg),
        ]
    )


# =====================
# Callback
# =====================
@app.callback(
    Output("status", "children"),
    Output("preview", "children"),
    Input("run-btn", "n_clicks"),
    State("method", "value"),
    State("upload", "contents"),
    State("upload", "filename"),
)
def run_inference(n_clicks, method, contents, filename):
    if not n_clicks:
        return "", ""

    if contents is None:
        return render_error("No file uploaded."), ""

    try:
        # ---- Decode upload ----
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)

        upload_path = UPLOADS_DIR / filename
        with open(upload_path, "wb") as f:
            f.write(decoded)

        # =====================
        # CIE (R)
        # =====================
        if method == "cie":
            out_csv = RESULTS_DIR / "cie_output.csv"

            cmd = [
                "Rscript",
                str(BASE_DIR / "run_cie.R"),
                str(upload_path),
                str(BASE_DIR / "networks" / "human_network.csv"),
                str(out_csv),
            ]

            subprocess.check_call(cmd)

            df = pd.read_csv(out_csv)
            return render_success(df, "CIE inference completed"), ""

        # =====================
        # ORNOR (stub)
        # =====================
        elif method == "ornor":
            from backend.ornor import run_ornor

            out_csv = RESULTS_DIR / "ornor_output.csv"
            run_ornor(upload_path, out_csv)

            df = pd.read_csv(out_csv)
            return render_success(df, "ORNOR inference completed"), ""

        else:
            return render_error("Unknown inference method."), ""

    except subprocess.CalledProcessError as e:
        return render_error(str(e)), ""

    except Exception as e:
        return render_error(str(e)), ""


# =====================
# Main
# =====================
if __name__ == "__main__":
    app.run(debug=True)
