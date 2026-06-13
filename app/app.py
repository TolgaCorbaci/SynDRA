"""
app.py
======
Phase 9: SynDRA Shiny web application.

Two tabs:
  Lookup     - resolve a drug name/ID to its canonical node + all synonyms + xrefs
  Enrichment - paste a drug list -> ORA across all 27 libraries + coverage

Runs against syndra_redistributable outputs.
Run: shiny run app/app.py  (from project root)
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_HERE = Path(__file__).parent
_ROOT = _HERE.parent
sys.path.insert(0, str(_ROOT / "build"))
sys.path.insert(0, str(_ROOT / "enrichment"))

import pandas as pd
from shiny import App, Inputs, Outputs, Session, render, reactive, ui

# ---------------------------------------------------------------------------
# Data loading (lazy; app still starts if outputs aren't built yet)
# ---------------------------------------------------------------------------

OUTPUT_DIR = _ROOT / "outputs"
DB_DIR = _ROOT / "synonyms" / "input" / "enrichment_databases"


def _try_load(parquet_name: str, csv_name: str | None = None) -> pd.DataFrame:
    p = OUTPUT_DIR / parquet_name
    if p.exists():
        return pd.read_parquet(p)
    if csv_name:
        c = OUTPUT_DIR / csv_name
        if c.exists():
            return pd.read_csv(c, dtype=str).fillna("")
    return pd.DataFrame()


def _load_all():
    compounds = _try_load("syndra_redistributable_compounds.parquet",
                          "syndra_redistributable_compounds.csv")
    xrefs = _try_load("syndra_redistributable_xrefs.parquet",
                      "syndra_redistributable_xrefs.csv")
    synonyms = _try_load("syndra_redistributable_synonyms.parquet",
                         "syndra_redistributable_synonyms.csv")
    return compounds, xrefs, synonyms


_compounds_df, _xrefs_df, _synonyms_df = _load_all()
_OUTPUTS_READY = not _compounds_df.empty


# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------

def _build_lookup_index(synonyms_df: pd.DataFrame) -> dict[str, str]:
    """synonym_norm -> syndra_id"""
    idx: dict[str, str] = {}
    for _, row in synonyms_df.iterrows():
        norm = str(row.get("synonym_norm", "")).strip()
        sid = str(row.get("syndra_id", "")).strip()
        if norm and sid and norm not in idx:
            idx[norm] = sid
    return idx


def _resolve(query: str, idx: dict[str, str]) -> str | None:
    from normalize import normalize_name
    norm = normalize_name(query.strip())
    return idx.get(norm)


def _get_compound(sid: str, compounds_df: pd.DataFrame) -> dict:
    rows = compounds_df[compounds_df["syndra_id"] == sid]
    return rows.iloc[0].to_dict() if not rows.empty else {}


def _get_xrefs(sid: str, xrefs_df: pd.DataFrame) -> pd.DataFrame:
    return xrefs_df[xrefs_df["syndra_id"] == sid][["id_type", "id_value", "source"]]


def _get_synonyms(sid: str, synonyms_df: pd.DataFrame) -> pd.DataFrame:
    return synonyms_df[synonyms_df["syndra_id"] == sid][
        ["synonym_raw", "synonym_type", "source"]
    ]


# ---------------------------------------------------------------------------
# Enrichment helpers
# ---------------------------------------------------------------------------

def _run_enrichment(drug_list: list[str]) -> tuple[pd.DataFrame, dict, list[str]]:
    from syndra_enrich import enrich
    return enrich(
        query_names=drug_list,
        synonyms_df=_synonyms_df,
        db_dir=str(DB_DIR),
        min_overlap=1,
        fdr_threshold=0.05,
    )


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

app_ui = ui.page_navbar(
    ui.nav_panel(
        "Lookup",
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Search"),
                ui.input_text("lookup_query", "Drug name or ID",
                              placeholder="e.g. aspirin, BRD-K12345678"),
                ui.input_action_button("lookup_btn", "Lookup", class_="btn-primary w-100"),
                ui.hr(),
                ui.p(ui.tags.small(
                    "Enter any synonym, INN, trade name, BRD ID, TTD ID, or PubChem CID."
                )),
            ),
            ui.card(
                ui.card_header("Result"),
                ui.output_ui("lookup_result"),
            ),
        ),
    ),
    ui.nav_panel(
        "Enrichment",
        ui.layout_sidebar(
            ui.sidebar(
                ui.h5("Drug Set Input"),
                ui.input_textarea(
                    "enrich_query",
                    "Drug names (one per line or comma-separated)",
                    placeholder="atorvastatin\nsimvastatin\npravastatin\nlovastatin",
                    rows=8,
                ),
                ui.input_action_button("enrich_btn", "Run Enrichment",
                                       class_="btn-primary w-100"),
                ui.hr(),
                ui.input_numeric("fdr_thresh", "FDR threshold", value=0.05,
                                 min=0.001, max=1.0, step=0.01),
                ui.input_numeric("min_overlap", "Min overlap", value=1, min=1, max=10),
            ),
            ui.card(
                ui.card_header("Enrichment Results"),
                ui.output_ui("enrich_status"),
                ui.output_data_frame("enrich_table"),
            ),
            ui.card(
                ui.card_header("Library Coverage"),
                ui.output_data_frame("coverage_table"),
            ),
        ),
    ),
    title="SynDRA",
    id="main_nav",
    bg="#1a1a2e",
    inverse=True,
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def server(input: Inputs, output: Outputs, session: Session):

    if _OUTPUTS_READY:
        _idx = _build_lookup_index(_synonyms_df)
    else:
        _idx = {}

    # ------------------------------------------------------------------
    # Lookup tab
    # ------------------------------------------------------------------

    @output
    @render.ui
    @reactive.event(input.lookup_btn)
    def lookup_result():
        if not _OUTPUTS_READY:
            return ui.div(
                ui.tags.div(
                    "SynDRA outputs not built yet. Run ",
                    ui.tags.code("python build/build_all.py"),
                    " from the project root first.",
                    class_="alert alert-warning",
                )
            )

        query = input.lookup_query().strip()
        if not query:
            return ui.p("Enter a drug name or ID above.")

        sid = _resolve(query, _idx)

        # Try direct syndra_id or xref lookup
        if sid is None and _xrefs_df is not None and not _xrefs_df.empty:
            q_upper = query.strip().upper()
            hits = _xrefs_df[_xrefs_df["id_value"].str.upper() == q_upper]
            if not hits.empty:
                sid = hits.iloc[0]["syndra_id"]

        if sid is None:
            return ui.div(
                ui.tags.div(f"No match found for: {query!r}", class_="alert alert-info")
            )

        compound = _get_compound(sid, _compounds_df)
        xrefs = _get_xrefs(sid, _xrefs_df)
        syns = _get_synonyms(sid, _synonyms_df)

        return ui.div(
            ui.h4(compound.get("preferred_name", sid)),
            ui.tags.table(
                ui.tags.tr(ui.tags.td(ui.tags.b("SynDRA ID")), ui.tags.td(sid)),
                ui.tags.tr(ui.tags.td(ui.tags.b("InChIKey")),
                           ui.tags.td(compound.get("inchikey", "—"))),
                ui.tags.tr(ui.tags.td(ui.tags.b("Has structure")),
                           ui.tags.td("Yes" if compound.get("has_structure") else "No")),
                ui.tags.tr(ui.tags.td(ui.tags.b("SMILES")),
                           ui.tags.td(ui.tags.small(compound.get("canonical_smiles", "—")))),
                class_="table table-sm table-bordered mb-3",
            ),
            ui.h6("External identifiers"),
            ui.output_data_frame("xref_table_inner") if False else
            _df_to_table(xrefs, ["id_type", "id_value", "source"]),
            ui.h6("Synonyms"),
            _df_to_table(syns.head(80), ["synonym_raw", "synonym_type", "source"]),
            ui.p(ui.tags.small(f"Showing up to 80 of {len(syns)} synonyms."))
            if len(syns) > 80 else ui.span(),
        )

    # ------------------------------------------------------------------
    # Enrichment tab
    # ------------------------------------------------------------------

    _enrich_data = reactive.Value(None)

    @reactive.effect
    @reactive.event(input.enrich_btn)
    def _do_enrich():
        if not _OUTPUTS_READY:
            _enrich_data.set(("error", "Run build/build_all.py first."))
            return

        raw = input.enrich_query().strip()
        if not raw:
            _enrich_data.set(("error", "Paste drug names above."))
            return

        # Split on newlines and commas
        import re
        drugs = [d.strip() for d in re.split(r"[\n,]+", raw) if d.strip()]
        if not drugs:
            _enrich_data.set(("error", "No drug names detected."))
            return

        try:
            coverage_df, results, unresolved = _run_enrichment(drugs)
            # Combine all library results
            rows = []
            for lib, df in results.items():
                if not df.empty:
                    df = df.copy()
                    df.insert(0, "library", lib)
                    rows.append(df)
            combined = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
            _enrich_data.set(("ok", combined, coverage_df, unresolved))
        except Exception as e:
            _enrich_data.set(("error", str(e)))

    @output
    @render.ui
    def enrich_status():
        data = _enrich_data()
        if data is None:
            return ui.p("Enter a drug list and click Run Enrichment.")
        if data[0] == "error":
            return ui.div(data[1], class_="alert alert-danger")
        _, _, _, unresolved = data
        if unresolved:
            return ui.div(
                f"Unresolved drugs ({len(unresolved)}): "
                + ", ".join(unresolved[:10])
                + ("..." if len(unresolved) > 10 else ""),
                class_="alert alert-warning",
            )
        return ui.span()

    @output
    @render.data_frame
    def enrich_table():
        data = _enrich_data()
        if data is None or data[0] != "ok":
            return pd.DataFrame()
        _, combined, _, _ = data
        if combined.empty:
            return combined
        fdr = input.fdr_thresh()
        ovl = input.min_overlap()
        filt = combined[(combined["qvalue"] <= fdr) &
                        (combined["overlap"] >= ovl)].sort_values(["pvalue"])
        return render.DataGrid(
            filt[["library", "term", "overlap", "term_size", "pvalue", "qvalue"]],
            row_selection_mode="none",
            filters=True,
        )

    @output
    @render.data_frame
    def coverage_table():
        data = _enrich_data()
        if data is None or data[0] != "ok":
            return pd.DataFrame()
        _, _, coverage_df, _ = data
        return render.DataGrid(coverage_df, row_selection_mode="none")


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _df_to_table(df: pd.DataFrame, cols: list[str]) -> ui.TagList:
    if df.empty:
        return ui.p(ui.tags.i("(none)"))
    rows = []
    for _, row in df.iterrows():
        rows.append(ui.tags.tr(*[ui.tags.td(str(row.get(c, ""))) for c in cols]))
    header = ui.tags.tr(*[ui.tags.th(c.replace("_", " ").title()) for c in cols])
    return ui.tags.table(
        ui.tags.thead(header),
        ui.tags.tbody(*rows),
        class_="table table-sm table-striped table-hover",
    )


app = App(app_ui, server)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
