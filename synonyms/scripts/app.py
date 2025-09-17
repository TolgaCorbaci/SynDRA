from shiny import App, render, ui
import pandas as pd
import re
import os

csv_path = os.path.join(os.path.dirname(__file__), "data", "merged_200K_drug_synonyms.csv")
synonyms_df = pd.read_csv(csv_path, dtype=str)


# Clean the synonyms column
synonyms_df["synonyms"] = synonyms_df["synonyms"].str.lower().str.strip()


# --- Function to find synonyms efficiently ---
def find_synonyms(drug_names):
    # Clean input
    drug_names = [d.lower().strip() for d in drug_names if d.strip()]
    if not drug_names:
        return pd.DataFrame()

    # Combine all drug names into a single regex pattern
    pattern = r"\b(" + "|".join(map(re.escape, drug_names)) + r")\b"
    matches = synonyms_df[synonyms_df["synonyms"].str.contains(pattern, regex=True, na=False)]

    return matches[["BROAD_drug_ID", "synonyms"]].drop_duplicates().reset_index(drop=True)


# --- UI ---
app_ui = ui.page_fluid(
    ui.h2("SynDRA: Synonym & BROAD_drug Finder"),
    ui.p("Enter drug names (newline or comma-separated)."),
    ui.input_text_area(
        "drugs", 
        "Drug names:", 
        placeholder="e.g., aspirin\nmetformin, ibuprofen"
    ),
    ui.input_action_button("go", "Find Synonyms"),
    ui.output_text_verbatim("match_count"),
    ui.output_table("results"),
    ui.download_button("download_results", "Download CSV"),
    ui.p(
        "GitHub repository: ",
        ui.a("SynDRA", href="https://github.com/hidelab/SynDRA", target="_blank")
    )
)


# --- Server ---
def server(input, output, session):
    cached_results = {}

    def get_results():
        # Check cache to avoid repeated computation
        drug_text = input.drugs()
        if drug_text in cached_results:
            return cached_results[drug_text]

        drug_list = re.split(r"[\n,]+", drug_text)
        df = find_synonyms(drug_list)
        cached_results[drug_text] = df
        return df

    @output
    @render.table
    def results():
        if input.go() == 0:
            return pd.DataFrame()
        return get_results()

    @output
    @render.text
    def match_count():
        if input.go() == 0:
            return ""
        df = get_results()
        if df.empty:
            return "No matches found."
        return f"{len(df)} matches found."

    @output
    @render.download
    def download_results():
        df = get_results()
        return df.to_csv(index=False)


# --- App object ---
app = App(app_ui, server)
