import os
import subprocess
from pathlib import Path
import pandas as pd
import streamlit as st


# ---------- Constants ----------
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data_set"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

PPI = DATA / "ppi_network/HIPPIE.tsv"
ONTO = DATA / "ontology/ontology_graph.txt"
ONCOKB = ROOT / "Dataset OncoKB.xlsx"

CANCERS = ["BRCA", "COAD", "LUAD", "THCA", "BLCA", "PRAD", "STAD"]

ABLATIONS = {
    "FULL": dict(use_c=True, use_de=True, use_onto=True),
    "noCOEXP": dict(use_c=False, use_de=True, use_onto=True),
    "noDE": dict(use_c=True, use_de=False, use_onto=True),
    "noONTO": dict(use_c=True, use_de=True, use_onto=False),
    "PPI_only": dict(use_c=False, use_de=False, use_onto=False),
}


# ---------- Helpers ----------
def build_cmd_main(ppi, seed, alpha, beta, restart_prob, out_txt, cfg, coexp=None, de=None, onto=None, diso=None):
    cmd = [
        "python",
        str(ROOT / "main.py"),
        "-p",
        str(ppi),
        "-s",
        str(seed),
        "-x",
        str(alpha),
        "-y",
        str(beta),
        "-r",
        str(restart_prob),
        "-o",
        str(out_txt),
    ]
    if cfg.get("use_c") and coexp is not None:
        cmd += ["-c", str(coexp)]
    if cfg.get("use_de") and de is not None:
        cmd += ["-de", str(de)]
    if cfg.get("use_onto") and onto is not None and diso is not None:
        cmd += ["-a", str(onto), "-do", str(diso)]
    return cmd


def run_once(cancer: str, tag: str, alpha: float, beta: float, cfg: dict) -> int:
    out_cancer = OUT / cancer
    out_cancer.mkdir(exist_ok=True)

    out_txt = out_cancer / f"results_{tag}.txt"
    out_tsv = out_cancer / f"top100_{tag}.tsv"

    seed = DATA / f"seed_set/TCGA-{cancer}_seed.txt"
    de = DATA / f"differentially_expressed_genes/TCGA-{cancer}_de_genes.tsv"
    coexp = DATA / f"co-expression_networks/TCGA-{cancer}__co_expression__t_70%.tsv"
    diso = DATA / f"disease_specific_ontologies/TCGA-{cancer}_disease_ontologies.txt"

    # sanity checks
    if not PPI.exists():
        raise FileNotFoundError(f"Missing PPI: {PPI}")
    if not seed.exists():
        raise FileNotFoundError(f"Missing seed set: {seed}")
    if cfg.get("use_c") and not coexp.exists():
        raise FileNotFoundError(f"Missing co-expression: {coexp}")
    if cfg.get("use_de") and not de.exists():
        raise FileNotFoundError(f"Missing DE genes: {de}")
    if cfg.get("use_onto"):
        if not ONTO.exists():
            raise FileNotFoundError(f"Missing ontology graph: {ONTO}")
        if not diso.exists():
            raise FileNotFoundError(f"Missing disease-specific ontologies: {diso}")

    cmd = build_cmd_main(
        PPI,
        seed,
        alpha,
        beta,
        restart_prob,
        out_txt=out_txt,
        cfg=cfg,
        coexp=coexp,
        de=de,
        onto=ONTO,
        diso=diso,
    )

    subprocess.run(cmd, cwd=str(ROOT), check=True)

    if not ONCOKB.exists():
        raise FileNotFoundError(f"Missing OncoKB Excel: {ONCOKB}")

    subprocess.run(
        [
            "python",
            str(ROOT / "check_genes.py"),
            "--input",
            str(out_txt),
            "--oncokb",
            str(ONCOKB),
            "--output",
            str(out_tsv),
        ],
        cwd=str(ROOT),
        check=True,
    )

    df = pd.read_csv(out_tsv, sep="\t")
    return int((df["In_OncoKB"] == "Yes").sum())


@st.cache_data
def load_summary(summary_path: Path) -> pd.DataFrame:
    if not summary_path.exists():
        return pd.DataFrame(columns=["Type", "Cancer", "Tag", "Alpha", "Beta", "OncoKB_hits"])
    df = pd.read_csv(summary_path)
    # enforce dtypes
    df["Alpha"] = pd.to_numeric(df["Alpha"], errors="coerce")
    df["Beta"] = pd.to_numeric(df["Beta"], errors="coerce")
    df["OncoKB_hits"] = pd.to_numeric(df["OncoKB_hits"], errors="coerce").fillna(0).astype(int)
    return df


# ---------- UI ----------
st.set_page_config(page_title="BRW GUI", layout="wide")
st.title("Biological Random Walks")

mode = st.sidebar.radio("Mode", ["Run BRW now", "Load from summary.csv"]) 

st.sidebar.markdown("---")
st.sidebar.markdown("**Input Settings**", help="Configure input data and algorithm parameters")
cancer = st.sidebar.selectbox("Cancer", CANCERS, index=0)
selected_opts = st.sidebar.multiselect(
    "Select options (ablation)", list(ABLATIONS.keys()), default=["FULL", "noCOEXP", "noDE", "noONTO", "PPI_only"],
)

# Restart probability
restart_prob = st.sidebar.number_input("Restart probability (r)", min_value=0.1, max_value=0.99, value=0.9, step=0.05)

# Separator between options and alpha-beta
st.sidebar.markdown("---")
st.sidebar.markdown("**Matrix Weights**", help="Configure alpha-beta weights for network combination")

# Initialize session state for alpha-beta pairs
if 'alpha_beta_pairs' not in st.session_state:
    st.session_state.alpha_beta_pairs = [(0.5, 0.5), (0.25, 0.75), (0.75, 0.25)]

# Quick presets

col1, col2 = st.sidebar.columns(2)

with col1:
    if st.button("Load", use_container_width=True, help="Load 7 common alpha-beta pairs"):
        st.session_state.alpha_beta_pairs = [(1.0, 1.0), (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.5), (0.75, 0.25), (0.25, 0.75)]

with col2:
    if st.button("Reset", use_container_width=True, help="Reset to 3 default pairs"):
        st.session_state.alpha_beta_pairs = [(0.5, 0.5), (0.25, 0.75), (0.75, 0.25)]

# Number of pairs control
num_pairs = st.sidebar.number_input("Number of alpha-beta pairs", min_value=1, max_value=10, value=len(st.session_state.alpha_beta_pairs), step=1, help="Adjust number of alpha-beta pairs")

# Ensure we have enough pairs
while len(st.session_state.alpha_beta_pairs) < num_pairs:
    st.session_state.alpha_beta_pairs.append((0.5, 0.5))

# Trim if too many
if len(st.session_state.alpha_beta_pairs) > num_pairs:
    st.session_state.alpha_beta_pairs = st.session_state.alpha_beta_pairs[:num_pairs]

# Display and edit pairs
for i in range(num_pairs):
    col1, col2 = st.sidebar.columns(2)
    with col1:
        alpha = col1.number_input(f"α{i+1}", min_value=0.0, max_value=1.0, value=st.session_state.alpha_beta_pairs[i][0], step=0.05, key=f"alpha_{i}", help=f"Alpha value for pair {i+1}")
    with col2:
        beta = col2.number_input(f"β{i+1}", min_value=0.0, max_value=1.0, value=st.session_state.alpha_beta_pairs[i][1], step=0.05, key=f"beta_{i}", help=f"Beta value for pair {i+1}")
    
    # Update session state
    st.session_state.alpha_beta_pairs[i] = (alpha, beta)

# Use the session state
alpha_beta_pairs = st.session_state.alpha_beta_pairs


if mode == "Run BRW now":
    st.subheader("Run BRW for 1 cancer with multiple options and alpha-beta")
    st.write(f"App will run {len(selected_opts)} options × {len(alpha_beta_pairs)} alpha-beta pairs = {len(selected_opts) * len(alpha_beta_pairs)} combinations")
    st.write(f"**Restart probability:** r = {restart_prob:.2f}")
    
    # Show selected pairs
    st.write("**Selected alpha-beta pairs:**")
    for i, (a, b) in enumerate(alpha_beta_pairs):
        st.write(f"Pair {i+1}: α={a:.2f}, β={b:.2f}")
    
    run_btn = st.button("Run")

    if run_btn:
        results = []
        total_combinations = len(selected_opts) * len(alpha_beta_pairs)
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        with st.spinner("Running..."):
            current = 0
            for name in selected_opts:
                cfg = ABLATIONS[name]
                for alpha, beta in alpha_beta_pairs:
                    current += 1
                    status_text.text(f"Running {name} with α={alpha:.2f}, β={beta:.2f}, r={restart_prob:.2f} ({current}/{total_combinations})")
                    progress_bar.progress(current / total_combinations)
                    
                    # Create tag based on option and alpha-beta
                    if name == "FULL":
                        tag = f"FULL_A{alpha}_B{beta}"
                    else:
                        tag = f"{name}_A{alpha}_B{beta}"
                    
                    try:
                        hits = run_once(cancer, tag, alpha, beta, cfg)
                        results.append(dict(
                            Cancer=cancer, 
                            Option=name,
                            Tag=tag, 
                            Alpha=round(alpha, 2), 
                            Beta=round(beta, 2), 
                            Restart_prob=round(restart_prob, 2),
                            OncoKB_hits=hits
                        ))
                    except Exception as e:
                        st.error(f"{name} (alpha={alpha:.2f}, beta={beta:.2f}): {e}")

        if results:
            df = pd.DataFrame(results)
            st.success("Completed!")
            
            # Display results table
            st.subheader("Detailed Results")
            st.dataframe(df, use_container_width=True)
            
            # Comparison charts
            st.subheader("OncoKB Hits Comparison")
            
            # Chart 1: By option (averaging across alpha-beta)
            option_avg = df.groupby('Option')['OncoKB_hits'].mean().reset_index()
            st.write("**Average OncoKB hits by ablation:**")
            st.bar_chart(option_avg.set_index("Option")["OncoKB_hits"])
            
            # Chart 2: By alpha-beta pair (averaging across options)
            pair_avg = df.groupby(['Alpha', 'Beta'])['OncoKB_hits'].mean().reset_index()
            pair_avg['Pair'] = pair_avg.apply(lambda x: f"α={x['Alpha']:.2f}, β={x['Beta']:.2f}", axis=1)
            st.write("**Average OncoKB hits by alpha-beta:**")
            st.bar_chart(pair_avg.set_index("Pair")["OncoKB_hits"])
            
            # Chart 3: Heatmap view
            st.write("**Heatmap OncoKB hits (Option × Alpha-Beta):**")
            pivot_df = df.pivot(index='Option', columns=['Alpha', 'Beta'], values='OncoKB_hits')
            st.dataframe(pivot_df, use_container_width=True)
            
            # Chart 4: Restart probability analysis
            if len(df['Restart_prob'].unique()) > 1:
                restart_avg = df.groupby('Restart_prob')['OncoKB_hits'].mean().reset_index()
                st.write("**Average OncoKB hits by restart probability:**")
                st.bar_chart(restart_avg.set_index("Restart_prob")["OncoKB_hits"])

else:
    st.subheader("Quick load from outputs/summary.csv")
    summary_path = OUT / "summary.csv"
    df_sum = load_summary(summary_path)
    if df_sum.empty:
        st.warning(f"No summary file found: {summary_path}")
    else:
        df_view = df_sum[df_sum["Cancer"] == cancer]
        
        if not df_view.empty:
            # Filter by selected options if any
            if selected_opts:
                df_view = df_view[df_view["Tag"].str.contains("|".join(selected_opts), case=False, na=False)]
            
            # Filter by alpha-beta pairs if any
            if alpha_beta_pairs:
                # Extract alpha-beta from tags like "FULL_A0.5_B0.5"
                def extract_alpha_beta(tag):
                    import re
                    match = re.search(r'A([\d.]+)_B([\d.]+)', tag)
                    if match:
                        return float(match.group(1)), float(match.group(2))
                    return None, None
                
                df_view['Extracted_Alpha'] = df_view['Tag'].apply(lambda x: extract_alpha_beta(x)[0])
                df_view['Extracted_Beta'] = df_view['Tag'].apply(lambda x: extract_alpha_beta(x)[1])
                
                # Filter by selected alpha-beta pairs
                mask = df_view.apply(lambda row: 
                    (row['Extracted_Alpha'], row['Extracted_Beta']) in alpha_beta_pairs, axis=1)
                df_view = df_view[mask]
            
            st.dataframe(df_view, use_container_width=True)
            
            if not df_view.empty:
                # Multiple visualization options
                st.subheader("Comparison Charts")
                
                # Chart 1: By option (if we have multiple options)
                if 'Option' in df_view.columns or df_view['Tag'].str.contains('FULL|noCOEXP|noDE|noONTO|PPI_only').any():
                    option_col = 'Option' if 'Option' in df_view.columns else 'Tag'
                    option_avg = df_view.groupby(option_col)['OncoKB_hits'].mean().reset_index()
                    st.write("**Average OncoKB hits by option:**")
                    st.bar_chart(option_avg.set_index(option_col)["OncoKB_hits"])
                
                # Chart 2: By alpha-beta pair
                if 'Alpha' in df_view.columns and 'Beta' in df_view.columns:
                    pair_avg = df_view.groupby(['Alpha', 'Beta'])['OncoKB_hits'].mean().reset_index()
                    pair_avg['Pair'] = pair_avg.apply(lambda x: f"α={x['Alpha']:.2f}, β={x['Beta']:.2f}", axis=1)
                    st.write("**Average OncoKB hits by alpha-beta pair:**")
                    st.bar_chart(pair_avg.set_index("Pair")["OncoKB_hits"])
                
                # Chart 3: Simple tag-based chart
                st.write("**OncoKB hits by tag:**")
                st.bar_chart(df_view.set_index("Tag")["OncoKB_hits"])
        else:
            st.info("No data for this cancer or with selected filters")

st.markdown("---")

# Summary statistics
if mode == "Run BRW now" and 'results' in locals() and results:
    st.subheader("📊 Summary Statistics")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        best_option = max(results, key=lambda x: x['OncoKB_hits'])
        st.metric("Best Option", f"{best_option['Option']} (α={best_option['Alpha']:.2f}, β={best_option['Beta']:.2f}, r={best_option['Restart_prob']:.2f})", 
                 f"{best_option['OncoKB_hits']} hits")
    
    with col2:
        avg_hits = sum(r['OncoKB_hits'] for r in results) / len(results)
        st.metric("Average OncoKB hits", f"{avg_hits:.1f}")
    
    with col3:
        total_runs = len(results)
        st.metric("Total runs", total_runs)
    
    # Export results
    st.subheader("💾 Export Results")
    csv_data = df.to_csv(index=False)
    st.download_button(
        label="Download CSV",
        data=csv_data,
        file_name=f"BRW_results_{cancer}_{len(alpha_beta_pairs)}pairs.csv",
        mime="text/csv"
    )

elif mode == "Load from summary.csv" and 'df_view' in locals() and not df_view.empty:
    st.subheader("📊 Statistics from summary.csv")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        best_tag = df_view.loc[df_view['OncoKB_hits'].idxmax()]
        st.metric("Best Result", best_tag['Tag'], f"{best_tag['OncoKB_hits']} hits")
    
    with col2:
        avg_hits = df_view['OncoKB_hits'].mean()
        st.metric("Average OncoKB hits", f"{avg_hits:.1f}")
    
    with col3:
        total_results = len(df_view)
        st.metric("Total results", total_results)

st.markdown("---")
st.caption("💡 **Tip:** Run with `streamlit run app_gui.py` in the `BiologicalRandomWalks` directory. Requirements: `pandas`, `streamlit`, `mygene`, `openpyxl`.")


