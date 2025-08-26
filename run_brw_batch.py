import subprocess, pathlib, csv
import pandas as pd

# ------- 1. Thư mục dữ liệu và đầu ra -------
DATA = pathlib.Path("data_set")
OUT  = pathlib.Path("outputs")
OUT.mkdir(exist_ok=True)

PPI   = DATA / "ppi_network/HIPPIE.tsv"
ONTO  = DATA / "ontology/ontology_graph.txt"
ONCOKB = pathlib.Path("Dataset OncoKB.xlsx")

# ------- 2. Danh sách các loại cancer -------
CANCERS = ["BRCA", "COAD", "LUAD", "THCA", "BLCA", "PRAD", "STAD"]

# ------- 3. Danh sách ablation -------
ablations = {
    "FULL"      : dict(use_c=True, use_de=True, use_onto=True),
    "noCOEXP"   : dict(use_c=False, use_de=True, use_onto=True),
    "noDE"      : dict(use_c=True,  use_de=False, use_onto=True),
    "noONTO"    : dict(use_c=True,  use_de=True, use_onto=False),
    "PPI_only"  : dict(use_c=False, use_de=False, use_onto=False),
}

# ------- 4. Lưới alpha-beta -------
alpha_beta_grid = [
    (1.0, 1.0), (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.5),
    (0.75, 0.25), (0.25, 0.75), (0.8, 0.2), (0.2, 0.8), (0.25, 0.25)
]

# ------- 5. Hàm chạy một lần -------
def run_once(cancer, tag, alpha, beta, cfg):
    OUT_CANCER = OUT / cancer
    OUT_CANCER.mkdir(exist_ok=True)

    out_txt = OUT_CANCER / f"results_{tag}.txt"
    out_tsv = OUT_CANCER / f"top100_{tag}.tsv"

    SEED  = DATA / f"seed_set/TCGA-{cancer}_seed.txt"
    DE    = DATA / f"differentially_expressed_genes/TCGA-{cancer}_de_genes.tsv"
    COEXP = DATA / f"co-expression_networks/TCGA-{cancer}__co_expression__t_70%.tsv"
    DISO  = DATA / f"disease_specific_ontologies/TCGA-{cancer}_disease_ontologies.txt"

    cmd = ["python", "main.py",
           "-p", str(PPI), "-s", str(SEED),
           "-x", str(alpha), "-y", str(beta), "-r", "0.9",
           "-o", str(out_txt)]

    if cfg["use_c"]:   cmd += ["-c", str(COEXP)]
    if cfg["use_de"]:  cmd += ["-de", str(DE)]
    if cfg["use_onto"]: cmd += ["-a", str(ONTO), "-do", str(DISO)]

    subprocess.run(cmd, check=True)

    subprocess.run([
        "python", "check_genes.py",
        "--input", str(out_txt),
        "--oncokb", str(ONCOKB),
        "--output", str(out_tsv)
    ], check=True)

    df = pd.read_csv(out_tsv, sep="\t")
    return (df["In_OncoKB"] == "Yes").sum()

# ------- 6. Chạy toàn bộ -------
summary = []

for cancer in CANCERS:
    print(f"\n🧬 Cancer: {cancer}")
    for name, cfg in ablations.items():
        hits = run_once(cancer, name, 0.5, 0.5, cfg)
        summary.append(("ABL", cancer, name, 0.5, 0.5, hits))

    for (a, b) in alpha_beta_grid:
        tag = f"FULL_A{a}_B{b}"
        hits = run_once(cancer, tag, a, b, ablations["FULL"])
        summary.append(("Tuning(alpha,beta)", cancer, tag, a, b, hits))

# ------- 7. Lưu tổng kết -------
with open(OUT/"summary.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["Type", "Cancer", "Tag", "Alpha", "Beta", "OncoKB_hits"])
    w.writerows(summary)

print("\n✅ Hoàn tất. Xem thư mục outputs/<Cancer>/ và outputs/summary.csv")