import pandas as pd
from pathlib import Path

RAW = Path("data/raw/luad_tcga_pan_can_atlas_2018")
OUT = Path("data/processed")
OUT.mkdir(exist_ok=True)

def load_cbioportal(filename, skip=4):
    return pd.read_csv(RAW / filename, sep="\t", skiprows=skip, low_memory=False)

# Clinical
patient = load_cbioportal("data_clinical_patient.txt")
sample  = load_cbioportal("data_clinical_sample.txt")

# Merge
df = patient.merge(sample, on="PATIENT_ID", how="inner")

# Columnas clave
keep = [
    "PATIENT_ID", "AGE", "SEX", "RACE", "ETHNICITY",
    "SUBTYPE",
    "AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_T_STAGE", "PATH_N_STAGE", "PATH_M_STAGE",
    "OS_STATUS", "OS_MONTHS",
    "DFS_STATUS", "DFS_MONTHS",
    "RADIATION_THERAPY", "HISTORY_NEOADJUVANT_TRTYN",
    "GENETIC_ANCESTRY_LABEL",
    "TMB_NONSYNONYMOUS", "MSI_SENSOR_SCORE", "ANEUPLOIDY_SCORE",
    "GRADE", "TUMOR_TYPE", "CANCER_TYPE_DETAILED"
]

df = df[keep].copy()

# OS: convertir status a binario
df["OS_EVENT"] = df["OS_STATUS"].apply(lambda x: 1 if "DECEASED" in str(x) else 0)
df["DFS_EVENT"] = df["DFS_STATUS"].apply(lambda x: 1 if "Recurred" in str(x) else 0)

# Stage simplificado
def simplify_stage(s):
    s = str(s).upper()
    if "IV" in s: return "IV"
    if "III" in s: return "III"
    if "II" in s: return "II"
    if "I" in s: return "I"
    return None

df["STAGE_SIMPLE"] = df["AJCC_PATHOLOGIC_TUMOR_STAGE"].apply(simplify_stage)

# Info
print("Shape:", df.shape)
print("OS events:", df["OS_EVENT"].sum(), "/", len(df))
print("DFS events:", df["DFS_EVENT"].sum(), "/", len(df))
print("Stage distribution:\n", df["STAGE_SIMPLE"].value_counts())
print("Missing OS_MONTHS:", df["OS_MONTHS"].isna().sum())

# ── Mutations ──────────────────────────────────────────────────────────────
mut = pd.read_csv(RAW / "data_mutations.txt", sep="\t", low_memory=False,
                  usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "HGVSp_Short"])

# Extract PATIENT_ID from barcode (first 12 chars)
mut["PATIENT_ID"] = mut["Tumor_Sample_Barcode"].str[:12]

# Genes of interest
GENES = ["KRAS", "TP53", "EGFR", "STK11", "KEAP1", "BRAF", "ALK", "ROS1", "MET"]

for gene in GENES:
    patients = mut[mut["Hugo_Symbol"] == gene]["PATIENT_ID"].unique()
    df[f"{gene}_MUT"] = df["PATIENT_ID"].isin(patients).astype(int)

# KRAS subtypes
kras_mut = mut[mut["Hugo_Symbol"] == "KRAS"][["PATIENT_ID", "HGVSp_Short"]].drop_duplicates()
kras_mut = kras_mut.rename(columns={"HGVSp_Short": "KRAS_SUBTYPE"})
df = df.merge(kras_mut, on="PATIENT_ID", how="left")

# Simplify KRAS subtype
def kras_group(s):
    if pd.isna(s): return "WT"
    if s == "p.G12C": return "G12C"
    if s == "p.G12V": return "G12V"
    if s == "p.G12D": return "G12D"
    if "G12" in s: return "G12other"
    return "other"

df["KRAS_GROUP"] = df["KRAS_SUBTYPE"].apply(kras_group)

print("\nKRAS group distribution:\n", df["KRAS_GROUP"].value_counts())
print("\nMutation prevalence:")
for gene in GENES:
    print(f"  {gene}: {df[f'{gene}_MUT'].sum()} ({df[f'{gene}_MUT'].mean():.1%})")



df.to_csv(OUT / "clinical_clean.csv", index=False)
print("\n✅ Saved to data/processed/clinical_clean.csv")

# ── RNA-seq data ───────────────────────────────────────────────────────────
print("\nLoading RNA-seq data...")
rna = pd.read_csv(RAW / "data_mrna_seq_v2_rsem.txt", sep="\t", index_col=0, low_memory=False)
rna = rna.drop(columns=["Entrez_Gene_Id"], errors="ignore")

# Clean gene symbols
rna = rna[rna.index.notna()]
rna = rna[rna.index != ""]
rna = rna[~rna.index.duplicated(keep="first")]

# Extract patient IDs from sample barcodes
rna.columns = [c[:12] for c in rna.columns]

# Filter low-expression genes
rna_filtered = rna[rna.median(axis=1) > 1]

# Top 1000 most variable genes
gene_var = rna_filtered.var(axis=1)
top_genes = gene_var.nlargest(1000).index
rna_top = rna_filtered.loc[top_genes]

# Transpose
rna_t = rna_top.T
rna_t.index.name = "PATIENT_ID"

# Merge with clinical
clinical_dedup = df.drop_duplicates(subset="PATIENT_ID", keep="first")
multimodal = clinical_dedup.merge(rna_t, on="PATIENT_ID", how="inner")
multimodal = multimodal[multimodal["OS_MONTHS"] > 0].dropna(subset=["OS_MONTHS", "OS_EVENT"])

print(f"Multimodal dataset: {multimodal.shape}")
print(f"OS events: {multimodal['OS_EVENT'].sum()} ({multimodal['OS_EVENT'].mean():.1%})")

multimodal.to_csv(OUT / "multimodal_clean.csv", index=False)
print("✓ Saved to data/processed/multimodal_clean.csv")
