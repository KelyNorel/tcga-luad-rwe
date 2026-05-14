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

df.to_csv(OUT / "clinical_clean.csv", index=False)
print("\n✅ Saved to data/processed/clinical_clean.csv")