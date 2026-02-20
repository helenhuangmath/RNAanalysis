"""
Build RNAseq_pipeline_tutorial.ipynb programmatically.
Run:  python notebooks/make_notebook.py
"""
import json, os, textwrap, uuid

def cell_id():
    return uuid.uuid4().hex[:16]

def md(source: str) -> dict:
    return {
        "cell_type": "markdown",
        "id": cell_id(),
        "metadata": {},
        "source": source.strip(),
    }

def code(source: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": cell_id(),
        "metadata": {},
        "outputs": [],
        "source": textwrap.dedent(source).strip(),
    }

# =============================================================================
cells = []
# =============================================================================

cells.append(md(r"""
# RNA-seq Data Analysis Pipeline Tutorial
### FastQC / MultiQC → STAR → featureCounts → TPM → DESeq2 → Pathway/GSEA → Gene Network → TF Binding

---

| Step | Tool | Purpose |
|------|------|---------|
| 1 | FastQC / MultiQC | Raw read quality control |
| 2 | STAR (2-pass) | Splice-aware genome alignment |
| 3 | featureCounts | Gene-level read quantification |
| 4 | TPM | Library-size + length normalization |
| 5 | pyDESeq2 | Differential expression testing |
| 6 | gseapy (Enrichr) | Over-representation analysis (ORA) |
| 7 | gseapy (Prerank) | Gene set enrichment analysis (GSEA) |
| 8 | STRING + networkx | Gene–gene interaction network |
| 9 | ChEA3 + decoupler | Transcription factor binding prediction |

**Input :** `data/raw_counts.csv`, `data/gene_lengths.csv`, `data/metadata.csv`
**Output:** `results/tables/` (DE results) and `results/figures/` (all plots)

> **Tutorial mode:** Section 3 simulates a realistic RNA-seq dataset so every cell
> runs end-to-end without real FASTQ files. Replace the simulated data with your
> own `raw_counts.csv` / `gene_lengths.csv` to analyse real data.
"""))

# ── Dependencies ──────────────────────────────────────────────────────────────
cells.append(md("## 0. Install & Import Dependencies"))

cells.append(code("""
# Run this cell once to install required packages
import subprocess, sys

pkgs = [
    "pydeseq2", "gseapy", "networkx", "decoupler",
    "matplotlib", "seaborn", "pandas", "numpy", "scipy",
    "scikit-learn", "adjustText", "requests",
]
for p in pkgs:
    subprocess.check_call([sys.executable, "-m", "pip", "install", p, "-q"],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
print("All packages installed.")
"""))

cells.append(code("""
import warnings
warnings.filterwarnings("ignore")

import os, json, itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import networkx as nx
import requests

# ── Plotting defaults ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})
sns.set_style("white")
sns.set_context("notebook")

COND_COLORS = {"Control": "#4C72B0", "Treatment": "#DD8452"}

# ── Output directories ────────────────────────────────────────────────────────
for d in ["data", "results/figures", "results/tables"]:
    os.makedirs(d, exist_ok=True)

print("Setup complete.")
"""))

# ── Section 1: Metadata ───────────────────────────────────────────────────────
cells.append(md("""
## 1. Sample Metadata

The metadata CSV maps each sample to its experimental condition and covariates.
Every downstream step uses this file to label plots and define contrasts.
"""))

cells.append(code("""
# ── Create / load metadata ────────────────────────────────────────────────────
meta_path = "data/metadata.csv"

if not os.path.exists(meta_path):
    metadata = pd.DataFrame({
        "sample_id" : ["ctrl_1","ctrl_2","ctrl_3","ctrl_4",
                        "trt_1", "trt_2", "trt_3", "trt_4"],
        "condition" : ["Control"]*4 + ["Treatment"]*4,
        "batch"     : ["A","A","B","B","A","A","B","B"],
        "sex"       : ["M","F","M","F","M","F","M","F"],
    })
    metadata.set_index("sample_id", inplace=True)
    metadata.to_csv(meta_path)
else:
    metadata = pd.read_csv(meta_path, index_col=0)

SAMPLES    = metadata.index.tolist()
N_SAMPLES  = len(SAMPLES)

print(f"Metadata loaded: {N_SAMPLES} samples")
print(metadata.to_string())
"""))

# ── Section 2: Upstream shell commands ────────────────────────────────────────
cells.append(md("""
## 2. Upstream Processing (Shell)

> **Skip if you already have `data/raw_counts.csv` and `data/gene_lengths.csv`.**
> The full shell pipeline is in `scripts/upstream_processing.sh`.

The cells below illustrate the key commands; run them in a terminal or Jupyter
`%%bash` cells after installing the required bioinformatics tools.

### 2.1 Quality Control — FastQC / MultiQC
```bash
# Per-sample FastQC
fastqc data/fastq/*.fastq.gz --outdir results/qc/fastqc_raw --threads 8

# Aggregate report
multiqc results/qc/fastqc_raw --outdir results/qc --filename multiqc_raw
```

### 2.2 Splice-Aware Alignment — STAR (2-pass)
```bash
# Build genome index (once)
STAR --runMode genomeGenerate \\
     --genomeDir genome/STAR_index \\
     --genomeFastaFiles genome/genome.fa \\
     --sjdbGTFfile genome/annotation.gtf \\
     --sjdbOverhang 149 --runThreadN 8

# Align each sample
STAR --runThreadN 8 \\
     --genomeDir genome/STAR_index \\
     --readFilesIn data/fastq/ctrl_1_R1.fastq.gz data/fastq/ctrl_1_R2.fastq.gz \\
     --readFilesCommand zcat \\
     --outSAMtype BAM SortedByCoordinate \\
     --quantMode GeneCounts \\
     --outFileNamePrefix results/aligned/ctrl_1/
samtools index results/aligned/ctrl_1/Aligned.sortedByCoord.out.bam
```

### 2.3 Gene-Level Quantification — featureCounts
```bash
featureCounts -T 8 -p --countReadPairs -s 2 -B -C \\
    -a genome/annotation.gtf \\
    -o results/counts/featureCounts_raw.txt \\
    results/aligned/*/Aligned.sortedByCoord.out.bam
```

**Key flag:** `-s 2` = reverse-stranded (dUTP/TruSeq kits).
Use `infer_experiment.py` (RSeQC) to verify strandedness before counting.
"""))

# ── Section 3: Simulate data ──────────────────────────────────────────────────
cells.append(md("""
## 3. Simulate Tutorial Dataset

We generate a realistic negative-binomial RNA-seq count matrix so every
downstream cell executes without real FASTQ files.

**Design:** 4 Control vs 4 Treatment samples, ~15 000 genes,
~160 true up-regulated DEGs (immune/apoptosis) and ~100 true down-regulated DEGs
(cell cycle/proliferation), remainder background.
"""))

cells.append(code("""
np.random.seed(42)

# ── True DEG lists — real HGNC symbols for meaningful pathway/TF results ──────
UP_GENES = [
    # Interferon / antiviral
    "IFNG","STAT1","STAT2","IRF1","IRF3","IRF7","IRF9",
    "ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","OASL",
    "IFIT1","IFIT2","IFIT3","IFIT5","IFITM1","IFITM2","IFITM3",
    "RSAD2","HERC5","HERC6","USP18","TRIM25","DDX58","IFIH1",
    # Cytokines / chemokines
    "IL6","IL1A","IL1B","IL12A","IL12B","IL15","IL18",
    "TNF","TNFSF10","IFNB1","CXCL9","CXCL10","CXCL11",
    "CCL2","CCL5","CCL7","CXCL8",
    # NF-κB pathway
    "NFKB1","NFKB2","RELA","RELB","REL","NFKBIA","NFKBIZ",
    "TLR3","TLR4","TLR7","MYD88","TICAM1","RIPK1",
    # Apoptosis / p53
    "TP53","BAX","BAK1","BBC3","PMAIP1","BID","CYCS","APAF1",
    "CASP3","CASP8","CASP9","CASP7","FAS","FASLG",
    "CDKN1A","GADD45A","GADD45B","MDM2",
    # Cytotoxic T-cell / NK markers
    "CD8A","GZMB","GZMK","PRF1","GNLY","NKG7","CD274",
    "PDCD1","HAVCR2","LAG3","TIGIT",
]

DOWN_GENES = [
    # Cell cycle — CDKs / Cyclins
    "CDK1","CDK2","CDK4","CDK6",
    "CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCNE1","CCNE2",
    "E2F1","E2F2","E2F3","MYBL2",
    # Replication / proliferation markers
    "MKI67","PCNA","MCM2","MCM4","MCM7","GINS1","CDC45",
    "RRM1","RRM2","TYMS","DHFR","TK1",
    # Mitosis regulators
    "AURKB","AURKA","PLK1","PLK4","BUB1","BUB1B",
    "MAD2L1","CDC20","CDH1","TOP2A","PTTG1",
    "CENPE","CENPF","KIF11","KIF20A","RACGAP1","PRC1","BIRC5",
    # Pro-survival / anti-apoptosis
    "BCL2","BCL2L1","MCL1","XIAP","BIRC3",
    # PI3K / AKT / mTOR
    "AKT1","AKT2","AKT3","PIK3CA","PIK3CB",
    "MTOR","RPS6KB1","EIF4EBP1",
    # MYC targets / ribosome biogenesis
    "MYC","MYCN","NPM1","NCL","TERT","NME1",
    # Warburg / hypoxia metabolism
    "LDHA","PKM","HK2","GPI","PFKL","ALDOA","ENO1",
    "SLC2A1","VEGFA","CA9","HIF1A",
]

# Remove any overlap
DOWN_GENES = [g for g in DOWN_GENES if g not in UP_GENES]

ALL_DEG  = list(dict.fromkeys(UP_GENES + DOWN_GENES))
BG_GENES = [f"ENSG{i:011d}" for i in range(1, 14601)]
ALL_GENES = ALL_DEG + BG_GENES
N_GENES  = len(ALL_GENES)

print(f"Up-regulated true DEGs  : {len(UP_GENES)}")
print(f"Down-regulated true DEGs: {len(DOWN_GENES)}")
print(f"Background genes        : {len(BG_GENES):,}")
print(f"Total genes             : {N_GENES:,}")
print(f"Total samples           : {N_SAMPLES}")
"""))

cells.append(code("""
# ── Negative-binomial count simulation ────────────────────────────────────────
print("Simulating count matrix …")

# Mean expression in Control (log-normal → realistic range)
base_mean = np.random.lognormal(mean=3.8, sigma=2.0, size=N_GENES).clip(0.5, 8000)

# Per-sample library-size factors
lib_factors = np.array([1.15, 0.90, 1.05, 0.92, 1.20, 0.85, 0.98, 1.10])

# Per-gene dispersion  (larger dispersion = noisier)
dispersion = np.random.exponential(scale=0.15, size=N_GENES).clip(0.01, 2.0)

# Simulate gene lengths (bp) — log-normal, clipped to realistic range
gene_lengths_bp = np.random.lognormal(mean=7.5, sigma=1.0, size=N_GENES).astype(int)
gene_lengths_bp = np.clip(gene_lengths_bp, 200, 100_000)

counts_mat = np.zeros((N_GENES, N_SAMPLES), dtype=np.int32)

for i, gene in enumerate(ALL_GENES):
    mu_ctrl = base_mean[i]
    mu_trt  = mu_ctrl

    if gene in UP_GENES:
        lfc = np.random.uniform(1.0, 3.5)
        mu_trt = mu_ctrl * (2 ** lfc)
    elif gene in DOWN_GENES:
        lfc = np.random.uniform(1.0, 3.5)
        mu_trt = mu_ctrl / (2 ** lfc)

    for j in range(N_SAMPLES):
        mu = (mu_ctrl if j < 4 else mu_trt) * lib_factors[j]
        d  = dispersion[i]
        r  = max(1.0 / d, 0.05)
        p  = r / (r + mu)
        p  = np.clip(p, 1e-8, 1 - 1e-8)
        counts_mat[i, j] = np.random.negative_binomial(r, p)

# ── DataFrames ────────────────────────────────────────────────────────────────
counts_df = pd.DataFrame(counts_mat, index=ALL_GENES, columns=SAMPLES)
counts_df.index.name = "gene_id"

gene_lengths = pd.Series(gene_lengths_bp, index=ALL_GENES, name="length")
gene_lengths.index.name = "gene_id"

# Save
counts_df.to_csv("data/raw_counts.csv")
gene_lengths.to_csv("data/gene_lengths.csv", header=True)

print(f"Saved: data/raw_counts.csv  ({counts_df.shape[0]:,} genes × {counts_df.shape[1]} samples)")
print(f"Saved: data/gene_lengths.csv")
print(f"Library sizes: {counts_df.sum().min()/1e6:.1f}M – {counts_df.sum().max()/1e6:.1f}M reads")
"""))

# ── Section 4: Load ───────────────────────────────────────────────────────────
cells.append(md("""
## 4. Load Count Matrix

```python
# To use your own data, replace the simulated files:
#   data/raw_counts.csv   — gene_id × sample matrix of integer counts
#   data/gene_lengths.csv — gene_id, length (bp)  [from featureCounts column 5]
#   data/metadata.csv     — sample_id, condition, [batch, ...]
```
"""))

cells.append(code("""
counts_df     = pd.read_csv("data/raw_counts.csv",   index_col=0)
gene_lengths  = pd.read_csv("data/gene_lengths.csv", index_col=0).squeeze()
metadata      = pd.read_csv("data/metadata.csv",     index_col=0)

# Align sample order
counts_df = counts_df[metadata.index]

print("Count matrix :", counts_df.shape)
print("Gene lengths :", gene_lengths.shape)
print("Metadata     :", metadata.shape)
print()
display(counts_df.iloc[:5, :])
"""))

# ── Section 5: QC ─────────────────────────────────────────────────────────────
cells.append(md("""
## 5. Count-Level Quality Control

Check three indicators before normalisation:
1. **Library sizes** — total mapped reads per sample
2. **Genes detected** — number of genes with count > 0
3. **Sample correlation** — are replicates more similar than cross-condition pairs?
"""))

cells.append(code("""
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

sample_colors = [COND_COLORS[metadata.loc[s, "condition"]] for s in SAMPLES]

# 1. Library sizes
lib_sizes = counts_df.sum() / 1e6
axes[0].bar(SAMPLES, lib_sizes, color=sample_colors, edgecolor="white", linewidth=0.5)
axes[0].axhline(lib_sizes.mean(), color="red", ls="--", alpha=0.7,
                label=f"Mean {lib_sizes.mean():.1f} M")
axes[0].set_ylabel("Mapped reads (millions)")
axes[0].set_title("Library Sizes", fontweight="bold")
axes[0].tick_params(axis="x", rotation=45)
axes[0].legend(fontsize=9)

# 2. Genes detected
detected = (counts_df > 0).sum()
axes[1].bar(SAMPLES, detected, color=sample_colors, edgecolor="white", linewidth=0.5)
axes[1].set_ylabel("Genes detected (count > 0)")
axes[1].set_title("Gene Detection", fontweight="bold")
axes[1].tick_params(axis="x", rotation=45)

# 3. Sample correlation (log2 counts)
log_raw = np.log2(counts_df + 1)
corr    = log_raw.corr()
sns.heatmap(corr, ax=axes[2], cmap="RdYlBu_r", vmin=0.85, vmax=1.0,
            annot=True, fmt=".3f", annot_kws={"size": 8},
            square=True, cbar_kws={"shrink": 0.8},
            xticklabels=SAMPLES, yticklabels=SAMPLES)
axes[2].set_title("Sample Pearson Correlation\n(log₂ counts)", fontweight="bold")
axes[2].tick_params(axis="x", rotation=45)
axes[2].tick_params(axis="y", rotation=0)

# Legend
for ax in axes[:2]:
    ax.legend(handles=[
        mpatches.Patch(color=COND_COLORS["Control"],   label="Control"),
        mpatches.Patch(color=COND_COLORS["Treatment"], label="Treatment"),
    ], fontsize=9)

plt.tight_layout()
plt.savefig("results/figures/01_qc_summary.pdf", bbox_inches="tight")
plt.savefig("results/figures/01_qc_summary.png", bbox_inches="tight", dpi=300)
plt.show()
print("Saved → results/figures/01_qc_summary.pdf/.png")
"""))

# ── Section 6: TPM ────────────────────────────────────────────────────────────
cells.append(md("""
## 6. TPM Normalization

**TPM (Transcripts Per Million)** normalises for:
1. **Gene length** — longer genes accumulate more fragments
2. **Sequencing depth** — samples with more total reads give larger raw counts

$$\\text{TPM}_i = \\frac{X_i / L_i}{\\sum_j X_j / L_j} \\times 10^6$$

TPM is suitable for visualisation and within-study comparisons.
Differential expression uses raw counts via DESeq2's own normalisation.
"""))

cells.append(code("""
def compute_tpm(counts: pd.DataFrame, lengths_bp: pd.Series) -> pd.DataFrame:
    \"\"\"
    Compute TPM.

    Parameters
    ----------
    counts     : genes × samples DataFrame of integer raw counts
    lengths_bp : per-gene effective length in base pairs (e.g. from featureCounts)

    Returns
    -------
    genes × samples DataFrame of TPM values
    \"\"\"
    lengths_kb = lengths_bp.reindex(counts.index) / 1_000          # bp → kb
    rpk        = counts.div(lengths_kb, axis=0)                    # reads per kb
    scale      = rpk.sum(axis=0) / 1e6                             # per-million factor
    tpm        = rpk.div(scale, axis=1)
    return tpm


tpm_df = compute_tpm(counts_df, gene_lengths)

# Filter lowly expressed genes: keep genes with TPM ≥ 1 in ≥ min_samples
MIN_TPM      = 1.0
MIN_SAMPLES  = 2
keep         = (tpm_df >= MIN_TPM).sum(axis=1) >= MIN_SAMPLES
tpm_filt     = tpm_df.loc[keep].copy()
counts_filt  = counts_df.loc[keep].copy()

# Log2-TPM for visualisation
log2_tpm = np.log2(tpm_filt + 1)

print(f"Genes before filtering : {counts_df.shape[0]:,}")
print(f"Genes after  filtering : {counts_filt.shape[0]:,}  "
      f"(TPM ≥ {MIN_TPM} in ≥ {MIN_SAMPLES} samples)")
print()
print("TPM summary (filtered):")
print(tpm_filt.describe().round(1).to_string())
"""))

cells.append(code("""
# Distribution plots: raw log2-counts vs log2-TPM
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for j, sample in enumerate(SAMPLES):
    c = COND_COLORS[metadata.loc[sample, "condition"]]

    raw_vals = np.log2(counts_filt[sample] + 1)
    raw_vals = raw_vals[raw_vals > 0]
    axes[0].hist(raw_vals, bins=60, density=True, alpha=0.35, color=c,
                 histtype="stepfilled", linewidth=0)

    tpm_vals = np.log2(tpm_filt[sample] + 1)
    tpm_vals = tpm_vals[tpm_vals > 0]
    axes[1].hist(tpm_vals, bins=60, density=True, alpha=0.35, color=c,
                 histtype="stepfilled", linewidth=0)

for ax, title, xlabel in zip(
    axes,
    ["Raw Counts (log₂ + 1)", "TPM Normalised (log₂ + 1)"],
    ["log₂(count + 1)",        "log₂(TPM + 1)"],
):
    ax.set_xlabel(xlabel); ax.set_ylabel("Density"); ax.set_title(title, fontweight="bold")
    ax.legend(handles=[
        mpatches.Patch(color=COND_COLORS["Control"],   alpha=0.6, label="Control"),
        mpatches.Patch(color=COND_COLORS["Treatment"], alpha=0.6, label="Treatment"),
    ], fontsize=9)

plt.tight_layout()
plt.savefig("results/figures/02_normalization.pdf", bbox_inches="tight")
plt.savefig("results/figures/02_normalization.png", bbox_inches="tight", dpi=300)
plt.show()
print("Saved → results/figures/02_normalization.pdf/.png")
"""))

# ── Section 7: PCA ────────────────────────────────────────────────────────────
cells.append(md("""
## 7. Principal Component Analysis (PCA)

PCA reveals global structure: replicates should cluster together and conditions
should separate along the first principal component if there is a transcriptional
response.
"""))

cells.append(code("""
# Prepare matrix: samples × genes, z-scored per gene
X_mat    = log2_tpm.T.values                          # (8, n_genes)
scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X_mat)

pca         = PCA(n_components=min(N_SAMPLES, 8))
pca_coords  = pca.fit_transform(X_scaled)
var_exp     = pca.explained_variance_ratio_ * 100

pca_df = pd.DataFrame(pca_coords[:, :4], index=SAMPLES,
                       columns=[f"PC{i+1}" for i in range(4)])
pca_df = pca_df.join(metadata[["condition", "batch"]])

print("Variance explained by top PCs:")
for i, v in enumerate(var_exp[:5]):
    print(f"  PC{i+1}: {v:.1f}%")
"""))

cells.append(code("""
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# ── PC1 vs PC2 coloured by condition ────────────────────────────────────────
for cond, grp in pca_df.groupby("condition"):
    axes[0].scatter(grp["PC1"], grp["PC2"],
                    c=COND_COLORS[cond], s=130, edgecolors="white",
                    linewidth=1.5, zorder=5, label=cond, alpha=0.9)
for idx, row in pca_df.iterrows():
    axes[0].annotate(idx, (row["PC1"], row["PC2"]),
                     fontsize=8, ha="left", xytext=(5, 4),
                     textcoords="offset points")
axes[0].set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
axes[0].set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
axes[0].set_title("PCA — Condition", fontweight="bold")
axes[0].legend(fontsize=10)
for ax in axes[:2]:
    ax.axhline(0, color="gray", ls="--", alpha=0.4, lw=0.8)
    ax.axvline(0, color="gray", ls="--", alpha=0.4, lw=0.8)

# ── PC1 vs PC2 coloured by batch ─────────────────────────────────────────────
batch_pal = {"A": "#2ca02c", "B": "#d62728"}
for batch, grp in pca_df.groupby("batch"):
    axes[1].scatter(grp["PC1"], grp["PC2"],
                    c=batch_pal[batch], s=130,
                    marker="o" if batch == "A" else "s",
                    edgecolors="white", linewidth=1.5,
                    zorder=5, label=f"Batch {batch}", alpha=0.9)
for idx, row in pca_df.iterrows():
    axes[1].annotate(idx, (row["PC1"], row["PC2"]),
                     fontsize=8, ha="left", xytext=(5, 4),
                     textcoords="offset points")
axes[1].set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
axes[1].set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
axes[1].set_title("PCA — Batch", fontweight="bold")
axes[1].legend(fontsize=10)

# ── Scree plot ────────────────────────────────────────────────────────────────
xs = range(1, len(var_exp) + 1)
axes[2].bar(xs, var_exp, color="steelblue", edgecolor="white")
axes[2].plot(xs, var_exp, "ko-", markersize=5)
axes[2].set_xlabel("Principal Component")
axes[2].set_ylabel("Variance Explained (%)")
axes[2].set_title("Scree Plot", fontweight="bold")
axes[2].set_xticks(list(xs))

plt.tight_layout()
plt.savefig("results/figures/03_pca.pdf", bbox_inches="tight")
plt.savefig("results/figures/03_pca.png", bbox_inches="tight", dpi=300)
plt.show()
print("Saved → results/figures/03_pca.pdf/.png")
"""))

# ── Section 8: DESeq2 ─────────────────────────────────────────────────────────
cells.append(md("""
## 8. Differential Expression — DESeq2

**pyDESeq2** is the official Python port of the R DESeq2 package. It uses:
- Negative-binomial GLM with shrunk dispersion estimates
- Wald test for contrasts
- Benjamini–Hochberg FDR correction
- Cook's distance filtering for outlier counts

> **Contrast:** Treatment vs Control
> Positive log₂FC = higher expression in Treatment.
"""))

cells.append(code("""
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds  import DeseqStats
    HAVE_PYDESEQ2 = True
    print("pyDESeq2 available — running DESeq2 pipeline.")
except ImportError:
    HAVE_PYDESEQ2 = False
    print("pyDESeq2 not found.  Install: pip install pydeseq2")
    print("Falling back to Welch t-test + BH correction.")
"""))

cells.append(code("""
if HAVE_PYDESEQ2:
    # pyDESeq2 expects counts as (samples × genes), integer dtype
    counts_deseq = counts_filt.T.astype(int)          # (8, n_filtered_genes)
    meta_deseq   = metadata[["condition"]].copy()
    counts_deseq = counts_deseq.loc[meta_deseq.index]

    print(f"Running DESeq2 on {counts_deseq.shape[1]:,} genes × {counts_deseq.shape[0]} samples …")

    dds = DeseqDataSet(
        counts        = counts_deseq,
        metadata      = meta_deseq,
        design_factors= "condition",
        refit_cooks   = True,
        n_cpus        = 4,
    )
    dds.deseq2()

    stat_res = DeseqStats(
        dds,
        contrast = ["condition", "Treatment", "Control"],
        alpha    = 0.05,
        cooks_filter         = True,
        independent_filter   = True,
    )
    stat_res.summary()

    de_results = stat_res.results_df.copy()
    de_results.index.name = "gene"
    de_results = de_results.reset_index()

else:
    # ── Fallback: Welch t-test + BH ──────────────────────────────────────────
    from statsmodels.stats.multitest import multipletests

    ctrl_s = [s for s in SAMPLES if metadata.loc[s,"condition"]=="Control"]
    trt_s  = [s for s in SAMPLES if metadata.loc[s,"condition"]=="Treatment"]
    lc     = log2_tpm

    rows = []
    for gene in lc.index:
        c, t = lc.loc[gene, ctrl_s].values, lc.loc[gene, trt_s].values
        lfc  = t.mean() - c.mean()
        tst, pv = stats.ttest_ind(t, c, equal_var=False)
        rows.append({"gene": gene, "baseMean": np.exp2(np.concatenate([c,t])).mean(),
                     "log2FoldChange": lfc, "stat": tst, "pvalue": pv})

    de_results = pd.DataFrame(rows)
    _, padj, *_ = multipletests(de_results["pvalue"].fillna(1), method="fdr_bh")
    de_results["padj"]  = padj
    de_results["lfcSE"] = np.abs(de_results["log2FoldChange"]) * 0.3

print(f"Total genes tested: {len(de_results):,}")
"""))

cells.append(code("""
# ── Annotate results ──────────────────────────────────────────────────────────
LFC_THR  = 1.0   # |log2FC| > 1 = 2-fold
PADJ_THR = 0.05  # FDR < 5 %

de_results["direction"] = "Not significant"
de_results.loc[
    (de_results.padj < PADJ_THR) & (de_results.log2FoldChange >  LFC_THR),
    "direction"] = "Up-regulated"
de_results.loc[
    (de_results.padj < PADJ_THR) & (de_results.log2FoldChange < -LFC_THR),
    "direction"] = "Down-regulated"

n_up   = (de_results.direction == "Up-regulated").sum()
n_down = (de_results.direction == "Down-regulated").sum()

print("=" * 50)
print("  DE Results (Treatment vs Control)")
print("=" * 50)
print(f"  Total tested  : {len(de_results):,}")
print(f"  Up-regulated  : {n_up:,}  (padj<{PADJ_THR}, log2FC>{LFC_THR})")
print(f"  Down-regulated: {n_down:,}  (padj<{PADJ_THR}, log2FC<-{LFC_THR})")
print("=" * 50)

# ── Save ──────────────────────────────────────────────────────────────────────
de_sorted = de_results.sort_values("padj")
de_sorted.to_csv("results/tables/DE_results_all.csv", index=False)
de_sorted[de_sorted.direction != "Not significant"] \
    .to_csv("results/tables/DE_results_significant.csv", index=False)

print("\nSaved → results/tables/DE_results_all.csv")
print("Saved → results/tables/DE_results_significant.csv")
display(de_sorted.head(10))
"""))

# Volcano
cells.append(md("### 8.1 Volcano Plot"))

cells.append(code("""
fig, ax = plt.subplots(figsize=(9, 7))

dir_pal  = {"Not significant": "#AAAAAA",
            "Up-regulated":    "#E74C3C",
            "Down-regulated":  "#3498DB"}
dir_size = {"Not significant": 12, "Up-regulated": 22, "Down-regulated": 22}
dir_zord = {"Not significant":  1, "Up-regulated":  5, "Down-regulated":  5}

for direction, df_sub in de_results.groupby("direction"):
    log_pv = -np.log10(df_sub["pvalue"].clip(lower=1e-300))
    ax.scatter(df_sub["log2FoldChange"], log_pv,
               c=dir_pal[direction], s=dir_size[direction],
               alpha=0.5 if direction == "Not significant" else 0.75,
               zorder=dir_zord[direction],
               label=f"{direction} (n={len(df_sub):,})")

# Label top 20 DEGs
top_label = de_results[de_results.direction != "Not significant"] \
               .nsmallest(20, "padj")
for _, row in top_label.iterrows():
    if pd.notna(row["pvalue"]) and row["pvalue"] > 0:
        ax.annotate(row["gene"],
                    (row["log2FoldChange"], -np.log10(row["pvalue"])),
                    fontsize=7.5, fontweight="bold", ha="center",
                    xytext=(0, 5), textcoords="offset points")

# Threshold lines
ax.axvline(-LFC_THR,  color="gray", ls="--", alpha=0.6, lw=1)
ax.axvline( LFC_THR,  color="gray", ls="--", alpha=0.6, lw=1)
ax.axhline(-np.log10(PADJ_THR), color="gray", ls="--", alpha=0.6, lw=1)

ax.set_xlabel("log₂(Fold Change)", fontsize=13)
ax.set_ylabel("−log₁₀(p-value)", fontsize=13)
ax.set_title("Volcano Plot: Treatment vs Control", fontsize=14, fontweight="bold")
ax.legend(fontsize=10, loc="upper left")

xlim = ax.get_xlim(); ylim = ax.get_ylim()
ax.text(xlim[1]*0.85, ylim[1]*0.96, f"↑ {n_up}", color="#E74C3C",
        fontsize=12, fontweight="bold", ha="center")
ax.text(xlim[0]*0.85, ylim[1]*0.96, f"↓ {n_down}", color="#3498DB",
        fontsize=12, fontweight="bold", ha="center")

plt.tight_layout()
plt.savefig("results/figures/04_volcano.pdf", bbox_inches="tight")
plt.savefig("results/figures/04_volcano.png", bbox_inches="tight", dpi=300)
plt.show()
print("Saved → results/figures/04_volcano.pdf/.png")
"""))

# ── Section 9: Heatmap ────────────────────────────────────────────────────────
cells.append(md("""
## 9. Heatmap of Top Differential Genes

Hierarchical clustering of log₂-TPM z-scores for the top 25 up- and 25
down-regulated genes reveals the expression patterns across samples.
"""))

cells.append(code("""
N_HM = 25   # top N per direction

top_up   = de_results[de_results.direction=="Up-regulated"  ].nsmallest(N_HM, "padj")["gene"].tolist()
top_down = de_results[de_results.direction=="Down-regulated"].nsmallest(N_HM, "padj")["gene"].tolist()
hm_genes = [g for g in top_up + top_down if g in log2_tpm.index]

hm_expr = log2_tpm.loc[hm_genes]

# Z-score across samples
hm_z = pd.DataFrame(
    zscore(hm_expr.values, axis=1),
    index=hm_expr.index, columns=hm_expr.columns,
)

col_colors = pd.Series(
    [COND_COLORS[metadata.loc[s, "condition"]] for s in hm_z.columns],
    index=hm_z.columns,
)
row_colors = pd.Series(
    ["#E74C3C" if g in top_up else "#3498DB" for g in hm_genes],
    index=hm_genes,
)

g = sns.clustermap(
    hm_z, col_colors=col_colors, row_colors=row_colors,
    cmap="RdBu_r", vmin=-3, vmax=3, center=0,
    figsize=(10, 14), linewidths=0,
    dendrogram_ratio=(0.1, 0.18),
    cbar_pos=(0.02, 0.82, 0.03, 0.12),
    xticklabels=True, yticklabels=True,
)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=9, rotation=45, ha="right")
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
g.cax.set_ylabel("Z-score", fontsize=9, rotation=90)
g.fig.suptitle(f"Top {len(hm_genes)} DEGs — Z-scored log₂-TPM",
               y=1.01, fontsize=13, fontweight="bold")

# Legend patches
leg_kw = dict(bbox_to_anchor=(1.15, -0.02), loc="upper right", fontsize=9)
g.ax_col_colors.legend(handles=[
    mpatches.Patch(color=COND_COLORS["Control"],   label="Control"),
    mpatches.Patch(color=COND_COLORS["Treatment"], label="Treatment"),
], title="Condition", **leg_kw)
g.ax_row_colors.legend(handles=[
    mpatches.Patch(color="#E74C3C", label="Up-regulated"),
    mpatches.Patch(color="#3498DB", label="Down-regulated"),
], title="Direction", loc="lower right", bbox_to_anchor=(0, -0.02), fontsize=9)

plt.savefig("results/figures/05_heatmap.pdf", bbox_inches="tight")
plt.savefig("results/figures/05_heatmap.png", bbox_inches="tight", dpi=300)
plt.show()
print("Saved → results/figures/05_heatmap.pdf/.png")
"""))

# ── Section 10: Boxplots ──────────────────────────────────────────────────────
cells.append(md("""
## 10. Boxplots of Top Individual Genes

Show the expression distribution of the most significantly changed genes,
with statistical significance annotations.
"""))

cells.append(code("""
# Top 5 up + 4 down = 9 panels in a 3×3 grid
box_up   = [g for g in de_results[de_results.direction=="Up-regulated"  ]
              .nsmallest(5,"padj")["gene"].tolist() if g in log2_tpm.index]
box_down = [g for g in de_results[de_results.direction=="Down-regulated"]
              .nsmallest(4,"padj")["gene"].tolist() if g in log2_tpm.index]
box_genes = (box_up + box_down)[:9]

fig, axes = plt.subplots(3, 3, figsize=(12, 10))
axes = axes.flatten()

for i, gene in enumerate(box_genes):
    ax = axes[i]
    df_box = pd.DataFrame({
        "Condition" : [metadata.loc[s, "condition"] for s in SAMPLES],
        "Expression": [log2_tpm.loc[gene, s] for s in SAMPLES],
    })
    sns.boxplot(data=df_box, x="Condition", y="Expression",
                palette=COND_COLORS, order=["Control","Treatment"],
                width=0.5, flierprops={"marker":"o","markersize":4}, ax=ax)
    sns.stripplot(data=df_box, x="Condition", y="Expression",
                  palette=COND_COLORS, order=["Control","Treatment"],
                  size=6, jitter=0.12, alpha=0.85,
                  edgecolor="white", linewidth=0.5, ax=ax)

    row = de_results[de_results.gene==gene]
    if not row.empty:
        lfc   = row["log2FoldChange"].values[0]
        padj  = row["padj"].values[0]
        dirn  = row["direction"].values[0]
        sig   = "***" if padj<0.001 else "**" if padj<0.01 else "*" if padj<0.05 else "ns"
        ymax  = df_box["Expression"].max()
        ax.text(0.5, 1.01, sig, transform=ax.transAxes,
                ha="center", fontsize=14 if sig!="ns" else 10, fontweight="bold")
        color = "#E74C3C" if dirn=="Up-regulated" else "#3498DB"
        ax.set_title(f"{gene}\nlog₂FC={lfc:+.2f}, FDR={padj:.1e}",
                     fontsize=9, fontweight="bold", color=color, pad=14)
    ax.set_xlabel("")
    ax.set_ylabel("log₂(TPM+1)" if i%3==0 else "", fontsize=9)

for i in range(len(box_genes), len(axes)):
    axes[i].set_visible(False)

plt.suptitle("Top Differentially Expressed Genes", fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig("results/figures/06_boxplots.pdf", bbox_inches="tight")
plt.savefig("results/figures/06_boxplots.png", bbox_inches="tight", dpi=300)
plt.show()
print("Saved → results/figures/06_boxplots.pdf/.png")
"""))

# ── Section 11: ORA ───────────────────────────────────────────────────────────
cells.append(md("""
## 11. Pathway Enrichment — Over-Representation Analysis (ORA)

Uses the **Enrichr** API (via gseapy) to ask: *are genes in my DEG list
over-represented in known pathways?*

**Tested databases:** KEGG 2021, GO Biological Process 2021, Reactome 2022
"""))

cells.append(code("""
try:
    import gseapy as gp
    HAVE_GSEAPY = True
    print(f"gseapy {gp.__version__} available.")
except ImportError:
    HAVE_GSEAPY = False
    print("gseapy not found. Install: pip install gseapy")
    print("Sections 11 and 12 require internet access.")
"""))

cells.append(code("""
GENE_SETS_ORA = ["KEGG_2021_Human", "GO_Biological_Process_2021", "Reactome_2022"]

# Use only real HGNC symbols (exclude ENSG placeholders)
def real_symbols(gene_list):
    return [g for g in gene_list if not g.startswith("ENSG")]

if HAVE_GSEAPY:
    up_genes   = real_symbols(de_results[de_results.direction=="Up-regulated"  ]["gene"].tolist())
    down_genes = real_symbols(de_results[de_results.direction=="Down-regulated"]["gene"].tolist())

    print(f"ORA input — Up: {len(up_genes)} genes | Down: {len(down_genes)} genes")
    print("Querying Enrichr (requires internet) …")

    ora_parts = []
    for direction, glist in [("Up-regulated", up_genes), ("Down-regulated", down_genes)]:
        if not glist:
            continue
        try:
            enr = gp.enrichr(gene_list=glist, gene_sets=GENE_SETS_ORA,
                              organism="Human", outdir=None, no_plot=True, cutoff=0.1)
            df  = enr.results.copy()
            df["Direction"] = direction
            ora_parts.append(df)
        except Exception as e:
            print(f"  {direction} ORA failed: {e}")

    if ora_parts:
        ora_df = pd.concat(ora_parts, ignore_index=True)
        ora_df = ora_df[ora_df["Adjusted P-value"] < 0.05]
        ora_df.to_csv("results/tables/ORA_pathway_enrichment.csv", index=False)
        print(f"Significant terms (FDR<0.05): {len(ora_df)}")
    else:
        print("No ORA results.")
        ora_df = pd.DataFrame()
"""))

cells.append(code("""
# ── ORA dot plot ──────────────────────────────────────────────────────────────
if HAVE_GSEAPY and "ora_df" in dir() and len(ora_df) > 0:

    top_terms = []
    for db in GENE_SETS_ORA:
        for dirn in ["Up-regulated", "Down-regulated"]:
            sub = ora_df[(ora_df.Gene_set==db) & (ora_df.Direction==dirn)].head(6)
            top_terms.append(sub)
    plot_df = pd.concat(top_terms).drop_duplicates("Term").reset_index(drop=True)
    plot_df["-log10FDR"] = -np.log10(plot_df["Adjusted P-value"].clip(1e-50))
    plot_df["GeneRatio"] = plot_df["Overlap"].apply(
        lambda x: int(x.split("/")[0]) / max(int(x.split("/")[1]), 1) if "/" in str(x) else 0)
    plot_df["Term_short"] = plot_df["Term"].str[:55]
    plot_df = plot_df.sort_values(["Direction", "-log10FDR"], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(11, max(6, len(plot_df)*0.38)))
    ax.scatter(
        plot_df["-log10FDR"],
        range(len(plot_df)),
        c=plot_df["Direction"].map({"Up-regulated":"#E74C3C","Down-regulated":"#3498DB"}),
        s=plot_df["GeneRatio"]*600,
        alpha=0.85, edgecolors="white", linewidth=0.5,
    )
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["Term_short"], fontsize=8)
    ax.axvline(-np.log10(0.05), color="gray", ls="--", alpha=0.6, lw=1, label="FDR=0.05")
    ax.set_xlabel("−log₁₀(Adjusted P-value)", fontsize=11)
    ax.set_title("ORA: Pathway Enrichment", fontsize=13, fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    ax.legend(handles=[
        mpatches.Patch(color="#E74C3C", label="Up-regulated genes"),
        mpatches.Patch(color="#3498DB", label="Down-regulated genes"),
    ], fontsize=10)

    plt.tight_layout()
    plt.savefig("results/figures/07_ORA_dotplot.pdf", bbox_inches="tight")
    plt.savefig("results/figures/07_ORA_dotplot.png", bbox_inches="tight", dpi=300)
    plt.show()
    print("Saved → results/figures/07_ORA_dotplot.pdf/.png")
"""))

# ── Section 12: GSEA ─────────────────────────────────────────────────────────
cells.append(md("""
## 12. Gene Set Enrichment Analysis (GSEA)

Unlike ORA, GSEA uses **all** tested genes ranked by a continuous metric
(`−log₁₀(p-value) × sign(log₂FC)`), so it is not limited to a binary DEG cutoff.

> Use `permutation_num=1000` for publication-quality results.
> We use 100 here for speed.
"""))

cells.append(code("""
if HAVE_GSEAPY:
    # Build ranked list (all genes, real symbols preferred at top)
    de_sym = de_results[~de_results.gene.str.startswith("ENSG")].copy()
    de_sym["rank_metric"] = (
        -np.log10(de_sym["pvalue"].clip(lower=1e-300)) *
        np.sign(de_sym["log2FoldChange"])
    )
    ranking = (de_sym.dropna(subset=["rank_metric"])
                     .sort_values("rank_metric", ascending=False)
                     .set_index("gene")["rank_metric"])

    print(f"Ranked gene list: {len(ranking):,} genes")
    print(f"  Max rank score : {ranking.max():.1f}")
    print(f"  Min rank score : {ranking.min():.1f}")

    print("Running GSEA prerank (requires internet) …")
    try:
        gsea_res = gp.prerank(
            rnk            = ranking,
            gene_sets      = ["KEGG_2021_Human", "Hallmark_2020"],
            outdir         = None,
            min_size       = 15,
            max_size       = 500,
            permutation_num= 100,
            seed           = 42,
            no_plot        = True,
            verbose        = False,
        )
        gsea_df = gsea_res.res2d.copy()
        gsea_df = gsea_df[gsea_df["FDR q-val"] < 0.25].sort_values("NES", ascending=False)
        gsea_df.to_csv("results/tables/GSEA_results.csv", index=False)
        print(f"Significant gene sets (FDR<0.25): {len(gsea_df)}")
        display(gsea_df.head(10))
    except Exception as e:
        print(f"GSEA failed: {e}")
        gsea_df = pd.DataFrame()
"""))

cells.append(code("""
# ── GSEA bar plot (NES) ───────────────────────────────────────────────────────
if HAVE_GSEAPY and "gsea_df" in dir() and len(gsea_df) > 0:

    top_gsea = pd.concat([
        gsea_df[gsea_df.NES > 0].head(10),
        gsea_df[gsea_df.NES < 0].tail(10),
    ]).reset_index(drop=True)
    top_gsea["Term_short"] = top_gsea["Term"].str[:60]
    top_gsea = top_gsea.sort_values("NES")

    fig, ax = plt.subplots(figsize=(10, max(5, len(top_gsea)*0.38)))
    colors = ["#3498DB" if n < 0 else "#E74C3C" for n in top_gsea["NES"]]
    ax.barh(top_gsea["Term_short"], top_gsea["NES"],
            color=colors, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Normalised Enrichment Score (NES)", fontsize=11)
    ax.set_title("GSEA: Gene Set Enrichment Analysis", fontsize=13, fontweight="bold")
    ax.legend(handles=[
        mpatches.Patch(color="#E74C3C", label="Enriched in Treatment"),
        mpatches.Patch(color="#3498DB", label="Enriched in Control"),
    ], fontsize=10)

    plt.tight_layout()
    plt.savefig("results/figures/08_GSEA_barplot.pdf", bbox_inches="tight")
    plt.savefig("results/figures/08_GSEA_barplot.png", bbox_inches="tight", dpi=300)
    plt.show()
    print("Saved → results/figures/08_GSEA_barplot.pdf/.png")
"""))

# ── Section 13: Gene–gene interaction ────────────────────────────────────────
cells.append(md("""
## 13. Gene–Gene Interaction Network (STRING)

Query the **STRING** protein interaction database for all significant DEGs,
then visualise the interaction network with genes coloured by expression direction.
Hub genes (high degree centrality) are candidate master regulators.
"""))

cells.append(code("""
def query_string(gene_list, species=9606, min_score=400, max_genes=100):
    \"\"\"
    Retrieve protein–protein interactions from STRING v12.

    Parameters
    ----------
    gene_list  : list of HGNC gene symbols
    species    : NCBI taxonomy ID (9606 = Homo sapiens)
    min_score  : combined interaction score threshold (0–1000)
    max_genes  : cap on query size

    Returns
    -------
    pd.DataFrame with columns gene1, gene2, score (0–1)
    \"\"\"
    url    = "https://string-db.org/api/json/network"
    params = {
        "identifiers"   : "%0d".join(gene_list[:max_genes]),
        "species"       : species,
        "required_score": min_score,
        "caller_identity": "rnaseq_tutorial",
    }
    try:
        r = requests.post(url, data=params, timeout=30)
        r.raise_for_status()
        data = r.json()
        if not data:
            return pd.DataFrame()
        df = pd.DataFrame(data)[["preferredName_A", "preferredName_B", "score"]]
        df.columns = ["gene1", "gene2", "score"]
        df["score"] = df["score"].astype(float) / 1000     # normalise to 0–1
        df = df[df.gene1 != df.gene2]
        return df
    except Exception as e:
        print(f"STRING query failed: {e}")
        return pd.DataFrame()
"""))

cells.append(code("""
sig_genes_real = real_symbols(
    de_results[de_results.direction != "Not significant"]["gene"].tolist())

print(f"Querying STRING for {len(sig_genes_real)} DEGs (requires internet) …")
interactions = query_string(sig_genes_real, min_score=400)

if interactions.empty:
    print("No interactions returned — using curated fallback network.")
    interactions = pd.DataFrame({
        "gene1": ["CDK1","CDK1","CDK2","CCNB1","PLK1","AURKB","TOP2A",
                  "TP53","CDKN1A","BAX","CASP3","CASP9","STAT1","IRF1",
                  "TNF","NFKB1","IL6","MYC","E2F1","BCL2"],
        "gene2": ["CCNB1","PCNA","CCNE1","CDK1","CDK1","CDK1","BIRC5",
                  "CDKN1A","CDK2","BCL2","CASP9","CYCS","IRF1","ISG15",
                  "NFKB1","RELA","STAT3","CDK4","CCND1","MCL1"],
        "score": [0.97,0.89,0.94,0.97,0.92,0.88,0.87,
                  0.98,0.96,0.92,0.95,0.91,0.93,0.90,
                  0.86,0.94,0.89,0.88,0.91,0.87],
    })
    print(f"Fallback: {len(interactions)} interactions loaded.")
else:
    print(f"Retrieved {len(interactions)} interactions among "
          f"{len(set(interactions.gene1)|set(interactions.gene2))} genes.")
"""))

cells.append(code("""
# ── Build networkx graph ──────────────────────────────────────────────────────
G = nx.Graph()
for _, row in interactions.iterrows():
    G.add_edge(row.gene1, row.gene2, weight=row.score)

# Attach DE metadata to nodes
dir_map   = de_results.set_index("gene")["direction"].to_dict()
lfc_map   = de_results.set_index("gene")["log2FoldChange"].to_dict()

for node in G.nodes():
    G.nodes[node]["direction"] = dir_map.get(node, "Not significant")
    G.nodes[node]["lfc"]       = lfc_map.get(node, 0.0)

# Restrict to largest connected component
largest_cc = max(nx.connected_components(G), key=len)
G_main     = G.subgraph(largest_cc).copy()

deg_cent = nx.degree_centrality(G_main)
btw_cent = nx.betweenness_centrality(G_main)

print(f"Network  — nodes: {G.number_of_nodes()}, edges: {G.number_of_edges()}")
print(f"Largest CC — nodes: {G_main.number_of_nodes()}, edges: {G_main.number_of_edges()}")
print()
print("Top 10 hub genes (degree centrality):")
for gene, cent in sorted(deg_cent.items(), key=lambda x: -x[1])[:10]:
    dirn = G_main.nodes[gene]["direction"]
    print(f"  {gene:<12} {cent:.3f}  [{dirn}]")
"""))

cells.append(code("""
node_pal = {"Up-regulated": "#E74C3C", "Down-regulated": "#3498DB",
            "Not significant": "#AAAAAA"}

fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# ── Network layout ────────────────────────────────────────────────────────────
pos   = nx.spring_layout(G_main,
                         k=2.5/np.sqrt(G_main.number_of_nodes()),
                         seed=42, iterations=120)
ncolors = [node_pal[G_main.nodes[n]["direction"]] for n in G_main.nodes()]
nsizes  = [350 + deg_cent[n]*4000 for n in G_main.nodes()]
ewidths = [G_main[u][v]["weight"]*3 for u, v in G_main.edges()]

nx.draw_networkx_edges(G_main, pos, ax=axes[0], alpha=0.35,
                       width=ewidths, edge_color="gray")
nx.draw_networkx_nodes(G_main, pos, ax=axes[0],
                       node_color=ncolors, node_size=nsizes,
                       alpha=0.9, edgecolors="white", linewidths=1)
# Label high-degree nodes
top50pct = np.percentile(list(deg_cent.values()), 50)
nx.draw_networkx_labels(G_main, pos,
    labels={n: n for n in G_main.nodes() if deg_cent[n] >= top50pct},
    font_size=7.5, font_weight="bold", ax=axes[0])

axes[0].set_title("Gene Interaction Network (STRING)", fontsize=13, fontweight="bold")
axes[0].axis("off")
axes[0].legend(handles=[
    mpatches.Patch(color=c, label=d)
    for d, c in node_pal.items()
], fontsize=10, loc="lower left")

# ── Top hub genes bar chart ───────────────────────────────────────────────────
hub_df = (pd.DataFrame({"Gene": list(deg_cent.keys()),
                         "Degree Centrality": list(deg_cent.values())})
            .sort_values("Degree Centrality", ascending=False)
            .head(15)
            .reset_index(drop=True))
hub_df["Direction"] = hub_df["Gene"].map(
    lambda g: G_main.nodes[g]["direction"] if g in G_main.nodes() else "Not significant")
hub_colors = hub_df["Direction"].map(node_pal).values

axes[1].barh(hub_df["Gene"][::-1], hub_df["Degree Centrality"][::-1],
             color=hub_colors[::-1], edgecolor="white", linewidth=0.5)
axes[1].set_xlabel("Degree Centrality", fontsize=11)
axes[1].set_title("Top Hub Genes", fontsize=13, fontweight="bold")

plt.tight_layout()
plt.savefig("results/figures/09_gene_network.pdf", bbox_inches="tight")
plt.savefig("results/figures/09_gene_network.png", bbox_inches="tight", dpi=300)
plt.show()
print("Saved → results/figures/09_gene_network.pdf/.png")
"""))

# ── Section 14: TF binding prediction ────────────────────────────────────────
cells.append(md("""
## 14. Transcription Factor Binding Prediction

Two complementary approaches:

| Method | Tool | What it asks |
|--------|------|-------------|
| **TF Enrichment** | ChEA3 API | Which TFs have targets over-represented among my DEGs? |
| **TF Activity** | decoupler-py + CollecTRI | Which TFs show differential *activity* across conditions? |

### 14.1 ChEA3 — TF–Target Enrichment

ChEA3 integrates ENCODE ChIP-seq, literature ChIP-seq, ARCHS4 co-expression,
and other resources to rank TFs whose known targets are enriched in a gene list.
"""))

cells.append(code("""
CHEA3_URL = "https://maayanlab.cloud/chea3/api/enrich/"

def query_chea3(gene_list: list, query_name: str = "DEG_query") -> dict:
    \"\"\"
    Submit a gene list to the ChEA3 API.

    Returns
    -------
    dict  keyed by library name; each value is a list of TF dicts with
          fields: TF, Rank, Score, Overlapping_Genes, ...
    \"\"\"
    payload = {"query_name": query_name, "gene_set": gene_list}
    try:
        r = requests.post(CHEA3_URL, json=payload, timeout=60)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        print(f"ChEA3 query failed: {e}")
        return {}


def parse_chea3(result: dict, library: str = "Integrated--meanRank",
                top_n: int = 25) -> pd.DataFrame:
    \"\"\"
    Parse ChEA3 JSON into a tidy DataFrame.
    Tries 'Integrated--meanRank' first, falls back to the first available library.
    \"\"\"
    if not result:
        return pd.DataFrame()
    if library not in result:
        library = next(iter(result))
        print(f"  Using library: {library}")
    rows = result[library][:top_n]
    df   = pd.DataFrame(rows)
    df["Rank"]  = pd.to_numeric(df["Rank"],  errors="coerce")
    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    return df.sort_values("Rank").reset_index(drop=True)
"""))

cells.append(code("""
up_real   = real_symbols(de_results[de_results.direction=="Up-regulated"  ]["gene"].tolist())
down_real = real_symbols(de_results[de_results.direction=="Down-regulated"]["gene"].tolist())

print(f"Querying ChEA3 with {len(up_real)} up and {len(down_real)} down DEGs …")

chea3_up   = query_chea3(up_real,   "Upregulated_DEGs")
chea3_down = query_chea3(down_real, "Downregulated_DEGs")

tf_up   = parse_chea3(chea3_up,   top_n=20)
tf_down = parse_chea3(chea3_down, top_n=20)

if not tf_up.empty:
    tf_up["Direction"]   = "Up-regulated targets"
if not tf_down.empty:
    tf_down["Direction"] = "Down-regulated targets"

tf_all = pd.concat([tf_up, tf_down], ignore_index=True)
if not tf_all.empty:
    tf_all.to_csv("results/tables/ChEA3_TF_enrichment.csv", index=False)
    print("Saved → results/tables/ChEA3_TF_enrichment.csv")
    display(tf_all.head(10))
"""))

cells.append(code("""
# ── ChEA3 bar chart ───────────────────────────────────────────────────────────
if not tf_up.empty or not tf_down.empty:

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    for ax, df_tf, title, color in zip(
        axes,
        [tf_up.head(15), tf_down.head(15)],
        ["TFs enriched in Up-regulated DEGs",
         "TFs enriched in Down-regulated DEGs"],
        ["#E74C3C", "#3498DB"],
    ):
        if df_tf.empty:
            ax.text(0.5, 0.5, "No results", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(title, fontweight="bold")
            continue

        df_plot = df_tf.sort_values("Score", ascending=True)
        ax.barh(df_plot["TF"], df_plot["Score"],
                color=color, edgecolor="white", linewidth=0.4, alpha=0.85)
        ax.set_xlabel("ChEA3 Score (higher = more enriched)", fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.tick_params(labelsize=9)

    plt.suptitle("ChEA3: Predicted Upstream Transcription Factors",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig("results/figures/10a_ChEA3_TF_enrichment.pdf", bbox_inches="tight")
    plt.savefig("results/figures/10a_ChEA3_TF_enrichment.png", bbox_inches="tight", dpi=300)
    plt.show()
    print("Saved → results/figures/10a_ChEA3_TF_enrichment.pdf/.png")
"""))

cells.append(md("""
### 14.2 TF Activity Inference — decoupler-py + CollecTRI

**CollecTRI** is a manually curated TF–target regulon.
**decoupler-py** uses a univariate linear model (ULM) to score how much each TF's
target programme is up- or down-regulated in each sample, then compares activities
across conditions.
"""))

cells.append(code("""
try:
    import decoupler as dc
    HAVE_DC = True
    print(f"decoupler {dc.__version__} available.")
except ImportError:
    HAVE_DC = False
    print("decoupler not found. Install: pip install decoupler")
    print("Section 14.2 will be skipped.")
"""))

cells.append(code("""
if HAVE_DC:
    print("Downloading CollecTRI regulon …")
    try:
        regulon = dc.get_collectri(organism="human", split_complexes=False)
        print(f"  {len(regulon):,} TF–target interactions, "
              f"{regulon['source'].nunique()} TFs")

        # log2-TPM matrix: genes × samples
        mat = np.log2(tpm_filt + 1)

        print("Inferring TF activities with ULM …")
        acts, pvals = dc.run_ulm(
            mat     = mat.T,         # samples × genes
            net     = regulon,
            source  = "source",
            target  = "target",
            weight  = "weight",
            verbose = False,
        )

        # Differential TF activity: t-test Control vs Treatment
        ctrl_s = [s for s in SAMPLES if metadata.loc[s,"condition"]=="Control"]
        trt_s  = [s for s in SAMPLES if metadata.loc[s,"condition"]=="Treatment"]

        tf_stats = []
        for tf in acts.columns:
            c_vals = acts.loc[ctrl_s, tf].values
            t_vals = acts.loc[trt_s,  tf].values
            t_stat, pv = stats.ttest_ind(t_vals, c_vals, equal_var=False)
            mean_diff  = t_vals.mean() - c_vals.mean()
            tf_stats.append({"TF": tf, "mean_diff": mean_diff,
                              "t_stat": t_stat, "pvalue": pv})

        tf_act_df = pd.DataFrame(tf_stats)
        from statsmodels.stats.multitest import multipletests
        _, tf_act_df["padj"], *_ = multipletests(
            tf_act_df["pvalue"].fillna(1), method="fdr_bh")
        tf_act_df = tf_act_df.sort_values("padj")
        tf_act_df.to_csv("results/tables/TF_activity_decoupler.csv", index=False)

        print(f"Saved → results/tables/TF_activity_decoupler.csv")
        sig_tfs = tf_act_df[tf_act_df.padj < 0.05]
        print(f"Significantly differential TFs (FDR<0.05): {len(sig_tfs)}")
        display(tf_act_df.head(15))

    except Exception as e:
        print(f"decoupler analysis failed: {e}")
        tf_act_df = pd.DataFrame()
"""))

cells.append(code("""
# ── TF activity heatmap ───────────────────────────────────────────────────────
if HAVE_DC and "acts" in dir() and len(acts) > 0 and "tf_act_df" in dir() \
        and len(tf_act_df) > 0:

    top_tf_act = tf_act_df.nsmallest(30, "padj")["TF"].tolist()
    top_tf_act = [t for t in top_tf_act if t in acts.columns]

    act_hm = acts[top_tf_act].T   # TFs × samples

    col_colors = pd.Series(
        [COND_COLORS[metadata.loc[s, "condition"]] for s in act_hm.columns],
        index=act_hm.columns)

    act_hm_z = pd.DataFrame(
        zscore(act_hm.values, axis=1),
        index=act_hm.index, columns=act_hm.columns)

    g = sns.clustermap(
        act_hm_z, col_colors=col_colors,
        cmap="RdBu_r", vmin=-2, vmax=2, center=0,
        figsize=(9, 10), linewidths=0, yticklabels=True,
        dendrogram_ratio=(0.1, 0.15),
        cbar_pos=(0.02, 0.82, 0.03, 0.12),
    )
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=9,
                                  rotation=45, ha="right")
    g.cax.set_ylabel("Z-score", fontsize=9)
    g.fig.suptitle("TF Activity (decoupler ULM + CollecTRI)\nZ-scored across samples",
                   y=1.01, fontsize=12, fontweight="bold")
    g.ax_col_colors.legend(handles=[
        mpatches.Patch(color=COND_COLORS["Control"],   label="Control"),
        mpatches.Patch(color=COND_COLORS["Treatment"], label="Treatment"),
    ], loc="upper right", bbox_to_anchor=(1.4, 0), fontsize=9)

    plt.savefig("results/figures/10b_TF_activity_heatmap.pdf", bbox_inches="tight")
    plt.savefig("results/figures/10b_TF_activity_heatmap.png", bbox_inches="tight", dpi=300)
    plt.show()
    print("Saved → results/figures/10b_TF_activity_heatmap.pdf/.png")
"""))

# ── Summary ───────────────────────────────────────────────────────────────────
cells.append(md("""
## 15. Summary

### Pipeline complete ✓

| Output file | Contents |
|------------|---------|
| `results/tables/DE_results_all.csv` | All tested genes with log₂FC, p-value, FDR |
| `results/tables/DE_results_significant.csv` | Significant DEGs only |
| `results/tables/ORA_pathway_enrichment.csv` | ORA enrichment results |
| `results/tables/GSEA_results.csv` | GSEA prerank results |
| `results/tables/ChEA3_TF_enrichment.csv` | TF enrichment (ChEA3) |
| `results/tables/TF_activity_decoupler.csv` | Differential TF activity |
| `results/figures/01_qc_summary` | Library sizes, gene detection, sample correlation |
| `results/figures/02_normalization` | Raw vs TPM count distributions |
| `results/figures/03_pca` | PCA plots coloured by condition and batch |
| `results/figures/04_volcano` | Volcano plot |
| `results/figures/05_heatmap` | Clustered heatmap of top DEGs |
| `results/figures/06_boxplots` | Boxplots of top individual genes |
| `results/figures/07_ORA_dotplot` | ORA pathway dotplot |
| `results/figures/08_GSEA_barplot` | GSEA NES barplot |
| `results/figures/09_gene_network` | STRING gene interaction network |
| `results/figures/10a_ChEA3_TF_enrichment` | ChEA3 TF bar charts |
| `results/figures/10b_TF_activity_heatmap` | TF activity heatmap |

### Next steps
- Validate top DEGs with qPCR or orthogonal data
- Overlay TF binding sites with ATAC-seq / ChIP-seq peaks
- Run multi-contrast or time-course DE analysis
- Integrate with proteomics or single-cell RNA-seq data
"""))

cells.append(code("""
print("=" * 60)
print("  OUTPUT SUMMARY")
print("=" * 60)
up_n   = (de_results.direction=="Up-regulated").sum()
down_n = (de_results.direction=="Down-regulated").sum()
print(f"  Samples          : {N_SAMPLES} ({', '.join(metadata.condition.value_counts().to_dict().keys())})")
print(f"  Genes analysed   : {counts_filt.shape[0]:,}")
print(f"  Up-regulated DEGs: {up_n}")
print(f"  Down-regulated   : {down_n}")
print("=" * 60)
print()
print("Results saved to:")
for f in sorted(
    [os.path.join(r, f) for r, _, fs in os.walk("results") for f in fs
     if f.endswith((".csv",".pdf",".png"))]
):
    print(f"  {f}")
"""))

# =============================================================================
# Assemble notebook JSON
# =============================================================================
notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {
            "name": "python",
            "version": "3.9.0",
        },
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out_path = os.path.join(os.path.dirname(__file__), "RNAseq_pipeline_tutorial.ipynb")
with open(out_path, "w") as fh:
    json.dump(notebook, fh, indent=1, ensure_ascii=False)

print(f"Notebook written: {out_path}")
print(f"  Cells: {len(cells)}  "
      f"({sum(1 for c in cells if c['cell_type']=='code')} code, "
      f"{sum(1 for c in cells if c['cell_type']=='markdown')} markdown)")
