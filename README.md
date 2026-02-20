# RNA-seq Data Analysis Pipeline

A complete, reproducible RNA-seq analysis pipeline from raw FASTQ files to biological insights.

## Pipeline Overview

```
FASTQ files
    │
    ▼
FastQC / MultiQC        ── raw read quality control
    │
    ▼
STAR (2-pass)           ── splice-aware genome alignment
    │
    ▼
featureCounts           ── gene-level read quantification
    │
    ▼
TPM Normalization       ── length + depth correction
    │
    ▼
DESeq2 (pyDESeq2)       ── differential expression testing
    │
    ├── Volcano plot
    ├── Heatmap of top DEGs
    └── Boxplots of individual genes
    │
    ▼
Pathway Analysis
    ├── ORA (Enrichr / gseapy)     ── over-representation analysis
    └── GSEA (prerank / gseapy)    ── gene set enrichment analysis
    │
    ▼
Gene–Gene Interaction Network
    └── STRING API + networkx      ── hub gene identification
    │
    ▼
TF Binding Prediction
    ├── ChEA3 API                  ── TF–target enrichment
    └── decoupler-py + CollecTRI   ── TF activity inference
```

**Input:**  `data/raw_counts.csv`, `data/gene_lengths.csv`, `data/metadata.csv`
**Output:** `results/figures/` (PDF + PNG), `results/tables/` (CSV)

---

## Repository Structure

```
RNAanalysis/
├── notebooks/
│   └── RNAseq_pipeline_tutorial.ipynb   ← main tutorial (start here)
├── scripts/
│   └── upstream_processing.sh           ← shell pipeline (FastQC→STAR→featureCounts)
├── data/                                ← auto-created by notebook
├── results/                             ← auto-created by notebook
└── README.md
```

---

## Quick Start

### Option A — Run the full tutorial notebook (simulated data included)

```bash
pip install jupyter pydeseq2 gseapy networkx decoupler \
            matplotlib seaborn pandas numpy scipy scikit-learn \
            adjustText requests statsmodels

jupyter notebook notebooks/RNAseq_pipeline_tutorial.ipynb
```

Run all cells. Simulated data is generated in **Section 3** so no real FASTQ
files are needed. Replace the simulated files with your own to analyse real data.

### Option B — Run the upstream shell pipeline on real data

```bash
# Install bioinformatics tools (conda recommended)
conda install -c bioconda fastqc star subread samtools
pip install multiqc

# Edit SAMPLES list and paths in the script, then:
bash scripts/upstream_processing.sh
```

Output: `data/raw_counts.csv` and `data/gene_lengths.csv`
Then open the notebook and skip Section 3.

---

## Notebook Sections

| # | Section | Tool(s) |
|---|---------|---------|
| 0 | Setup & Dependencies | pip |
| 1 | Sample Metadata | pandas |
| 2 | Upstream Processing (shell) | FastQC, STAR, featureCounts |
| 3 | Simulate Tutorial Dataset | numpy (NB distribution) |
| 4 | Load Count Matrix | pandas |
| 5 | Count-Level QC | matplotlib, seaborn |
| 6 | TPM Normalization | custom + pandas |
| 7 | PCA Analysis | scikit-learn |
| 8 | Differential Expression | pyDESeq2 |
| 9 | Heatmap of Top DEGs | seaborn clustermap |
| 10 | Boxplots of Top Genes | seaborn |
| 11 | ORA Pathway Enrichment | gseapy (Enrichr) |
| 12 | GSEA | gseapy (prerank) |
| 13 | Gene–Gene Interaction Network | STRING API, networkx |
| 14 | TF Binding Prediction | ChEA3 API, decoupler-py |
| 15 | Summary | — |

---

## Output Files

| File | Description |
|------|-------------|
| `results/tables/DE_results_all.csv` | All tested genes: log₂FC, p-value, FDR |
| `results/tables/DE_results_significant.csv` | Significant DEGs only |
| `results/tables/ORA_pathway_enrichment.csv` | ORA (Enrichr) results |
| `results/tables/GSEA_results.csv` | GSEA prerank results |
| `results/tables/ChEA3_TF_enrichment.csv` | TF enrichment (ChEA3) |
| `results/tables/TF_activity_decoupler.csv` | Differential TF activity |
| `results/figures/01_qc_summary` | Library sizes, detection, correlation |
| `results/figures/02_normalization` | Raw vs TPM distributions |
| `results/figures/03_pca` | PCA (condition + batch) + scree |
| `results/figures/04_volcano` | Volcano plot |
| `results/figures/05_heatmap` | Hierarchically clustered heatmap |
| `results/figures/06_boxplots` | Top gene expression boxplots |
| `results/figures/07_ORA_dotplot` | ORA dotplot |
| `results/figures/08_GSEA_barplot` | GSEA NES barplot |
| `results/figures/09_gene_network` | STRING interaction network + hub genes |
| `results/figures/10a_ChEA3_TF_enrichment` | ChEA3 TF bar charts |
| `results/figures/10b_TF_activity_heatmap` | TF activity heatmap |

All figures saved as both **PDF** (vector, publication-ready) and **PNG** (300 dpi).

---

## Requirements

| Package | Version | Purpose |
|---------|---------|---------|
| pydeseq2 | ≥ 0.4 | Differential expression |
| gseapy | ≥ 1.0 | ORA + GSEA |
| decoupler | ≥ 1.4 | TF activity inference |
| networkx | ≥ 3.0 | Network analysis |
| scikit-learn | ≥ 1.0 | PCA |
| seaborn | ≥ 0.12 | Visualisation |
| pandas | ≥ 1.5 | Data handling |
| numpy | ≥ 1.23 | Numerics |
| requests | any | STRING + ChEA3 APIs |

**Bioinformatics tools** (for upstream shell script only):

| Tool | Purpose |
|------|---------|
| FastQC | Per-read quality metrics |
| MultiQC | Aggregate QC reports |
| STAR ≥ 2.7 | RNA-seq alignment |
| featureCounts (subread) | Gene-level counting |
| samtools | BAM indexing / flagstat |

---

## Notes on Key Parameters

- **Strandedness (`-s` in featureCounts):** use `infer_experiment.py` (RSeQC) to determine; most dUTP / TruSeq kits are reverse-stranded (`-s 2`).
- **DESeq2 contrast:** `["condition", "Treatment", "Control"]` — positive log₂FC = higher in Treatment.
- **TPM vs raw counts:** TPM is used for visualisation; DESeq2 always operates on raw integer counts with its own internal normalisation.
- **GSEA permutations:** set `permutation_num=1000` for publication results (100 used here for speed).

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **STAR:** Dobin A et al. *Bioinformatics* 2013
- **featureCounts:** Liao Y et al. *Bioinformatics* 2014
- **DESeq2:** Love MI et al. *Genome Biology* 2014
- **gseapy / Enrichr:** Chen EY et al. *BMC Bioinformatics* 2013
- **STRING:** Szklarczyk D et al. *Nucleic Acids Res* 2023
- **ChEA3:** Keenan AB et al. *Nucleic Acids Res* 2019
- **decoupler:** Badia-i-Montpas P et al. *Bioinformatics* 2022
- **CollecTRI:** Müller-Dott S et al. *Nucleic Acids Res* 2023
