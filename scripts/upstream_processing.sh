#!/usr/bin/env bash
# =============================================================================
#  RNA-seq Upstream Processing Pipeline
#  FastQC / MultiQC  →  STAR  →  featureCounts
# =============================================================================
#
#  Tool installation (conda recommended):
#    conda install -c bioconda fastqc multiqc star subread samtools
#    pip install multiqc
#
#  Expected input layout:
#    data/fastq/{sample}_R1.fastq.gz   (R2 for paired-end)
#    genome/genome.fa                  (reference FASTA, e.g. GRCh38)
#    genome/annotation.gtf             (Ensembl / GENCODE GTF)
#
#  Usage:
#    bash scripts/upstream_processing.sh
# =============================================================================
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
FASTQ_DIR="data/fastq"
GENOME_FA="genome/genome.fa"
GTF="genome/annotation.gtf"
STAR_INDEX="genome/STAR_index"
RESULTS="results"
THREADS=8

# Adjust strandedness for featureCounts:
#   0 = unstranded  |  1 = stranded  |  2 = reverse-stranded (most dUTP kits)
STRANDEDNESS=2

SAMPLES=(ctrl_1 ctrl_2 ctrl_3 ctrl_4 trt_1 trt_2 trt_3 trt_4)

# Create directory tree
mkdir -p "${RESULTS}"/{qc/fastqc_raw,qc/fastqc_trimmed,trimmed,aligned,counts}

# =============================================================================
#  STEP 0  –  Build STAR genome index  (run once, ~30 min, needs ≥32 GB RAM)
# =============================================================================
build_star_index() {
    echo ">>> Building STAR genome index …"
    mkdir -p "${STAR_INDEX}"
    STAR \
        --runMode            genomeGenerate \
        --runThreadN         "${THREADS}" \
        --genomeDir          "${STAR_INDEX}" \
        --genomeFastaFiles   "${GENOME_FA}" \
        --sjdbGTFfile        "${GTF}" \
        --sjdbOverhang       149          # read_length − 1  (adjust!)
    echo ">>> Index built: ${STAR_INDEX}"
}
# Uncomment to build index on first run:
# build_star_index

# =============================================================================
#  STEP 1  –  FastQC on raw reads
# =============================================================================
echo ">>> [1/5] FastQC – raw reads"
for SAMPLE in "${SAMPLES[@]}"; do
    fastqc \
        "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
        "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" \
        --outdir "${RESULTS}/qc/fastqc_raw" \
        --threads "${THREADS}" \
        --quiet
done

multiqc \
    "${RESULTS}/qc/fastqc_raw" \
    --outdir "${RESULTS}/qc" \
    --filename multiqc_raw \
    --title  "RNA-seq: Raw Read QC" \
    --quiet
echo "    Report → ${RESULTS}/qc/multiqc_raw.html"

# =============================================================================
#  STEP 2  –  Align with STAR  (2-pass mode for better splice junction detection)
# =============================================================================
echo ">>> [2/5] STAR alignment (2-pass)"

# --- Pass 1: collect novel splice junctions ---
SJ_FILES=()
for SAMPLE in "${SAMPLES[@]}"; do
    OUT_PREFIX="${RESULTS}/aligned/${SAMPLE}_pass1/"
    mkdir -p "${OUT_PREFIX}"
    STAR \
        --runThreadN            "${THREADS}" \
        --genomeDir             "${STAR_INDEX}" \
        --readFilesIn           "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
                                "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" \
        --readFilesCommand      zcat \
        --outSAMtype            BAM SortedByCoordinate \
        --outSAMattributes      NH HI AS NM MD \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin    8 \
        --alignSJDBoverhangMin  1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin        20 \
        --alignIntronMax        1000000 \
        --alignMatesGapMax      1000000 \
        --outFileNamePrefix     "${OUT_PREFIX}" \
        --outSAMmode            None          # skip BAM writing in pass 1
    SJ_FILES+=("${OUT_PREFIX}SJ.out.tab")
done

# --- Pass 2: re-align using all junctions ---
for SAMPLE in "${SAMPLES[@]}"; do
    OUT_PREFIX="${RESULTS}/aligned/${SAMPLE}/"
    mkdir -p "${OUT_PREFIX}"
    echo "    Aligning ${SAMPLE} (pass 2) …"
    STAR \
        --runThreadN            "${THREADS}" \
        --genomeDir             "${STAR_INDEX}" \
        --sjdbFileChrStartEnd   "${SJ_FILES[@]}" \
        --readFilesIn           "${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz" \
                                "${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz" \
        --readFilesCommand      zcat \
        --outSAMtype            BAM SortedByCoordinate \
        --outSAMattributes      NH HI AS NM MD \
        --outSAMunmapped        Within \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin    8 \
        --alignSJDBoverhangMin  1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin        20 \
        --alignIntronMax        1000000 \
        --alignMatesGapMax      1000000 \
        --quantMode             GeneCounts \
        --outFileNamePrefix     "${OUT_PREFIX}" \
        --outReadsUnmapped      Fastx

    # Index BAM
    samtools index "${OUT_PREFIX}Aligned.sortedByCoord.out.bam"

    # Flagstat for QC
    samtools flagstat \
        "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" \
        > "${RESULTS}/qc/${SAMPLE}_flagstat.txt"

    MAPPED=$(grep "mapped (" "${RESULTS}/qc/${SAMPLE}_flagstat.txt" \
             | head -1 | awk '{print $5}' | tr -d '()')
    echo "    ${SAMPLE}: mapped ${MAPPED}"
done

# Collect STAR alignment summaries with MultiQC
multiqc \
    "${RESULTS}/aligned" \
    --outdir "${RESULTS}/qc" \
    --filename multiqc_star \
    --title  "RNA-seq: STAR Alignment QC" \
    --quiet

# =============================================================================
#  STEP 3  –  featureCounts  (gene-level counts)
# =============================================================================
echo ">>> [3/5] featureCounts"

BAM_FILES=()
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILES+=("${RESULTS}/aligned/${SAMPLE}/Aligned.sortedByCoord.out.bam")
done

featureCounts \
    -T  "${THREADS}" \
    -a  "${GTF}" \
    -o  "${RESULTS}/counts/featureCounts_raw.txt" \
    -p  --countReadPairs \
    -s  "${STRANDEDNESS}" \
    -B  \
    -C  \
    --fracOverlap 0.2 \
    "${BAM_FILES[@]}"

echo "    Raw counts → ${RESULTS}/counts/featureCounts_raw.txt"

# =============================================================================
#  STEP 4  –  Reformat featureCounts output to tidy CSV
# =============================================================================
echo ">>> [4/5] Reformatting count matrix"

python3 - <<'PYEOF'
import os, re
import pandas as pd

fc = pd.read_csv("results/counts/featureCounts_raw.txt", sep="\t", comment="#")
fc = fc.set_index("Geneid")

# Gene lengths (used for TPM computation downstream)
gene_lengths = fc["Length"].copy()
gene_lengths.to_csv("data/gene_lengths.csv", header=True)

# Raw counts only (drop annotation columns)
count_cols = [c for c in fc.columns if c.endswith(".bam")]
counts = fc[count_cols].copy()

# Clean sample names:  results/aligned/ctrl_1/Aligned... → ctrl_1
counts.columns = [re.search(r"/([^/]+)/Aligned", c).group(1) for c in counts.columns]
counts.index.name = "gene_id"
counts.to_csv("data/raw_counts.csv")

print(f"Count matrix : {counts.shape[0]:,} genes × {counts.shape[1]} samples  →  data/raw_counts.csv")
print(f"Gene lengths : {len(gene_lengths):,} genes  →  data/gene_lengths.csv")
PYEOF

# =============================================================================
#  STEP 5  –  Final MultiQC  (everything in one report)
# =============================================================================
echo ">>> [5/5] Final MultiQC report"
multiqc \
    "${RESULTS}/qc" \
    "${RESULTS}/counts" \
    --outdir "${RESULTS}/qc" \
    --filename multiqc_final \
    --title  "RNA-seq: Full Pipeline QC" \
    --quiet

echo ""
echo "============================================================"
echo "  Upstream processing complete!"
echo "============================================================"
echo "  Count matrix    :  data/raw_counts.csv"
echo "  Gene lengths    :  data/gene_lengths.csv"
echo "  Final QC report :  results/qc/multiqc_final.html"
echo ""
echo "  Next → open notebooks/RNAseq_pipeline_tutorial.ipynb"
echo "============================================================"
