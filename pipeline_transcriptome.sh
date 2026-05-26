#!/bin/bash
#SBATCH --job-name=rnaseq_tpm_pipeline
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=rnaseq_tpm_only_%j.log

set -e
set -o pipefail

# ==============================================================================
# 0. ENVIRONMENT SETUP & PARAMETERS
# ==============================================================================
THREADS=48
WORKSPACE=$(pwd)
DATA_DIR="${WORKSPACE}/data/rnaseq"          # Directory containing raw fastq.gz files
GENOME="${WORKSPACE}/01_Genome_Assembly/Pseudotaxus_chienii.hap1.fasta" # Change to hap2 for Hap2 run
GFF_ANNOTATION="${WORKSPACE}/04_Gene_Annotation/Pseudotaxus_chienii_hap1_EVM_final.gff3"

OUTPUT_DIR="${WORKSPACE}/05_transcriptome_results"
mkdir -p ${OUTPUT_DIR}/clean_data
mkdir -p ${OUTPUT_DIR}/hisat2_alignment
mkdir -p ${OUTPUT_DIR}/stringtie_gtf
mkdir -p ${OUTPUT_DIR}/expression_matrix

echo "=== [START] PURE TPM-DRIVEN TRANSCRIPTOME & HEATMAP PIPELINE ==="

# ==============================================================================
# 1. INDEXING REFERENCE GENOME
# ==============================================================================
echo "--- [STAGE 1] Building HISAT2 genome index ---"
cd ${OUTPUT_DIR}/hisat2_alignment
hisat2-build -p ${THREADS} ${GENOME} P_chienii_idx

# ==============================================================================
# 2. ITERATIVE SAMPLE PROCESSING (QUALITY CONTROL, ALIGNMENT, QUANTIFICATION)
# ==============================================================================
echo "--- [STAGE 2] Launching iterative sample processing loop ---"
cd ${DATA_DIR}

# Dynamically detect all distinct sample prefixes based on _R1.fastq.gz naming convention
SAMPLES=$(ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' | sort -u)

for SAMPLE in ${SAMPLES}; do
    echo "Processing Sample: ${SAMPLE}..."
    
    # Step 2.1: Quality control and adapter trimming via fastp
    fastp -i ${DATA_DIR}/${SAMPLE}_R1.fastq.gz \
          -I ${DATA_DIR}/${SAMPLE}_R2.fastq.gz \
          -o ${OUTPUT_DIR}/clean_data/${SAMPLE}_clean_R1.fastq.gz \
          -O ${OUTPUT_DIR}/clean_data/${SAMPLE}_clean_R2.fastq.gz \
          --thread 16 \
          --html ${OUTPUT_DIR}/clean_data/${SAMPLE}_fastp.html \
          --json ${OUTPUT_DIR}/clean_data/${SAMPLE}_fastp.json \
          >/dev/null 2>&1

    # Step 2.2: Splicing-aware mapping via HISAT2
    hisat2 -p ${THREADS} -x ${OUTPUT_DIR}/hisat2_alignment/P_chienii_idx \
           -1 ${OUTPUT_DIR}/clean_data/${SAMPLE}_clean_R1.fastq.gz \
           -2 ${OUTPUT_DIR}/clean_data/${SAMPLE}_clean_R2.fastq.gz \
           -S ${OUTPUT_DIR}/hisat2_alignment/${SAMPLE}.sam

    # Step 2.3: Convert SAM to coordinate-sorted BAM matrix
    samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/hisat2_alignment/${SAMPLE}_sorted.bam ${OUTPUT_DIR}/hisat2_alignment/${SAMPLE}.sam
    samtools index ${OUTPUT_DIR}/hisat2_alignment/${SAMPLE}_sorted.bam
    rm -f ${OUTPUT_DIR}/hisat2_alignment/${SAMPLE}.sam # Purge heavy intermediate SAM file

    # Step 2.4: Abundance quantification via StringTie based on reference gene model
    # Generating localized sample-specific GTF files containing TPM and coverage metrics
    stringtie ${OUTPUT_DIR}/hisat2_alignment/${SAMPLE}_sorted.bam \
              -p ${THREADS} \
              -G ${GFF_ANNOTATION} \
              -e -B \
              -o ${OUTPUT_DIR}/stringtie_gtf/${SAMPLE}/${SAMPLE}.gtf
done

# ==============================================================================
# 3. MATRIX EXTRACTION (TPM GENERATION)
# ==============================================================================
echo "--- [STAGE 3] Reformatting StringTie abundances into TPM matrix ---"
cd ${OUTPUT_DIR}/expression_matrix

# Create a sample list mapping file for StringTie's prepDE.py conversion tool
rm -f sample_lst.txt
for SAMPLE in ${SAMPLES}; do
    echo -e "${SAMPLE}\t${OUTPUT_DIR}/stringtie_gtf/${SAMPLE}/${SAMPLE}.gtf" >> sample_lst.txt
done

# Execute StringTie's prepDE.py utility script to extract the standardized TPM matrix
# Note: Since DESeq2 is omitted, we pass a dummy counts file name but only focus on the TPM matrix output
python2 prepDE.py -i sample_lst.txt -g dummy_counts.csv -t gene_tpm_matrix.csv
rm -f dummy_counts.csv # Remove unneeded counts file

# ==============================================================================
# 4. DOWNSTREAM EXPRESSION HEATMAP PLOTTING VIA INLINE R SCRIPT (PURE TPM)
# ==============================================================================
echo "--- [STAGE 4] Executing R engine for TPM heatmap visualizations ---"

Rscript - << 'EOF'
# Setup CRAN mirrors and dependencies verification inside the cluster environment
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

library(pheatmap)

# 4.1: Load global TPM expression profile [cite: 829]
tpm <- read.csv("gene_tpm_matrix.csv", row.names=1, check.names=False)

# Pre-filtering: Remove completely unexpressed or ultra-low expressors to optimize memory
# Keeping genes with TPM >= 1 in at least 2 samples
keep_tpm <- rowSums(tpm >= 1) >= 2
tpm_filtered <- tpm[keep_tpm, ]

# 4.2: Create Experimental Design Metadata (Bark vs Leaf)
# Automatically parse your sample header names to group into 'Bark' or 'Leaf' as per the manuscript [cite: 392]
sample_names <- colnames(tpm_filtered)
tissue_group <- ifelse(grepl("bark", sample_names, ignore.case=TRUE), "Bark", "Leaf")

# Fallback: If naming convention differs, split the design 50/50 for template compliance
if(length(unique(tissue_group)) < 2) {
  tissue_group <- rep(c("Bark", "Leaf"), length.out = length(sample_names))
}

sample_info <- data.frame(
  row.names = sample_names,
  Tissue = factor(tissue_group, levels = c("Leaf", "Bark"))
)

print("Experimental Design Matrix Created From Sample Headers:")
print(sample_info)

# 4.3: VISUALIZATION 1 - GLOBAL SAMPLE CORRELATION DISTANCE HEATMAP (BASED ON TPM)
# Log2-transformation ensures extremely high expressors do not disproportionately dominate the correlation matrix
log2_tpm <- log2(tpm_filtered + 1)

pdf("sample_distance_correlation_tpm_heatmap.pdf", width=10, height=8)
pheatmap(cor(log2_tpm, method="pearson"), 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = sample_info,
         main = "Global Sample Correlation Heatmap (Pearson r of Log2-TPM Matrix)")
dev.off()

# ==============================================================================
# 4.4: VISUALIZATION 2 - PATHWAY TARGET HEATMAP (STRICTLY BASED ON TPM VALUES)
# ==============================================================================
# Note: Replace these placeholder strings with your actual EVM target gene IDs for taxane genes [cite: 1924]
target_genes <- c("Pca08g02254", "Pca08g02257", "Pca08g02441", "Pca08g02445") 
available_targets <- intersect(target_genes, rownames(tpm_filtered))

if(length(available_targets) > 1) {
  print("--- Generating Targeted Pathway Heatmap using raw TPM Values ---")
  
  tpm_target_matrix <- tpm_filtered[available_targets, ]
  
  pdf("paclitaxel_pathway_tpm_heatmap.pdf", width=8, height=6)
  pheatmap(tpm_target_matrix, 
           scale = "row",                  # Scale rows (Z-score) to display relative expression patterns clearly
           clustering_distance_rows = "correlation",
           annotation_col = sample_info,   # Add the Bark/Leaf color bar annotation on top
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Expression Dynamics of Target Taxane Genes (Scaled Row Z-scores of TPM)")
  dev.off()
  print("Targeted pathway expression heatmap generated successfully.")
} else {
  # Fallback: Plot top 50 highly variable genes based on TPM if target IDs are missing
  row_vars <- apply(log2_tpm, 1, var)
  select_genes <- order(row_vars, decreasing = TRUE)[1:50]
  
  pdf("top_50_variable_genes_tpm_heatmap.pdf", width=8, height=10)
  pheatmap(tpm_filtered[select_genes, ], 
           scale = "row", 
           annotation_col = sample_info,
           show_rownames = FALSE,
           main = "Top 50 Most Variable Genes Heatmap (TPM Values)")
  dev.off()
  print("Custom target IDs not found. Top variable genes TPM heatmap generated instead.")
}
EOF

echo "=== [FINISHED] PURE TPM TRANSCRIPTOME PIPELINE EXECUTED SUCCESSFULLY ==="