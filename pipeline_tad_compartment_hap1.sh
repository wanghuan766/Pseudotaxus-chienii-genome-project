#!/bin/bash
#SBATCH --job-name=3d_hic_hap1
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=hic_3d_hap1_%j.log

set -e
set -o pipefail

# ==============================================================================
# 0. ENVIRONMENT SETUP & PATHS
# ==============================================================================
THREADS=48
GENOME_HAP1="../01_Genome_Assembly/Pseudotaxus_chienii.hap1.fasta"
HIC_R1="../data/hic/P_chienii_HiC_R1.fastq.gz"
HIC_R2="../data/hic/P_chienii_HiC_R2.fastq.gz"
REGION="Pca08:1019626736-1100348071"

echo "=== [START] ALL-IN-ONE TAD & COMPARTMENT PIPELINE FOR HAP1 ==="

# ==============================================================================
# 1. HIC-PRO: ALIGNMENT & MATRIX GENERATION (50kb and 500kb)
# ==============================================================================
echo "--- [STAGE 1] Running HiC-Pro mapping and matrix normalization ---"
# Note: You need to configure the config-hicpro_hap1.txt template before running.
# Ensure restriction fragment sites are set to DpnII as stated in the manuscript.
HiC-Pro -i ../data/hic/ -o hicpro_hap1_out -c config-hicpro_hap1.txt -p

# Convert the normalized ICE/KR matrix from HiC-Pro to HiCExplorer .h5 format
# We will extract 50 kb resolution for TADs and 500 kb resolution for Compartments
hicConvertFormat --matrix ./hicpro_hap1_out/hic_results/matrix/raw/50000/matrix_50000.matrix \
                 --bed ./hicpro_hap1_out/hic_results/matrix/raw/50000/matrix_50000_abs.bed \
                 --inputFormat hicpro --outputFormat h5 \
                 -o matrix_hap1_50kb.h5

hicConvertFormat --matrix ./hicpro_hap1_out/hic_results/matrix/raw/500000/matrix_500000.matrix \
                 --bed ./hicpro_hap1_out/hic_results/matrix/raw/500000/matrix_500000_abs.bed \
                 --inputFormat hicpro --outputFormat h5 \
                 -o matrix_hap1_500kb.h5

# ==============================================================================
# 2. HICEXPLORER: A/B COMPARTMENT ANALYSIS (hicPCA)
# ==============================================================================
echo "--- [STAGE 2] Delineating A/B Chromatin Compartments via hicPCA ---"
# Differentiating transcriptionally active A from repressed B compartments using eigenvectors
hicPCA --matrix matrix_hap1_500kb.h5 \
       --numberOfEigenvectors 1 \
       --format bigwig \
       --outputFileName P_chienii_hap1_pca.bw

# ==============================================================================
# 3. HICEXPLORER: TAD & SUB-TAD BOUNDARY DETECTION (hicFindTADs)
# ==============================================================================
echo "--- [STAGE 3] Identifying TAD and hierarchical sub-TAD boundaries ---"
# Using strict parameters based on the manuscript: thresholdComparisons 0.01
# Defining major TADs (delta 0.25)
hicFindTADs --matrix matrix_hap1_50kb.h5 \
            --thresholdComparisons 0.01 \
            --delta 0.25 \
            --numberOfProcessors ${THREADS} \
            --outPrefix P_chienii_hap1_TADs

# Defining hierarchical sub-TADs (delta 0.05)
hicFindTADs --matrix matrix_hap1_50kb.h5 \
            --thresholdComparisons 0.01 \
            --delta 0.05 \
            --numberOfProcessors ${THREADS} \
            --outPrefix P_chienii_hap1_subTADs





# ==============================================================================
# 3. Generating Spatial Chromatin Visualization for Hap1
# ==============================================================================
echo "--- [STAGE 4] Generating Spatial Chromatin Visualization for Hap1 ---"
pyGenomeTracks --tracks tracks_hap1.ini \
               --region ${REGION} \
               --output P_chienii_hap1_TBA_3D_structure.pdf \
               --title "Taxane Biosynthetic Archipelago (TBA) Chromatin Architecture - Hap1" \
               --dpi 300

echo "=== [FINISHED] Plotting Pipeline for Hap1 Completed Successfully ==="



echo "=== [FINISHED] COMPLETE 3D GENOMICS and PLOTTING PIPELINE FOR HAP1 ==="