#!/bin/bash
#SBATCH --job-name=3ddna_scaffold
#SBATCH --partition=SMP
#SBATCH --cpus-per-task=64
#SBATCH --output=3ddna_%j.log

set -e

# Setup environment variables for software directories
JUICER_DIR="../software/juicer"
THREE_DDNA_DIR="../software/3d-dna"
REFERENCE="P_chienii_asm.hap1.contigs.fa"
HIC_R1="../data/hic/P_chienii_HiC_R1.fastq.gz"
HIC_R2="../data/hic/P_chienii_HiC_R2.fastq.gz"

echo "=== [STEP 1] Generating Juicer restriction site positions and genome index ==="
bwa index ${REFERENCE}
python ${JUICER_DIR}/misc/generate_site_positions.py DpnII P_chienii_asm ${REFERENCE}
cut -f1,2 ${REFERENCE}.fai > P_chienii_asm.chrom.sizes

echo "=== [STEP 2] Running Juicer pipeline for Hi-C read alignment ==="
${JUICER_DIR}/scripts/juicer.sh \
    -g P_chienii_asm \
    -z ${REFERENCE} \
    -y P_chienii_asm_DpnII.txt \
    -p P_chienii_asm.chrom.sizes \
    -t 64

echo "=== [STEP 3] Running 3D-DNA pipeline for chromosome-level scaffolding ==="
# Anchoring contigs into 12 target chromosomes as specified in the manuscript
${THREE_DDNA_DIR}/run-asm-pipeline.sh \
    -m input \
    -r 0 \
    -i 12 \
    ${REFERENCE} \
    aligned/merged_nodups.txt