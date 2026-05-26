#!/bin/bash
#SBATCH --job-name=merqury_3runs
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=merqury_3runs_%j.log

set -e

# Setup environment variables for Merqury
export MERQURY="../software/merqury"

# Define input assembly filenames and read directories
HAP1_GENOME="Pseudotaxus_chienii.hap1.fasta"
HAP2_GENOME="Pseudotaxus_chienii.hap2.fasta"
MGI_READS_DIR="../data/mgi_short_reads"

# ==============================================================================
# STEP 1: Generate the 21-mer database from MGI short reads (Shared across runs)
# ==============================================================================
if [ ! -d "short_reads.meryl" ]; then
    echo "=== [STEP 1] Generating 21-mer database using meryl ==="
    meryl k=21 count output short_reads.meryl ${MGI_READS_DIR}/*.fastq.gz
fi

# ==============================================================================
# STEP 2: Run 1 - Independent Evaluation of Haplotype 1 (Hap1)
# ==============================================================================
echo "=== [STEP 2] Running Run 1: Independent Evaluation for Hap1 ==="
mkdir -p run1_hap1
cd run1_hap1
$MERQURY/merqury.sh \
    ../short_reads.meryl \
    ../${HAP1_GENOME} \
    output_hap1
cd ..

# ==============================================================================
# STEP 3: Run 2 - Independent Evaluation of Haplotype 2 (Hap2)
# ==============================================================================
echo "=== [STEP 3] Running Run 2: Independent Evaluation for Hap2 ==="
mkdir -p run2_hap2
cd run2_hap2
$MERQURY/merqury.sh \
    ../short_reads.meryl \
    ../${HAP2_GENOME} \
    output_hap2
cd ..

# ==============================================================================
# STEP 4: Run 3 - Combined Evaluation for Phased Diploid Assembly (Hap1 + Hap2)
# ==============================================================================
echo "=== [STEP 4] Running Run 3: Combined Evaluation for Phased Hap1 + Hap2 ==="
mkdir -p run3_combined
cd run3_combined
$MERQURY/merqury.sh \
    ../short_reads.meryl \
    ../${HAP1_GENOME} \
    ../${HAP2_GENOME} \
    output_combined
cd ..

echo "=== [FINISHED] All 3 Independent Merqury Evaluations Completed Successfully ==="