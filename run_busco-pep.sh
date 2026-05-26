#!/bin/bash
#SBATCH --job-name=busco_protein_eval
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=busco_prot_%j.log

set -e

# Define software and lineage parameters
BUSCO="busco"
LINEAGE="embryophyta_odb10"

# Define paths to your predicted protein FASTA files
HAP1_PROTEINS="../04_Gene_Annotation/P_chienii_hap1.proteins.fa"
HAP2_PROTEINS="../04_Gene_Annotation/P_chienii_hap2.proteins.fa"
COMBINED_PROTEINS="../04_Gene_Annotation/P_chienii_combined_hap1_hap2.proteins.fa"

# 1. Evaluation for Haplotype 1 (Hap1) protein set
echo "=== Running BUSCO protein mode evaluation for Hap1 ==="
$BUSCO -m proteins \
       -i ${HAP1_PROTEINS} \
       -o P_chienii_Hap1_Prot_BUSCO \
       -l ${LINEAGE} \
       --cpu 32 \
       --offline

# 2. Evaluation for Haplotype 2 (Hap2) protein set
echo "=== Running BUSCO protein mode evaluation for Hap2 ==="
$BUSCO -m proteins \
       -i ${HAP2_PROTEINS} \
       -o P_chienii_Hap2_Prot_BUSCO \
       -l ${LINEAGE} \
       --cpu 32 \
       --offline

# 3. Evaluation for the merged Haplotype 1 and Haplotype 2 combined protein set
echo "=== Running BUSCO protein mode evaluation for Merged Hap1+Hap2 ==="
$BUSCO -m proteins \
       -i ${COMBINED_PROTEINS} \
       -o P_chienii_Combined_Prot_BUSCO \
       -l ${LINEAGE} \
       --cpu 32 \
       --offline

echo "=== [FINISHED] BUSCO Protein Completeness Assessment Completed Successfully ==="