#!/bin/bash
#SBATCH --job-name=repeat_anno_hap1
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=repeat_hap1_%j.log

set -e

GENOME="Pseudotaxus_chienii.hap1.fasta"

echo "=== [STEP 1] Building RepeatModeler database for Hap1 ==="
BuildDatabase -name P_chienii_hap1_db ${GENOME}

echo "=== [STEP 2] Running RepeatModeler de novo TE discovery (Hap1) ==="
# Using 48 threads globally via -pa parameter split
RepeatModeler -database P_chienii_hap1_db -pa 24 -LTRStruct


echo "=== [STEP 3] Executing RepeatMasker for softmasking (Hap1) ==="
RepeatMasker -lib  P_chienii_hap1_db-families.fa \
             -xsmall \
             -gff \
             -pa 48 \
             ${GENOME}

echo "=== [FINISHED] Repeat Annotation for Hap1 Completed Successfully ==="