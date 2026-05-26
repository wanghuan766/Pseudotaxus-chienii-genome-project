#!/bin/bash
#SBATCH --job-name=repeat_anno_hap2
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=repeat_hap2_%j.log

set -e

GENOME="Pseudotaxus_chienii.hap2.fasta"

echo "=== [STEP 1] Building RepeatModeler database for Hap2 ==="
BuildDatabase -name P_chienii_hap2_db ${GENOME}

echo "=== [STEP 2] Running RepeatModeler de novo TE discovery (Hap2) ==="
RepeatModeler -database P_chienii_hap2_db -pa 24 -LTRStruct

echo "=== [STEP 3] Executing RepeatMasker for softmasking (Hap2) ==="
RepeatMasker -lib P_chienii_hap2_db-families.fa \
             -xsmall \
             -gff \
             -pa 48 \
             ${GENOME}

echo "=== [FINISHED] Repeat Annotation for Hap2 Completed Successfully ==="