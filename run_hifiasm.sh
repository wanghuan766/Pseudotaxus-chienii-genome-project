#!/bin/bash
#SBATCH --job-name=hifiasm_Pc
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --output=hifiasm_%j.log

set -e

# Define software and resource parameters
HIFIASM="hifiasm"
THREADS=80
OUTPUT_PREFIX="P_chienii_asm"

# Define input data paths
HIFI_READS="../data/pacbio_hifi/P_chienii_hifi.fastq.gz"
HIC_R1="../data/hic/P_chienii_HiC_R1.fastq.gz"
HIC_R2="../data/hic/P_chienii_HiC_R2.fastq.gz"

echo "=== Started hifiasm integrated Haplotype-phased Assembly ==="
$HIFIASM -o ${OUTPUT_PREFIX} \
         -t ${THREADS} \
         --h1 ${HIC_R1} \
         --h2 ${HIC_R2} \
         ${HIFI_READS}

echo "=== Converting GFA to FASTA ==="
# 提取单倍型 1 (Hap1) 和 单倍型 2 (Hap2) 的 contigs
awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.hic.hap1.p_ctg.gfa > ${OUTPUT_PREFIX}.hap1.contigs.fa
awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.hic.hap2.p_ctg.gfa > ${OUTPUT_PREFIX}.hap2.contigs.fa

echo "=== Hifiasm Pipeline Completed Successfully ==="