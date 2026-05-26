#!/bin/bash
#SBATCH --job-name=func_pipeline_hap1
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=func_anno_hap1_%j.log

set -e
set -o pipefail

# ==============================================================================
# 0. ENVIRONMENT SETUP & PATHS
# ==============================================================================
THREADS=48
PROTEINS="../04_Gene_Annotation/Pseudotaxus_chienii_hap1_EVM_final.proteins.fa"

# Database Paths (Please update these to your cluster's actual database locations)
NR_DB="../data/databases/nr/nr.dmnd"
SWISSPROT_DB="../data/databases/swissprot/swissprot.dmnd"
EGGNOG_DB="../data/databases/eggnog/eggnog_5.0.dmnd"
EGGNOG_DATA_DIR="../data/databases/eggnog/data"

echo "=== [START] ALL-IN-ONE FUNCTIONAL ANNOTATION PIPELINE FOR HAP1 ==="

# ==============================================================================
# 1. DIAMOND BLASTP AGAINST NCBI NR DATABASE
# ==============================================================================
echo "--- [STAGE 1] Running DIAMOND blastp against NCBI NR Database ---"
diamond blastp --db ${NR_DB} \
               --query ${PROTEINS} \
               --out P_chienii_hap1_diamond_nr.txt \
               --evalue 1e-5 \
               --top 5 \
               --threads ${THREADS} \
               --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

# ==============================================================================
# 2. DIAMOND BLASTP AGAINST SWISS-PROT DATABASE
# ==============================================================================
echo "--- [STAGE 2] Running DIAMOND blastp against Swiss-Prot Database ---"
diamond blastp --db ${SWISSPROT_DB} \
               --query ${PROTEINS} \
               --out P_chienii_hap1_diamond_swissprot.txt \
               --evalue 1e-5 \
               --top 5 \
               --threads ${THREADS} \
               --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

# ==============================================================================
# 3. EGGNOG-MAPPER FOR EGGNOG ORTHOLOGY & PATHWAY ANNOTATION
# ==============================================================================
echo "--- [STAGE 3] Running eggNOG-mapper ---"
emapper.py -i ${PROTEINS} \
           --output P_chienii_hap1_eggnog \
           -m diamond \
           --dmnd_db ${EGGNOG_DB} \
           --data_dir ${EGGNOG_DATA_DIR} \
           --cpu ${THREADS} \
           --override

# ==============================================================================
# 4. INTERPROSCAN FOR INTERPRO, PFAM, GO, AND KEGG MAPPING
# ==============================================================================
echo "--- [STAGE 4] Running InterProScan (Includes Pfam, GO terms, and Pathways) ---"
# InterProScan natively supports multi-database integration. 
# --applications Pfam specifies Pfam-A domains. -goterms and -pa capture GO and KEGG pathways.
interproscan.sh -i ${PROTEINS} \
                -f TSV \
                -b P_chienii_hap1_interproscan_out \
                --applications Pfam \
                -goterms \
                -pa \
                -cpu ${THREADS}

echo "=== [FINISHED] COMPLETE FUNCTIONAL ANNOTATION PIPELINE FOR HAP1 ==="