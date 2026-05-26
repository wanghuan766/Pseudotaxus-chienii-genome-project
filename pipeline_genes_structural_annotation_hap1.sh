#!/bin/bash
#SBATCH --job-name=gff_pipeline_hap1
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=gene_anno_hap1_%j.log

set -e
set -o pipefail

# ==============================================================================
# 0. ENVIRONMENT SETUP & PATHS
# ==============================================================================
THREADS=48
GENOME="../01_Genome_Assembly/Pseudotaxus_chienii.hap1.fasta"
MASKED_GENOME="../03_Repeat_Annotation/Pseudotaxus_chienii.hap1.fasta.masked"
FASTQ_DIR="../data/rnaseq"
PROTEINS="../data/orthodb/lineage_proteins.fa"
EVM_DIR="../software/EvidenceModeler"

# Reference genomes for GeMoMa (e.g., Arabidopsis thaliana)
REF_GFF3="../data/references/Arabidopsis_thaliana.TAIR10.gff3"
REF_FASTA="../data/references/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

echo "=== [START] ALL-IN-ONE STRUCTURAL ANNOTATION PIPELINE FOR HAP1 ==="

# ==============================================================================
# 1. TRANSCRIPTOME ASSEMBLY VIA FASTP, HISAT2, AND STRINGTIE
# ==============================================================================
echo "--- [STAGE 1] Running fastp quality control ---"
fastp -i ${FASTQ_DIR}/RNA_R1.fastq.gz -I ${FASTQ_DIR}/RNA_R2.fastq.gz \
      -o RNA_clean_R1.fastq.gz -O RNA_clean_R2.fastq.gz \
      --thread 16 --html fastp_hap1.html

echo "--- [STAGE 1] Building HISAT2 Index ---"
hisat2-build -p ${THREADS} ${GENOME} P_chienii_hap1_idx

echo "--- [STAGE 1] Aligning RNA-seq reads via HISAT2 ---"
hisat2 -p ${THREADS} -x P_chienii_hap1_idx \
       -1 RNA_clean_R1.fastq.gz -2 RNA_clean_R2.fastq.gz \
       -S P_chienii_hap1_rnaseq.sam

echo "--- [STAGE 1] Sorting and indexing BAM file ---"
samtools sort -@ ${THREADS} -o aligned_merged_hap1.bam P_chienii_hap1_rnaseq.sam
samtools index aligned_merged_hap1.bam
rm P_chienii_hap1_rnaseq.sam # Clean up large intermediate file

echo "--- [STAGE 1] Transcriptome assembly via StringTie ---"
stringtie aligned_merged_hap1.bam \
          -p ${THREADS} \
          -o P_chienii_hap1_stringtie.gtf \
          -v

# ==============================================================================
# 2. AB INITIO PREDICTION VIA BRAKER3 (INTEGRATING RNA-SEQ & PROTEIN EVIDENCE)
# ==============================================================================
echo "--- [STAGE 2] Executing BRAKER3 ab initio prediction ---"
braker.pl --genome=${MASKED_GENOME} \
          --bam=aligned_merged_hap1.bam \
          --prot_seq=${PROTEINS} \
          --threads=${THREADS} \
          --workingdir=braker3_hap1_out \
          --gff3

# ==============================================================================
# 3. HOMOLOGY-BASED PREDICTION VIA GEMOMA
# ==============================================================================
echo "--- [STAGE 3] Running GeMoMa homology prediction ---"
java -Xmx64G -jar GeMoMa-1.9.jar CLI GeMoMa \
    g=${MASKED_GENOME} \
    s=own \
    a=${REF_GFF3} \
    f=${REF_FASTA} \
    outdir=gemoma_hap1_Ath_out \
    threads=${THREADS}

# ==============================================================================
# 4. CONVERGING EVIDENCES VIA EVIDENCEMODELER (EVM)
# ==============================================================================
echo "--- [STAGE 4] Preparing configuration for evidence weights ---"
cat << EOF > weights_hap1.txt
ABINITIO_PREDICTION	BRAKER	8
TRANSCRIPT	StringTie	9
PROTEIN	GeMoMa	5
EOF

echo "--- [STAGE 4] Converting StringTie GTF to EVM compatible GFF3 ---"
# Utilizing the utility script included within the EVM package
$EVM_DIR/EvmUtils/write_EVM_inputs.pl \
    --genome ${MASKED_GENOME} \
    --gene_predictions ./braker3_hap1_out/braker.gff3 \
    --protein_alignments ./gemoma_hap1_Ath_out/gemoma_predictions.gff3 \
    --transcript_alignments ./P_chienii_hap1_stringtie.gtf \
    --output_file_pre evm_input_hap1

echo "--- [STAGE 4] Running final EvidenceModeler consensus integration ---"
$EVM_DIR/EvidenceModeler --model_org_type OTHER \
                         --genome ${MASKED_GENOME} \
                         --weights ./weights_hap1.txt \
                         --gene_predictions ./braker3_hap1_out/braker.gff3 \
                         --protein_alignments ./gemoma_hap1_Ath_out/gemoma_predictions.gff3 \
                         --transcript_alignments ./P_chienii_hap1_stringtie.gtf \
                         --segmentSize 500000 --overlapSize 10000 --write_executions \
                         --output_file Pseudotaxus_chienii_hap1_EVM_final.gff3

echo "=== [FINISHED] COMPLETE STRUCTURAL ANNOTATION PIPELINE FOR HAP1 ==="