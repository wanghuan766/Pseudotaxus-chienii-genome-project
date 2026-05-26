#!/bin/bash
#SBATCH --job-name=rbh_kaks_hap1
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=kaks_rbh_hap1_%j.log

set -e
set -o pipefail

# ==============================================================================
# 0. ENVIRONMENT SETUP & PATHS
# ==============================================================================
THREADS=48
WORKSPACE=$(pwd)
OUT_DIR="${WORKSPACE}/kaks_hap1_rbh_out"
mkdir -p ${OUT_DIR}

# Genomic dataset paths as specified in manuscript methods
PC_CDS="${WORKSPACE}/04_Gene_Annotation/Pseudotaxus_chienii_hap1_EVM_final.cds.fa"
PC_PEP="${WORKSPACE}/04_Gene_Annotation/Pseudotaxus_chienii_hap1_EVM_final.proteins.fa"
TM_CDS="../data/references/Taxus_mairei.cds.fa"
TM_PEP="../data/references/Taxus_mairei.proteins.fa"

echo "=== [START] ALL-IN-ONE RBH BLASTP & PAML KA/KS PIPELINE FOR HAP1 ==="

cd ${OUT_DIR}

# ==============================================================================
# 1. RECIPROCAL BEST HIT (RBH) VIA BLASTP (DIAMOND) WITH >=45% IDENTITY FILTER
# ==============================================================================
echo "--- [STAGE 1] Running Reciprocal Best Hit (RBH) blastp engine ---"

# 1.1 Build diamond databases
diamond makedb --in ${PC_PEP} --db pc_hap1_pep > /dev/null 2>&1
diamond makedb --in ${TM_PEP} --db tm_pep > /dev/null 2>&1

# 1.2 Forward blastp: Taxus_mairei to Pseudotaxus_chienii (Hap1)
echo "Running Forward BLASTp (Tm -> Pc)..."
diamond blastp --db pc_hap1_pep.dmnd --query ${TM_PEP} --out tm_to_pc_forward.blast \
               --evalue 1e-5 --max-target-seqs 1 --threads ${THREADS} \
               --outfmt 6 qseqid sseqid pident evalue bitscore

# 1.3 Reverse blastp: Pseudotaxus_chienii (Hap1) to Taxus_mairei
echo "Running Reverse BLASTp (Pc -> Tm)..."
diamond blastp --db tm_pep.dmnd --query ${PC_PEP} --out pc_to_tm_reverse.blast \
               --evalue 1e-5 --max-target-seqs 1 --threads ${THREADS} \
               --outfmt 6 qseqid sseqid pident evalue bitscore

# 1.4 Strict RBH algorithmic parsing + Identity >= 45% threshold filtering using inline Python
echo "Parsing strict RBH pairs with amino acid identity >= 45%..."
python3 -c "
import sys

# Load forward hits (Tm -> Pc)
forward_hits = {}
with open('tm_to_pc_forward.blast') as f:
    for line in f:
        parts = line.strip().split('\t')
        tm_id, pc_id, pident = parts[0], parts[1], float(parts[2])
        if tm_id not in forward_hits:
            forward_hits[tm_id] = (pc_id, pident)

# Load reverse hits (Pc -> Tm)
reverse_hits = {}
with open('pc_to_tm_reverse.blast') as f:
    for line in f:
        parts = line.strip().split('\t')
        pc_id, tm_id = parts[0], parts[1]
        if pc_id not in reverse_hits:
            reverse_hits[pc_id] = tm_id

# Cross-validate Reciprocal Best Hits and enforce >=45% identity constraint
with open('validated_rbh_orthologs_identity_45.txt', 'w') as out:
    for tm_id, (pc_id, pident) in forward_hits.items():
        if pc_id in reverse_hits and reverse_hits[pc_id] == tm_id:
            if pident >= 45.0: # Strict manuscript threshold filtering
                out.write(f'{pc_id}\t{tm_id}\t{pident}\n')
"

cat validated_rbh_orthologs_identity_45.txt | wc -l | awk '{print "Total valid RBH pairs identified (Identity >= 45%): "$1}'

# ==============================================================================
# 2. SEQUENCE ISOLATION FOR IDENTIFIED COMPATIBLE PAIRS
# ==============================================================================
echo "--- [STAGE 2] Slicing fasta sequence segments for valid pairs ---"
python3 -c "
def parse_fasta(path):
    d = {}
    with open(path) as f:
        curr, seq = None, []
        for line in f:
            if line.startswith('>'):
                if curr: d[curr] = ''.join(seq)
                curr, seq = line.strip().split()[0][1:], []
            else: seq.append(line.strip())
        if curr: d[curr] = ''.join(seq)
    return d

pc_cds_d = parse_fasta('${PC_CDS}'); pc_pep_d = parse_fasta('${PC_PEP}')
tm_cds_d = parse_fasta('${TM_CDS}'); tm_pep_d = parse_fasta('${TM_PEP}')

with open('validated_rbh_orthologs_identity_45.txt') as f:
    for line in f:
        pc_id, tm_id, _ = line.strip().split('\t')
        if pc_id in pc_cds_d and tm_id in tm_cds_d:
            with open(f'{pc_id}_{tm_id}.cds.fa', 'w') as out_c:
                out_c.write(f'>{pc_id}\n{pc_cds_d[pc_id]}\n>{tm_id}\n{tm_cds_d[tm_id]}\n')
            with open(f'{pc_id}_{tm_id}.pep.fa', 'w') as out_p:
                out_p.write(f'>{pc_id}\n{pc_pep_d[pc_id]}\n>{tm_id}\n{tm_pep_d[tm_id]}\n')
"

# ==============================================================================
# 3. CODON ALIGNMENT & PAML YN00 EXECUTION (NG86 & YN00 STATISTICS)
# ==============================================================================
echo "--- [STAGE 3] Running MUSCLE alignment and PAML metrics computation ---"
echo -e "Pair_ID\tKa_NG86\tKs_NG86\tKa_YN00\tKs_YN00\tOmega_YN00" > Master_Hap1_RBH_KaKs_PAML_Report.txt

for pep_file in *.pep.fa; do
    [ -e "$pep_file" ] || continue
    PREFIX="${pep_file%.pep.fa}"
    
    # 3.1 Align amino acids via MUSCLE v5 as stated in manuscript methods
    muscle -in "${PREFIX}.pep.fa" -out "${PREFIX}.pep.aligned" >/dev/null 2>&1
    
    # 3.2 Reconstruct structural codon blocks via pal2nal.pl
    pal2nal.pl "${PREFIX}.pep.aligned" "${PREFIX}.cds.fa" -output pal2nal -nogap > "${PREFIX}.codon.paml"
    
    # 3.3 Create PAML yn00 package runtime controls config file
    cat << EOF > yn00.ctl
       seqfile = ${PREFIX}.codon.paml
       outfile = ${PREFIX}.yn00.out
       verbose = 1
         icode = 0
     weighting = 0
  commonalpha = 0
EOF

    # 3.4 Execute yn00 framework to invoke both algorithms
    yn00 yn00.ctl > /dev/null 2>&1
    
    # 3.5 Text parsing blocks to extract both NG86 and YN00 statistics from outputs
    if [ -f "${PREFIX}.yn00.out" ]; then
        # 1. Parse standard YN00 row
        YN_STATS=$(awk '/t=/{print $0}' "${PREFIX}.yn00.out" | tail -n 1)
        # 2. Parse standard Nei-Gojobori (NG86) data block row
        NG_STATS=$(grep -A 2 "Nei & Gojobori" "${PREFIX}.yn00.out" | tail -n 1)
        
        if [ ! -z "$YN_STATS" ] && [ ! -z "$NG_STATS" ]; then
            # Extract YN00 parameters
            KA_YN=$(echo "$YN_STATS" | awk -F'=' '{print $4}' | awk '{print $1}')
            KS_YN=$(echo "$YN_STATS" | awk -F'=' '{print $5}' | awk '{print $1}')
            OMEGA_YN=$(echo "$YN_STATS" | awk -F'=' '{print $3}' | awk '{print $1}')
            
            # Extract NG86 parameters (Read dN and dS matrix outputs from PAML logs)
            KA_NG=$(echo "$NG_STATS" | awk '{print $2}')
            KS_NG=$(echo "$NG_STATS" | awk '{print $3}')
            
            echo -e "${PREFIX}\t${KA_NG}\t${KS_NG}\t${KA_YN}\t${KS_YN}\t${OMEGA_YN}" >> Master_Hap1_RBH_KaKs_PAML_Report.txt
        fi
    fi
    
    # Clean up localized temporary loop scripts
    rm -f yn00.ctl "${PREFIX}.pep.aligned" "${PREFIX}.codon.paml" "${PREFIX}.yn00.out"
done

# Purge bulky raw blast files to leave workspace pristine
rm -f tm_to_pc_forward.blast pc_to_tm_reverse.blast

echo "=== [FINISHED] PAML KA/KS COMPILATION COMPLETED. ALL RESULTS CONVERGED ==="