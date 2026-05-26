#!/bin/bash
#SBATCH --job-name=comp_genomics_pipeline
#SBATCH --partition=SMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=comp_genomics_%j.log

set -e
set -o pipefail

# ==============================================================================
# 0. ENVIRONMENT SETUP & CONSTANTS
# ==============================================================================
THREADS=48
WORKSPACE=$(pwd)
DATA_DIR="${WORKSPACE}/data"
PEPTIDES_DIR="${DATA_DIR}/peptides"
SOFTWARE_DIR="${WORKSPACE}/software"

# Define Tool Paths (Ensure these are in your $PATH or update accordingly)
PROTTEST_JAR="${SOFTWARE_DIR}/prottest-3.4.2/prottest-3.4.2.jar"
MCMCTREE_BIN="${SOFTWARE_DIR}/paml4.9j/src/mcmctree"

echo "=== [START] ALL-IN-ONE COMPARATIVE GENOMICS PIPELINE ==="

# ==============================================================================
# STAGE 1: ORTHOFINDER ORTHOLOGY CLUSTERING
# ==============================================================================
echo "--- [STAGE 1] Running OrthoFinder clustering ---"
# Expects 12 species .pep.fa files inside ${PEPTIDES_DIR}
orthofinder -t ${THREADS} -a ${THREADS} -f ${PEPTIDES_DIR} -S blast

# Dynamically locate the latest OrthoFinder results directory
ORTHO_OUT_DIR=$(find ${PEPTIDES_DIR}/OrthoFinder/ -maxdepth 1 -type d -name "Results_*" | sort | tail -n 1)
echo "OrthoFinder results directory located at: ${ORTHO_OUT_DIR}"

# ==============================================================================
# STAGE 2: SINGLE-COPY ORTHOLOG EXTRACTION & MATRIX CONSTRUCTION
# ==============================================================================
echo "--- [STAGE 2] Extracting and trimming single-copy orthologs ---"
SINGLE_COPY_DIR="${ORTHO_OUT_DIR}/Single_Copy_Orthologue_Sequences"
ALIGN_OUT_DIR="${WORKSPACE}/02_alignment_processing"
mkdir -p ${ALIGN_OUT_DIR}

# Copy raw single-copy sequences to the local tracking folder
cp ${SINGLE_COPY_DIR}/*.fa ${ALIGN_OUT_DIR}/
cd ${ALIGN_OUT_DIR}

# Execute serialized trimming pipeline across all single-copy groups
for fa_file in *.fa; do
    [ -e "$fa_file" ] || continue
    
    # Step 2.1: Multiple sequence alignment via muscle
    muscle -in "${fa_file}" -out "${fa_file}.aligned" > /dev/null 2>&1
    
    # Step 2.2: Conserved block selection via Gblocks
    Gblocks "${fa_file}.aligned" -b4=5 -b5=h -t=p -e=.2 > /dev/null 2>&1
    
    # Step 2.3: Order sorting by sequence IDs to ensure correct concatenation coordinate
    seqkit sort "${fa_file}.aligned-gb" -o "${fa_file}.sorted"
    
    # Step 2.4: Convert multi-line sequences to single-line format
    seqkit seq "${fa_file}.sorted" -w 0 -o "${fa_file}.final_seq"
done

# Super-matrix concatenation step
echo "--- [STAGE 2] Concatenating single-copy groups into supermatrix ---"
paste -d " " *.final_seq > combined_raw.fa
sed "s/ //g" combined_raw.fa > all_concatenated.fa

# Clean up structural intermediate files
rm -f *.aligned *.aligned-gb* *.sorted *.final_seq combined_raw.fa

# ==============================================================================
# STAGE 3: RAXML MAXIMUM LIKELIHOOD PHYLOGENETIC RECONSTRUCTION
# ==============================================================================
echo "--- [STAGE 3] Preparing file alignment for RAxML ---"
cd ${WORKSPACE}
mkdir -p 03_phylogeny
mv ${ALIGN_OUT_DIR}/all_concatenated.fa 03_phylogeny/all.fa
cd 03_phylogeny

# Step 3.1: Inline Python conversion execution (FASTA to PHYLIP format)
python3 -c "
import re, argparse
with open('all.fa', 'r') as fin:
    sequences = [(m.group(1), ''.join(m.group(2).split())) for m in re.finditer(r'(?m)^>([^ \n]+)[^\n]*([^>]*)', fin.read())]
with open('all.phy', 'w') as fout:
    fout.write('%d %d\n' % (len(sequences), len(sequences[0][1])))
    for item in sequences: fout.write('%-20s %s\n' % item)
"

# Step 3.2: Best amino acid substitution model selection via ProtTest
echo "--- [STAGE 3] Evaluating substitution model via ProtTest ---"
java -jar ${PROTTEST_JAR} \
     -i all.phy \
     -all-distributions \
     -F \
     -AIC -BIC \
     -tc 0.5 \
     -threads ${THREADS} \
     -o prottest.out

# Step 3.3: Constructing Maximum Likelihood Tree with 1000 Rapid Bootstraps
echo "--- [STAGE 3] Executing core RAxML tree construction ---"
# Rooted by P.patens as specified in the manuscript
raxmlHPC-PTHREADS-SSE3 \
  -T 20 \
  -f a \
  -x 123 \
  -p 123 \
  -N 1000 \
  -m PROTGAMMAIJTTF \
  -k \
  -O \
  -o P.patens \
  -n all.tree \
  -s all.fa

# ==============================================================================
# STAGE 4: MCMCTREE DIVRGENCE TIME ESTIMATION (PAML)
# ==============================================================================
echo "--- [STAGE 4] Starting divergence time analysis via MCMCTree ---"
cd ${WORKSPACE}
mkdir -p 04_mcmctree
cd 04_mcmctree

# Link sequences and prepared tree inputs
cp ../03_phylogeny/all.phy ./tmp.phy

# CRITICAL MANUAL DEPENDENCY: User must provide a calibrated tree 'tmp.tree'
# including fossil constraints based on the TimeTree database.
if [ ! -f "tmp.tree" ]; then
    echo "WARNING: 'tmp.tree' with fossil calibrations not found in 04_mcmctree!"
    echo "Please create 'tmp.tree' manually to run MCMCTree."
else
    # Writing automated mcmctree control configuration file
    cat << EOF > mcmctree.ctl
       seed = -1
       seqfile = ./tmp.phy
       treefile = ./tmp.tree
       outfile = ./mcmctree.out
       ndata = 1
       seqtype = 2
       usedata = 1
       clock = 3
       RootAge = <5.09
       model = 0
       alpha = 0
       ncatG = 4
       cleandata = 0
       BDparas = 1 1 0
   kappa_gamma = 6 2
   alpha_gamma = 1 1
   rgene_gamma = 2 2
  sigma2_gamma = 1 10
      finetune = 0.05 0.1 0.12 0.1 0.1 0.12
         print = 1
        burnin = 4000
      sampfreq = 2
       nsample = 20000
EOF

    echo "--- [STAGE 4] Executing MCMCTree binary ---"
    ${MCMCTREE_BIN} mcmctree.ctl
fi

# ==============================================================================
# STAGE 5: CAFE5 GENE FAMILY EXPANSION & CONTRACTION ANALYSIS
# ==============================================================================
echo "--- [STAGE 5] Starting CAFE5 expansion and contraction analysis ---"
cd ${WORKSPACE}
mkdir -p 05_cafe5
cd 05_cafe5

# Extract original counts from OrthoFinder table
GENE_COUNTS_RAW="${ORTHO_OUT_DIR}/Orthogroups/Orthogroups.GeneCount.tsv"

if [ -f "${GENE_COUNTS_RAW}" ]; then
    echo "--- [STAGE 5] Formatting gene family counts matrix ---"
    # Step 5.1: Parse and reformat matrix headers via awk and sed
    awk -v OFS="\t" '{$NF=null;print $1,$0}' ${GENE_COUNTS_RAW} | \
    sed -E -e 's/Orthogroup/desc/' -e 's/_[^\t]+//g' > gene_families.txt

    # Step 5.2: Filtering out ultra-high copy clusters exceeding 100 copies
    awk 'NR==1 || ($3<100 && $4<100 && $5<100 && $6<100 && $7<100 && $8<100 && $9<100 && $10<100 && $11<100 && $12<100 && $13<100 && $14<100) {print $0}' \
    gene_families.txt > gene_families_filter.txt

    # CRITICAL MANUAL DEPENDENCY: User must provide a lambda-calibrated 'modified_tree.txt'
    if [ ! -f "modified_tree.txt" ]; then
        echo "WARNING: 'modified_tree.txt' not found in 05_cafe5!"
        echo "Please provide the branch-length calibrated tree to execute CAFE5."
    else
        echo "--- [STAGE 5] Executing CAFE5 engine ---"
        # Executing CAFE5 utilizing the Poisson distribution and discrete gamma categories
        cafe5 -i gene_families_filter.txt \
              -t modified_tree.txt \
              -p \
              -k 3 \
              -o k3p
    fi
else
    echo "ERROR: Orthogroups.GeneCount.tsv not found in OrthoFinder outputs."
fi

echo "=== [FINISHED] COMPLETE COMPARATIVE GENOMICS PIPELINE EXECUTED ==="