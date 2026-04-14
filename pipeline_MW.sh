#!/bin/bash

# =============================================
# INTERACTIVE PIPELINE: MW Polyomavirus Consensus Sequence
# =============================================
# Authors: Thuany Giovana Daniel & Gabriel de Campos
# Usage:
# ./pipeline_MW.sh <R1.fastq.gz> <R2.fastq.gz> <SAMPLE_ID> <OUTPUT_DIR> [THREADS]
# Example:
# ./pipeline_MW.sh R1.fastq.gz R2.fastq.gz MW57 ./results 8
# =============================================

set -e          # Exit immediately if a command exits with a non-zero status
set -o pipefail # The return value of a pipeline is the status of the last command to exit with a non-zero status

usage() {
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <SAMPLE_ID> <OUTPUT_DIR> [OPTIONS]"
    echo ""
    echo "Required Arguments:"
    echo "  <R1.fastq.gz>      Forward reads file (R1)"
    echo "  <R2.fastq.gz>      Reverse reads file (R2)"
    echo "  <SAMPLE_ID>        Unique identifier for the sample"
    echo "  <OUTPUT_DIR>       Directory to save all results"
    echo ""
    echo "Options:"
    echo "  --threads INT      Number of threads to use (default: 4)"
    echo "  -h, --help         Show this help message"
    exit 1
}

if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
    echo "Error: Incorrect number of arguments."
    usage
fi

# --- Input Variables ---
FASTQ_R1="$1"
FASTQ_R2="$2"
SAMPLE_ID="$3"
OUTPUT_DIR="$4"
THREADS="${5:-4}"  # Default to 4 threads if not specified
VAPOR_DB="MWPyV_complete.fasta"

# --- Internal Directory Structure ---
QC_DIR="$OUTPUT_DIR/01_fastp"
REF_DIR="$OUTPUT_DIR/02_reference"
ALIGN_DIR="$OUTPUT_DIR/03_alignment"
VAR_DIR="$OUTPUT_DIR/04_variants"
CONS_DIR="$OUTPUT_DIR/05_consensus"
MEGAHIT_DIR="$OUTPUT_DIR/06_deNovoAssembly"

mkdir -p "$QC_DIR" "$REF_DIR" "$ALIGN_DIR" "$VAR_DIR" "$CONS_DIR"

echo "MW Polyomavirus Pipeline"
echo "Sample: $SAMPLE_ID"
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "---------------------------------------------------"

# =============================================
# 1️⃣ QUALITY CONTROL & TRIMMING (FASTP)
# =============================================
TRIMMED_R1="$QC_DIR/${SAMPLE_ID}_R1.trimmed.fastq.gz"
TRIMMED_R2="$QC_DIR/${SAMPLE_ID}_R2.trimmed.fastq.gz"
FASTP_HTML="$QC_DIR/${SAMPLE_ID}_fastp.html"
FASTP_JSON="$QC_DIR/${SAMPLE_ID}_fastp.json"

echo "[STEP 1/7] Trimming reads with fastp..."
fastp \
    -i "$FASTQ_R1" -I "$FASTQ_R2" \
    -o "$TRIMMED_R1" -O "$TRIMMED_R2" \
    --qualified_quality_phred 20 -l 15 \
    -w "$THREADS" \
    --cut_front --cut_tail \
    -h "$FASTP_HTML" -j "$FASTP_JSON"

# =============================================
# 2️⃣ BEST REFERENCE IDENTIFICATION (VAPOR)
# =============================================
# Check if the database exists before running VAPOR
if [ ! -f "$VAPOR_DB" ]; then
    echo "Error: Database file not found at $VAPOR_DB"
    echo "Please create the 'databases' folder and add MWPyV_complete.fasta"
    exit 1
fi

REF_RESULT="$REF_DIR/${SAMPLE_ID}_best_ref.fasta"

echo "[STEP 2/7] Identifying best reference using VAPOR..."
python vapor.py \
    -fa "$VAPOR_DB" -fq "$TRIMMED_R1" "$TRIMMED_R2" \
    | cut -f 6 | sed 's/>//g' > "$REF_DIR/${SAMPLE_ID}_ref_name.txt"

# Extract the specific reference sequence found by VAPOR
seqtk subseq "$VAPOR_DB" "$REF_DIR/${SAMPLE_ID}_ref_name.txt" > "$REF_RESULT"
REF_FASTA="$REF_RESULT"

# =============================================
# 3️⃣ REFERENCE INDEXING (BWA)
# =============================================
if [ ! -f "${REF_FASTA}.bwt" ]; then
    echo "[STEP 3/7] Indexing reference with BWA..."
    bwa index "$REF_FASTA"
else
    echo "[STEP 3/7] Index already exists. Skipping."
fi

# =============================================
# 4️⃣ READ MAPPING (BWA MEM)
# =============================================
SORTED_BAM="$ALIGN_DIR/${SAMPLE_ID}.sorted.bam"
echo "[STEP 4/7] Aligning reads with BWA..."
bwa mem -t "$THREADS" "$REF_FASTA" "$TRIMMED_R1" "$TRIMMED_R2" | \
    samtools view -@ "$THREADS" -Sb - | \
    samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
samtools index "$SORTED_BAM"

# =============================================
# 5️⃣ VARIANT CALLING (FREEBAYES): Not obligatory
# =============================================
VCF_FILE="$VAR_DIR/${SAMPLE_ID}.freebayes.vcf"
echo "[STEP 5/7] Calling variants with FreeBayes..."
freebayes -f "$REF_FASTA" "$SORTED_BAM" > "$VCF_FILE"

# =============================================
# 6️⃣ GENOME POLISHING (PILON)
# =============================================
PILON_OUT="$ALIGN_DIR/${SAMPLE_ID}_pilon"
PILON_FASTA="${PILON_OUT}.fasta"
PILON_BAM="$ALIGN_DIR/${SAMPLE_ID}_pilon_mapped.sorted.bam"

echo "[STEP 6/7] Polishing reference with Pilon..."

java -Xmx16G -jar /usr/local/src/pilon/1.24/pilon-1.24.jar \
    --genome "$REF_FASTA" \
    --frags "$SORTED_BAM" \
    --output "$PILON_OUT" \
    --fix "gaps,indels" \
    --mindepth 5 --minqual 20 --minmq 10

# Re-mapping reads to the polished (Pilon) assembly
bwa index "$PILON_FASTA"
bwa mem -t "$THREADS" "$PILON_FASTA" "$TRIMMED_R1" "$TRIMMED_R2" | \
    samtools view -@ "$THREADS" -Sb - | \
    samtools sort -o "$PILON_BAM" -
samtools index "$PILON_BAM"

# =============================================
# 7️⃣ FINAL CONSENSUS GENERATION (IVAR)
# =============================================
CONSENSUS_FA="$CONS_DIR/${SAMPLE_ID}_consensus.fasta"
IVAR_CONS="$CONS_DIR/${SAMPLE_ID}_ivar" # iVar adds .fa suffix automatically
IVAR_BAM="$CONS_DIR/${SAMPLE_ID}_final_mapped.sorted.bam"
IVAR_VARIANTS="$CONS_DIR/${SAMPLE_ID}_iVar_variants"

echo "[STEP 7/7] Generating final consensus with iVar..."
samtools mpileup -aa -A -d 0 -B -Q 0 "$PILON_BAM" | \
    ivar consensus -q 20 -m 5 -n N -p "$IVAR_CONS"

# Validate final consensus via re-mapping and variant detection
bwa index "${IVAR_CONS}.fa"
bwa mem -t "$THREADS" "${IVAR_CONS}.fa" "$TRIMMED_R1" "$TRIMMED_R2" | \
    samtools sort -o "$IVAR_BAM"
samtools index "$IVAR_BAM"

samtools mpileup -aa -A -d 0 -B -Q 0 "$IVAR_BAM" | \
    ivar variants -p "$IVAR_VARIANTS" -q 30 -t 0.03 -m 100 -r "${IVAR_CONS}.fa"

# Count confirmed SNPs (excluding Ns)
grep 'TRUE' "$IVAR_VARIANTS.tsv" | grep -vP '\tN\t' | wc -l > "$CONS_DIR/${SAMPLE_ID}.SNPsCount"

# =============================================
# 8️⃣ QUALITY METRICS CALCULATION
# =============================================
SUMMARY_FILE="$OUTPUT_DIR/${SAMPLE_ID}_metrics.txt"

# Statistical calculations
NUM_TOTAL_READS=$(samtools view -c "$SORTED_BAM")
NUM_MAPPED_READS=$(samtools view -c -F 4 "$SORTED_BAM")
MEAN_DEPTH=$(samtools depth "$SORTED_BAM" | awk '{sum+=$3} END {if(NR>0) print sum/NR; else print 0}')
MEDIAN_DEPTH=$(samtools depth "$SORTED_BAM"  |  awk '{print $3}' | sort -n | awk 'NF{a[NR]=$1;c++}END {print (c%2==0)?(a[int(c/2)+1]+a[int(c/2)])/2:a[int(c/2)+1]}')
REF_SIZE=$(awk '/^>/ {if (l) print l; l=0; next} {l+=length} END {print l}' "$REF_FASTA")
BASES_10X=$(samtools depth "$SORTED_BAM" | awk '$3>=10' | wc -l)
DEPTH_10X_PCT=$(awk -v b="$BASES_10X" -v r="$REF_SIZE" 'BEGIN{if(r>0) printf "%.2f\n",(b/r)*100; else print 0}')
N_COUNT=$(grep -v '^>' "${IVAR_CONS}.fa" | tr -d '\n' | tr -cd 'Nn' | wc -c)
REF_COVERAGE_PCT=$(awk -v r="$REF_SIZE" -v n="$N_COUNT" 'BEGIN{if(r>0) printf "%.2f\n",(r-n)/r*100; else print 0}')
N_PCT=$(awk -v r="$REF_SIZE" -v n="$N_COUNT" 'BEGIN{if(r>0) printf "%.2f\n", (n/r)*100; else print 0}')
MAPPING_PCT=$(awk -v t="$NUM_TOTAL_READS" -v m="$NUM_MAPPED_READS" 'BEGIN{if(t>0) printf "%.2f\n",(m/t)*100; else print 0}')

# Generate TSV Header and Data
echo -e "sample\ttotal_reads\tmapped_reads\tmean_depth\tmedian_depth\tdepth_10X_pct\tcoverage_pct\tN_count\tN_pct\tmapping_pct" > "$SUMMARY_FILE"
echo -e "${SAMPLE_ID}\t${NUM_TOTAL_READS}\t${NUM_MAPPED_READS}\t${MEAN_DEPTH}\t${MEDIAN_DEPTH}\t${DEPTH_10X_PCT}\t${REF_COVERAGE_PCT}\t${N_COUNT}\t${N_PCT}\t${MAPPING_PCT}" >> "$SUMMARY_FILE"

# =============================================
# 9️⃣ ALIGNMENT OF CONSENSUS VS REFERENCE
# =============================================
echo "Finalizing: Aligning consensus with reference using MAFFT..."
mafft --quiet --auto --keeplength --inputorder --6merpair --leavegappyregion \
  --addfragments "${IVAR_CONS}.fa" "$REF_FASTA" > "$OUTPUT_DIR/${SAMPLE_ID}_ivar_aligned.fasta"

echo "Pipeline completed successfully!"
echo "Final Consensus: ${IVAR_CONS}.fa"
echo "Post-Pilon BAM: $PILON_BAM"
echo "FreeBayes VCF: $VCF_FILE"
echo "Metrics: $SUMMARY_FILE"
