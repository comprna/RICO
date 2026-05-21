#!/usr/bin/env bash
# RICO local runner
#
# This script runs the core RICO steps without Nextflow and access to NCI Australia:
#   1. Look for unaligned modBAM input
#   2. Align reads to rDNA-augmented reference
#   3. Estimate rDNA copy number from 45S and SCG coverage
#   4. Generate CpG bedMethyl output
#   5. Extract per-read rDNA methylation calls
#
# You need to install the tools and update the tool paths below.

set -euo pipefail

# -------------------------
# Tool paths (replace with your tool paths after intalling)
# -------------------------
MINIMAP2="/g/data/xc17/zaka/lib/minimap2-2.30/minimap2"
SAMTOOLS="/g/data/xc17/zaka/lib/samtools-1.23/bin/samtools"
BEDTOOLS="/g/data/xc17/zaka/lib/bedtools-2.31.1/bin/bedtools"
MODKIT="/g/data/xc17/zaka/lib/modkit-0.6.1/modkit"

# -------------------------
# Default settings
# -------------------------
SPECIES="human"
SCG="2"
SAMPLESHEET="samples.tsv"
OUTDIR="results"
THREADS="4"
SORT_MEM="2G"
FORCE="false"

log() {
  echo "[$(date '+%F %T')] $*" >&2
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

usage() {
  cat <<USAGE
Usage:
  bash rico_local.sh --samplesheet samples.tsv --outdir results [options]

Required:
  --samplesheet FILE       Tab-delimited samplesheet with a file_path column

Options:
  --species human|mouse    Default: human
  --scg 1|2|3              Human SCG panel. Default: 2. Ignored for mouse sample
  --outdir DIR             Default: results
  --threads N              Default: 4
  --sort-mem SIZE          Memory per samtools sort thread. (Default: 2G)
  --force                  Re-run steps even if output files already exist. (Default: false)
  -h, --help               help message

Example:
  bash rico_local.sh --samplesheet samples.tsv --species human --scg 2 --outdir result_test --threads 8
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --scg) SCG="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --sort-mem) SORT_MEM="$2"; shift 2 ;;
    --force) FORCE="true"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) fail "Unknown option: $1" ;;
  esac
done

# -------------------------
# Basic checks
# -------------------------
SPECIES="$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')"
[[ "$SPECIES" == "human" || "$SPECIES" == "mouse" ]] || fail "--species must be human or mouse"
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -ge 1 ]] || fail "--threads must be a positive integer"

[[ -x "$MINIMAP2" ]] || fail "minimap2 not found or not executable: $MINIMAP2"
[[ -x "$SAMTOOLS" ]] || fail "samtools not found or not executable: $SAMTOOLS"
[[ -x "$BEDTOOLS" ]] || fail "bedtools not found or not executable: $BEDTOOLS"
[[ -x "$MODKIT" ]] || fail "modkit not found or not executable: $MODKIT"
[[ -s "$SAMPLESHEET" ]] || fail "samplesheet not found: $SAMPLESHEET"

# Assume this script is placed in the RICO repo root or in RICO/scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -d "$SCRIPT_DIR/ref" ]]; then
  RICO_DIR="$SCRIPT_DIR"
else
  RICO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
fi

REF_DIR="$RICO_DIR/ref/$SPECIES"

if [[ "$SPECIES" == "human" ]]; then
  REF="$REF_DIR/GRCh38_rDNAx5.fa"
  RDNA_CONTIG="KY962518.1"
  [[ "$SCG" =~ ^[123]$ ]] || fail "For human, --scg must be 1, 2, or 3"
  SCG_BED="$REF_DIR/SCG-${SCG}.bed"
else
  REF="$REF_DIR/GRCm39_rDNAx5.fa"
  RDNA_CONTIG="BK000964.3_p1"
  SCG_BED="$REF_DIR/mSCG.bed"
fi

RDNA_GENE_BED="$REF_DIR/rDNAgene.bed"
REF_BASE="$(basename "$REF" .fa)"
REF_MMI="$REF_DIR/${REF_BASE}.mmi"
CG_MOTIFS_BED="$REF_DIR/${REF_BASE}.CG_motifs.bed"
SCG_NAME="$(basename "$SCG_BED" .bed)"

[[ -s "$REF" ]] || fail "reference FASTA not found: $REF"
[[ -s "$RDNA_GENE_BED" ]] || fail "rDNAgene BED not found: $RDNA_GENE_BED"
[[ -s "$SCG_BED" ]] || fail "SCG BED not found: $SCG_BED"

mkdir -p "$OUTDIR" "$OUTDIR/mean_cn" "$OUTDIR/mod"

log "RICO local run started"
log "species=$SPECIES"
log "samplesheet=$SAMPLESHEET"
log "outdir=$OUTDIR"
log "threads=$THREADS"
log "reference=$REF"
log "rDNA_contig=$RDNA_CONTIG"
log "SCG_bed=$SCG_BED"

# -------------------------
# Prepare reference files
# -------------------------
if [[ "$FORCE" == "false" && -s "$REF_MMI" ]]; then
  log "[index] Found $REF_MMI; skipping"
else
  log "[index] Building minimap2 index"
  "$MINIMAP2" -d "$REF_MMI" "$REF"
  [[ -s "$REF_MMI" ]] || fail "index was not created: $REF_MMI"
fi

if [[ "$FORCE" == "false" && -s "$CG_MOTIFS_BED" ]]; then
  log "[CG motifs] Found $CG_MOTIFS_BED; skipping"
else
  log "[CG motifs] Creating CG motif BED"
  "$MODKIT" motif bed "$REF" CG 0 > "$CG_MOTIFS_BED"
  [[ -s "$CG_MOTIFS_BED" ]] || fail "CG motif BED was not created: $CG_MOTIFS_BED"
fi

# -------------------------
# Read samplesheet
# -------------------------
FILE_COL=$(head -n 1 "$SAMPLESHEET" | awk -F '\t' '{for(i=1;i<=NF;i++) if($i=="file_path") print i}')
[[ -n "$FILE_COL" ]] || fail "samplesheet must contain a column called file_path"

mapfile -t BAM_FILES < <(awk -F '\t' -v c="$FILE_COL" 'NR > 1 && $c != "" {print $c}' "$SAMPLESHEET")
[[ "${#BAM_FILES[@]}" -gt 0 ]] || fail "no BAM files found in samplesheet"

# -------------------------
# Process each sample
# -------------------------
for UNALIGNED_BAM in "${BAM_FILES[@]}"; do

  if [[ ! -s "$UNALIGNED_BAM" ]]; then
    log "[sample] Missing unaligned BAM; skipping: $UNALIGNED_BAM"
    continue
  fi

  SAMPLE="$(basename "$UNALIGNED_BAM" .bam)"
  log "[sample] Starting $SAMPLE"

  ALIGNED_BAM="$OUTDIR/${SAMPLE}_${REF_BASE}.bam"
  ALIGNED_BAI="$ALIGNED_BAM.bai"

  CN_SUMMARY="$OUTDIR/mean_cn/${SAMPLE}_${REF_BASE}.CN.mean.txt"
  RDNA_COV="$OUTDIR/mean_cn/${SAMPLE}_${REF_BASE}.45S.coverage.tsv"
  SCG_COV="$OUTDIR/mean_cn/${SAMPLE}_${REF_BASE}_${SCG_NAME}_perBase.coverage"
  RDNA_ONLY_BAM="$OUTDIR/mean_cn/${SAMPLE}_${REF_BASE}_rDNAonly.bam"

  CPG_BED="$OUTDIR/mod/${SAMPLE}_${REF_BASE}.CpG.bed"
  RDNA_CALLS="$OUTDIR/mod/${SAMPLE}_${REF_BASE}.rDNA.calls.tsv"

  # 1. Align unaligned modBAM
  if [[ "$FORCE" == "false" && -s "$ALIGNED_BAM" && -s "$ALIGNED_BAI" ]]; then
    log "[align] Found aligned BAM and index; skipping"
  else
    log "[align] Aligning unaligned modBAM"
    "$SAMTOOLS" fastq -T "MM,ML" "$UNALIGNED_BAM" \
      | "$MINIMAP2" -t "$THREADS" -ax map-ont -Y -y --secondary=no "$REF_MMI" - \
      | "$SAMTOOLS" sort -@ "$THREADS" -m "$SORT_MEM" -o "$ALIGNED_BAM" -

    [[ -s "$ALIGNED_BAM" ]] || fail "aligned BAM was not created: $ALIGNED_BAM"

    log "[align] Indexing aligned BAM"
    "$SAMTOOLS" index "$ALIGNED_BAM"
    [[ -s "$ALIGNED_BAI" ]] || fail "BAM index was not created: $ALIGNED_BAI"
  fi

  # 2. Estimate rDNA copy number
  if [[ "$FORCE" == "false" && -s "$CN_SUMMARY" ]]; then
    log "[mean CN] Found CN summary; skipping"
  else
    log "[mean CN] Extracting rDNA-only alignments"
    "$SAMTOOLS" view -b "$ALIGNED_BAM" "$RDNA_CONTIG" > "$RDNA_ONLY_BAM"
    [[ -s "$RDNA_ONLY_BAM" ]] || fail "rDNA-only BAM was not created: $RDNA_ONLY_BAM"

    log "[mean CN] Calculating 45S mean coverage"
    "$BEDTOOLS" bamtobed -i "$RDNA_ONLY_BAM" \
      | "$BEDTOOLS" coverage -a "$RDNA_GENE_BED" -b - -mean > "$RDNA_COV"
    [[ -s "$RDNA_COV" ]] || fail "rDNA coverage file was not created: $RDNA_COV"

    log "[mean CN] Calculating SCG per-base coverage"
    "$SAMTOOLS" depth -b "$SCG_BED" "$ALIGNED_BAM" > "$SCG_COV"
    [[ -s "$SCG_COV" ]] || fail "SCG coverage file was not created: $SCG_COV"

    RDNA_SUM=$(awk '{sum += $NF} END {printf "%.2f", sum + 0}' "$RDNA_COV")
    SCG_MEAN=$(awk '{sum += $3; n += 1} END {if(n > 0) printf "%.2f", sum / n; else printf "0.00"}' "$SCG_COV")
    CN=$(awk -v r="$RDNA_SUM" -v s="$SCG_MEAN" 'BEGIN {if(s > 0) printf "%.2f", r / s; else printf "NA"}')

    {
      echo -e "bam\trDNA_contig\trDNAgene_bed\tSCG_bed\trDNA_SUM\tSCG_MEAN\tCN"
      echo -e "$(basename "$ALIGNED_BAM")\t$RDNA_CONTIG\t$RDNA_GENE_BED\t$SCG_BED\t$RDNA_SUM\t$SCG_MEAN\t$CN"
    } > "$CN_SUMMARY"

    [[ -s "$CN_SUMMARY" ]] || fail "CN summary was not created: $CN_SUMMARY"
  fi

  # 3. Create CpG bedMethyl output
  if [[ "$FORCE" == "false" && -s "$CPG_BED" ]]; then
    log "[mod pileup] Found CpG bedMethyl; skipping"
  else
    log "[mod pileup] Creating CpG bedMethyl"
    "$MODKIT" pileup "$ALIGNED_BAM" "$CPG_BED" \
      --ref "$REF" \
      --modified-bases 5mC \
      --cpg \
      --combine-strands \
      --threads "$THREADS"

    [[ -s "$CPG_BED" ]] || fail "CpG bedMethyl file was not created: $CPG_BED"
  fi

  # 4. Extract per-read rDNA methylation calls
  if [[ "$FORCE" == "false" && -s "$RDNA_CALLS" ]]; then
    log "[mod extract] Found per-read rDNA calls; skipping"
  else
    log "[mod extract] Extracting per-read rDNA methylation calls"
    "$MODKIT" extract calls "$ALIGNED_BAM" "$RDNA_CALLS" \
      -t "$THREADS" \
      --region "$RDNA_CONTIG" \
      --ref "$REF" \
      --include-bed "$CG_MOTIFS_BED" \
      --mapped-only

    [[ -s "$RDNA_CALLS" ]] || fail "per-read rDNA calls were not created: $RDNA_CALLS"
  fi

  log "[sample] Finished $SAMPLE"
done

log "RICO local run finished"
