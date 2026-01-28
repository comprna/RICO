nextflow.enable.dsl = 2

params.ref         = "/g/data/xc17/zaka/res/human/rDNAx5/GRCh38_rDNAx5.fa"
params.index_dir   = "/g/data/xc17/zaka/res/human/rDNAx5"
params.results_dir = "/g/data/xc17/zaka/nextflow/rDNA-CN-pipeline/results"

// TSV contains sample path
params.samplesheet = "/g/data/xc17/zaka/nextflow/rDNA-CN-pipeline/samples.tsv"


params.minimap2    = "/g/data/xc17/zaka/lib/minimap2_2.30/minimap2"
params.samtools    = "/g/data/xc17/zaka/lib/samtools-1.23/bin/samtools"


process index {
  tag { genome.baseName }

  publishDir "${params.index_dir}", mode: 'copy', overwrite: true

  input:
    path genome

  output:
    path "*.mmi"

  script:
    def out_mmi   = "${genome.baseName}.mmi"
    def exist_mmi = "${params.index_dir}/${out_mmi}"

    """
    set -euo pipefail
    ts(){ date +"%F %T"; }

    echo "[\$(ts)] [index] genome=${genome}" 1>&2
    echo "[\$(ts)] [index] out=${out_mmi}" 1>&2

    if [[ -s "${exist_mmi}" ]]; then
      echo "[\$(ts)] [index] FOUND ${exist_mmi} -> skipping build" 1>&2
      ln -sf "${exist_mmi}" "${out_mmi}"
    else
      echo "[\$(ts)] [index] MISSING ${exist_mmi} -> building" 1>&2
      ${params.minimap2} --version 1>&2
      ${params.minimap2} -d "${out_mmi}" "${genome}"
      echo "[\$(ts)] [index] DONE build ${out_mmi}" 1>&2
    fi
    """
}


process align {
  tag { unalign_bam.baseName }

  publishDir "${params.results_dir}", mode: 'copy', overwrite: true

  input:
    path unalign_bam
    path ref_mmi

  output:
    path "*.bam"
    path "*.bam.bai"

  script:
    // Output naming:
    //   <inputBasename>_<refBasename>.bam
    // Example:
    //   hg005_1_test.bam + GRCh38_rDNAx5.mmi -> hg005_1_test_GRCh38_rDNAx5.bam
    def sample     = unalign_bam.baseName
    def refTag     = ref_mmi.baseName
    def out_bam    = "${sample}_${refTag}.bam"
    def out_bai    = "${out_bam}.bai"
    def exist_bam  = "${params.results_dir}/${out_bam}"
    def exist_bai  = "${params.results_dir}/${out_bai}"

    """
    set -euo pipefail
    ts(){ date +"%F %T"; }

    echo "[\$(ts)] [align] input=${unalign_bam}" 1>&2
    echo "[\$(ts)] [align] index=${ref_mmi}" 1>&2
    echo "[\$(ts)] [align] out=${out_bam}" 1>&2

    if [[ -s "${exist_bam}" && -s "${exist_bai}" ]]; then
      echo "[\$(ts)] [align] FOUND results -> skipping alignment" 1>&2
      ln -sf "${exist_bam}" "${out_bam}"
      ln -sf "${exist_bai}" "${out_bai}"
      exit 0
    fi

    echo "[\$(ts)] [align] START bam->fastq (MM,ML) | minimap2 | sort" 1>&2

    ${params.samtools} fastq "${unalign_bam}" -T "MM,ML" \
      | ${params.minimap2} -t ${task.cpus} -ax map-ont "${ref_mmi}" -Y -y --secondary=no - \
      | ${params.samtools} sort -@ 8 -m 2G -o "${out_bam}"

    echo "[\$(ts)] [align] START samtools index ${out_bam}" 1>&2
    ${params.samtools} index "${out_bam}"

    echo "[\$(ts)] [align] DONE ${out_bam} (+.bai)" 1>&2
    """
}


workflow {

  ref_ch = Channel.fromPath(params.ref)
  mmi_ch = index(ref_ch)

  // Validate samplesheet param
  if( !params.samplesheet ) {
    log.warn "params.samplesheet not set; skipping alignment"
    return
  }

  // Read BAMs from TSV and fan out
  bam_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
      if( !row.file_path ) {
        throw new IllegalArgumentException("samplesheet missing 'file_path' value in a row")
      }
      row.file_path as String
    }
    .map { p -> file(p) }
    .filter { f ->
      if( !f.exists() ) {
        log.warn "Missing BAM, skipping: ${f}"
        return false
      }
      true
    }

  // Run alignments in parallel (one task per BAM)
  align(bam_ch, mmi_ch)
}

