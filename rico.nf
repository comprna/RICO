nextflow.enable.dsl = 2

def refRoot = "${projectDir}/ref"

//-------------------------------------------------------------------------
// Run-level options
//-------------------------------------------------------------------------
params.species = params.species ?: 'human'   // human or  mouse
params.scg     = params.scg ?: 2             // only meaningful for human (1/2/3)
params.outdir  = params.outdir ?: "${projectDir}/results"
params.samplesheet = params.samplesheet ?: "${projectDir}/samples.tsv"

params.base = projectDir


 // Results folder
 params.results_dir = params.outdir
 params.mod_results = "${params.results_dir}/mod"
  

//-------------------------------------------------------------------------
// Derived paths - human or mouse
//-------------------------------------------------------------------------
def species = params.species.toString().trim().toLowerCase()

if( !(species in ['human','mouse']) ) {
  error "Invalid --species '${params.species}'. Allowed values = [human, mouse]"
}

def refDir = "${refRoot}/${species}"

// Reference
params.ref       = params.ref ?: "${refDir}/${ species=='human' ? 'GRCh38_rDNAx5.fa' : 'GRCm39_rDNAx5.fa' }"
params.index_dir = params.index_dir ?: refDir
// rDNA bed
params.rDNAgene_bed = params.rDNAgene_bed ?: "${refDir}/rDNAgene.bed"

//-------------------------------------------------------------------------
// SCG bed
// Human: 3 SCG panels to choose from [default: SCG-2]
// Mouse: only one single SCG
//-------------------------------------------------------------------------
if( species == 'human' ) {

  def scgVal = params.scg.toString().trim().replace('SCG-','')  

  if( !(scgVal in ['1','2','3']) ) {
    error "Invalid --scg value given: ${params.scg}. Allowed values = [1,2,3] (human sample only)"
  }

  params.SCG_bed = "${refDir}/SCG-${scgVal}.bed"

} else {
  // mouse
  params.SCG_bed = "${refDir}/mSCG.bed"
}

//-------------------------------------------------------------------------
// rDNA contig
//-------------------------------------------------------------------------
params.rDNA_contig = params.rDNA_contig ?: (species == 'human' ? 'KY962518.1' : 'BK000964.3_p1')

//-------------------------------------------------------------------------
// Validate if reference, rDNA bed and SCG bed exist?
//-------------------------------------------------------------------------
['ref','rDNAgene_bed','SCG_bed'].each { k ->
  def p = params[k]
  if( !p || !file(p).exists() ) error "Missing ${k}: ${p}"
}

//-------------------------------------------------------------------------
// Logging
//-------------------------------------------------------------------------
log.info "[params] species=${species}"
log.info "[params] ref=${params.ref}"
log.info "[params] SCG_bed=${params.SCG_bed}"

// Hard-coded params 
// params.ref          = "/g/data/xc17/zaka/nextflow/RICO/ref/human/GRCh38_rDNAx5.fa"
// params.index_dir    = "/g/data/xc17/zaka/nextflow/RICO/ref/human/"
// params.rDNAgene_bed = "/g/data/xc17/zaka/nextflow/RICO/ref/human/rDNAgene.bed"
// params.rDNA_contig  = "KY962518.1"
// params.results_dir  = "/g/data/xc17/zaka/nextflow/RICO/results"
// params.samplesheet  = "/g/data/xc17/zaka/nextflow/RICO/samples.tsv"

//-------------------------------------------------------------------------
// Tools
//-------------------------------------------------------------------------
params.minimap2 = "/g/data/xc17/zaka/lib/minimap2-2.30/minimap2"
params.samtools = "/g/data/xc17/zaka/lib/samtools-1.23/bin/samtools"
params.bedtools = "/g/data/xc17/zaka/lib/bedtools-2.31.1/bin/bedtools"
params.modkit   = "/g/data/xc17/zaka/lib/modkit-0.6.1/modkit"

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
      echo "[\$(ts)] [index] FOUND ${exist_mmi} -> skipping" 1>&2
      ln -sf "${exist_mmi}" "${out_mmi}"
    else
      echo "[\$(ts)] [index] MISSING ${exist_mmi} -> building" 1>&2
      ${params.minimap2} --version 1>&2
      ${params.minimap2} -d "${out_mmi}" "${genome}"
      echo "[\$(ts)] [index] DONE -> ${out_mmi}" 1>&2
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
    tuple path("${unalign_bam.baseName}_${ref_mmi.baseName}.bam"),
          path("${unalign_bam.baseName}_${ref_mmi.baseName}.bam.bai")

  script:
    // Output naming:
    //   <inputBasename>_<refBasename>.bam
    // Example:
    //   hg005_1_test.bam + GRCh38_rDNAx5.mmi -> hg005_1_test_GRCh38_rDNAx5.bam
    def sample    = unalign_bam.baseName
    def refTag    = ref_mmi.baseName
    def out_bam   = "${sample}_${refTag}.bam"
    def out_bai   = "${out_bam}.bai"
    def exist_bam = "${params.results_dir}/${out_bam}"
    def exist_bai = "${params.results_dir}/${out_bai}"

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

    echo "[\$(ts)] [align] Start doing bam->fastq (MM,ML) | minimap2 | sort" 1>&2

    ${params.samtools} fastq "${unalign_bam}" -T "MM,ML" \
      | ${params.minimap2} -t ${task.cpus} -ax map-ont "${ref_mmi}" -Y -y --secondary=no - \
      | ${params.samtools} sort -@ ${task.cpus} -m 1G -o "${out_bam}"

    echo "[\$(ts)] [align] START samtools index ${out_bam}" 1>&2
    ${params.samtools} index "${out_bam}"

    echo "[\$(ts)] [align] DONE -> ${out_bam} (+.bai)" 1>&2
    """
}

// Mean-coverage approachh to calculate rDNA CN (rDNA_SUM / SCG_MEAN)
process mean_cn {
  tag { bam.baseName }

  publishDir "${params.results_dir}/mean_cn", mode: 'copy', overwrite: true

  input:
    tuple path(bam), path(bai)

  output:
    path "*.CN.mean.txt"
    path "*.45S.coverage.tsv"
    path "*_perBase.coverage"

  script:
    def bamBase  = bam.baseName
    def out_sum  = "${bamBase}.CN.mean.txt"
    def out_rDNA = "${bamBase}.45S.coverage.tsv"
    def scgName = file(params.SCG_bed).baseName   // SCG-1, SCG-2, SCG-3
    def out_scg = "${bamBase}_${scgName}_perBase.coverage"

    def exist_sum = "${params.results_dir}/mean_cn/${out_sum}"

    """
    set -euo pipefail
    ts(){ date +"%F %T"; }

    if [[ -s "${exist_sum}" ]]; then
      echo "[\$(ts)] [mean_cn] FOUND ${exist_sum} -> skipping" 1>&2
      ln -sf "${exist_sum}" "${out_sum}"

      if [[ -s "${params.results_dir}/mean_cn/${out_rDNA}" ]]; then
        ln -sf "${params.results_dir}/mean_cn/${out_rDNA}" "${out_rDNA}"
      else
        : > "${out_rDNA}"
      fi

      if [[ -s "${params.results_dir}/mean_cn/${out_scg}" ]]; then
        ln -sf "${params.results_dir}/mean_cn/${out_scg}" "${out_scg}"
      else
        : > "${out_scg}"
      fi

      exit 0
    fi

    echo "[\$(ts)] [mean_cn] Using ${bam}" 1>&2
    echo "[\$(ts)] [mean_cn] rDNA_contig=${params.rDNA_contig}" 1>&2
    echo "[\$(ts)] [mean_cn] rDNAgene_bed=${params.rDNAgene_bed}" 1>&2
    echo "[\$(ts)] [mean_cn] SCG_bed=${params.SCG_bed}" 1>&2

    ############################ 1) Filter rDNA-only reads
    RONLY_BAM="${bamBase}_rDNAonly.bam"

    echo "[\$(ts)] [mean_cn] Filtering rDNA-only reads (${params.rDNA_contig})..." 1>&2

    ${params.samtools} view -b "${bam}" "${params.rDNA_contig}" > "\${RONLY_BAM}"

    ############################ 2) rDNA_SUM: sum of mean coverages across mutiple repeats
    echo "[\$(ts)] [mean_cn] Computing rDNA_SUM across rDNAgene features..." 1>&2

    ${params.bedtools} bamtobed -i "\${RONLY_BAM}" \\
      | ${params.bedtools} coverage -a "${params.rDNAgene_bed}" -b - -mean > "${out_rDNA}"

    RDNA_SUM=\$(awk '{sum += \$NF} END {printf "%.2f", (sum+0)}' "${out_rDNA}")

    echo "[\$(ts)] [mean_cn] RDNA_SUM=\${RDNA_SUM}" 1>&2

    ############################ 3) SCG mean per-base coverage
    echo "[\$(ts)] [mean_cn] Computing per-base SCG depth..." 1>&2

    ${params.samtools} depth -b "${params.SCG_bed}" "${bam}" > "${out_scg}"

    SCG_MEAN=\$(awk '{sum+=\$3; n+=1} END {if(n>0) printf "%.2f", sum/n; else printf "0.00"}' "${out_scg}")

    echo "[\$(ts)] [mean_cn] SCG_MEAN=\${SCG_MEAN}" 1>&2

    ############################ 4) CN = RDNA_SUM / SCG_MEAN
    CN=\$(awk -v r="\${RDNA_SUM}" -v m="\${SCG_MEAN}" 'BEGIN{ if(m>0) printf "%.2f", r/m; else printf "NA" }')
 
    echo "[\$(ts)] [mean_cn] CN=\${CN}" 1>&2

    ############################ 5) Summary record
    SCG_realpath=\$(realpath "${params.SCG_bed}")

    {
      echo -e "bam\trDNA_contig\trDNAgene_bed\tSCG_bed\trDNA_SUM\tSCG_MEAN\tCN"
      echo -e "${bamBase}.bam\t${params.rDNA_contig}\t${params.rDNAgene_bed}\t\${SCG_realpath}\t\${RDNA_SUM}\t\${SCG_MEAN}\t\${CN}"
    } > "${out_sum}"

    echo "[\$(ts)] [mean_cn] DONE -> ${out_sum}" 1>&2
    """
}

process mod_pileup {
  tag { bam.baseName }

  publishDir "${params.mod_results}", mode: 'copy', overwrite: true

  input:
    tuple path(bam), path(bai)
    path genome

  output:
    path "*.CpG.bed"

  script:
    def bamBase   = bam.baseName
    def out_bed   = "${bamBase}.CpG.bed"
    def exist_out = "${params.mod_results}/${out_bed}"

    """
    set -euo pipefail
    ts(){ date +"%F %T"; }

    if [[ -s "${exist_out}" ]]; then
      echo "[\$(ts)] [mod_pileup] FOUND ${exist_out} -> skipping" 1>&2
      ln -sf "${exist_out}" "${out_bed}"
      exit 0
    fi

    echo "[\$(ts)] [mod_pileup] Start constructing a bedMethyl..." 1>&2
    echo "[\$(ts)] [mod_pileup] bam=${bam}" 1>&2
    echo "[\$(ts)] [mod_pileup] genome=${genome}" 1>&2
    echo "[\$(ts)] [mod_pileup] out=${out_bed}" 1>&2

    ${params.modkit} pileup "${bam}" "${out_bed}" \
      --ref "${genome}" \
      --modified-bases 5mC \
      --cpg \
      --combine-strands \
      --threads ${task.cpus}

    echo "[\$(ts)] [mod_pileup] DONE -> ${out_bed}" 1>&2
    """
}



process cg_motifs {
  tag { genome.baseName }

  publishDir "${params.index_dir}", mode: 'copy', overwrite: true

  input:
    path genome

  output:
    path "*.CG_motifs.bed"

  script:
    def refBase   = genome.baseName
    def out_bed   = "${refBase}.CG_motifs.bed"
    def exist_bed = "${params.index_dir}/${out_bed}"

    """
    set -euo pipefail
    ts(){ date +"%F %T"; }

    echo "[\$(ts)] [cg_motifs] genome=${genome}" 1>&2
    echo "[\$(ts)] [cg_motifs] out=${out_bed}" 1>&2

    if [[ -s "${exist_bed}" ]]; then
      echo "[\$(ts)] [cg_motifs] FOUND ${exist_bed} -> skipping" 1>&2
      ln -sf "${exist_bed}" "${out_bed}"
      exit 0
    fi

    echo "[\$(ts)] [cg_motifs] Extract only sites aligned to a CG motif..." 1>&2
    
    ${params.modkit} motif bed "${genome}" CG 0 > "${out_bed}"
    echo "[\$(ts)] [cg_motifs] DONE -> ${out_bed}" 1>&2
    """
}


process mod_extract_read {
  tag { bam.baseName }

  publishDir "${params.mod_results}", mode: 'copy', overwrite: true

  input:
    tuple path(bam), path(bai)
    path genome
    path cg_motifs_bed

  output:
    path "*.rDNA.calls.tsv"

  script:
    def bamBase   = bam.baseName
    def out_tsv   = "${bamBase}.rDNA.calls.tsv"
    def exist_out = "${params.mod_results}/${out_tsv}"

    """
    set -euo pipefail
    ts(){ date +"%F %T"; }

    if [[ -s "${exist_out}" ]]; then
      echo "[\$(ts)] [mod_extract_read] FOUND ${exist_out} -> skipping" 1>&2
      ln -sf "${exist_out}" "${out_tsv}"
      exit 0
    fi

    echo "[\$(ts)] [mod_extract_read] Start extracting per-read base modification from rDNA only..." 1>&2
    echo "[\$(ts)] [mod_extract_read] bam=${bam}" 1>&2
    echo "[\$(ts)] [mod_extract_read] genome=${genome}" 1>&2
    echo "[\$(ts)] [mod_extract_read] region=${params.rDNA_contig}" 1>&2
    echo "[\$(ts)] [mod_extract_read] include-bed=${cg_motifs_bed}" 1>&2
    echo "[\$(ts)] [mod_extract_read] out=${out_tsv}" 1>&2

    ${params.modkit} extract calls "${bam}" "${out_tsv}" \
      -t ${task.cpus} \
      --region "${params.rDNA_contig}" \
      --ref "${genome}" \
      --include-bed "${cg_motifs_bed}" \
      --mapped-only

    echo "[\$(ts)] [mod_extract_read] DONE -> ${out_tsv}" 1>&2
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

  aligned_ch = align(bam_ch, mmi_ch)
  mean_cn(aligned_ch)
  cg_bed_ch = cg_motifs(ref_ch)
  mod_pileup(aligned_ch, ref_ch)
  mod_extract_read(aligned_ch, ref_ch, cg_bed_ch)

}
