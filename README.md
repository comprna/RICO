# RICO

RICO (**Ri**bosomal DNA **Co**py number) is long-read sequencing–based pipeline for estimating **ribosomal DNA (rDNA) copy number (CN)** using coverage normalization against **single-copy genes (SCGs)**, with optional **DNA methylation analysis**.

The pipeline is designed for **Oxford Nanopore Technologies (ONT)** whole-genome sequencing data and supports both **rDNA CN estimation** and **CpG methylation profiling**.

![workflow](images/rDNA-CN-workflow.png)

---

## Overview

This pipeline performs the following major steps:

1. Alignment of unaligned ONT BAMs to a custom reference genome containing rDNA repeats  
2. Coverage quantification over rDNA loci and genome-wide SCGs  
3. Estimation of rDNA copy number by normalizing rDNA coverage against SCG coverage  
4. DNA methylation analysis at CpG sites within rDNA regions

The workflow is implemented in Nextflow (DSL2) and is optimised for high-performance computing (HPC) environments.

---

## Assumptions and current scope

- The pipeline is currently **designed and tested only on NCI Australia (Gadi)**.
- Job submission uses **PBSPro**.
- Tool paths are specified as absolute paths to ensure reproducibility.
- The reference genome is currently fixed to a custom rDNA-augmented GRCh38 build.

## Tools used in the pipeline
* Nextflow (25.04.6)
* Minimap2 (2.30)
* Bedtools (2.31.1)
* Modkit (0.6.1)

---

## Get the pipeline

Clone the repository on Gadi:

```
git clone https://github.com/comprna/rDNA-CN-pipeline.git
cd rDNA-CN-pipeline
```

## Set up

Load Nextflow on Gadi
```
module load nextflow/25.04.6
```

Check that Nextflow is available:
```
nextflow -version
```
## Configure project and storage

Edit `nextflow.config` and update the following fields to match your NCI project:
```bash
process {
  executor = 'pbspro'
  project  = 'jd21' # replace with your NCI project
  storage  = 'gdata/jd21+scratch/jd21+gdata/xc17' # replace with your NCI project
}
```
* project: your NCI project code
* storage: file systems used by the pipeline (where data and results live)

## Inputs

### Reference genome

Currently supported reference:
```
params.ref = "/g/data/xc17/zaka/res/human/rDNAx5/GRCh38_rDNAx5.fa"
```
This reference contains 5 copies of the human rDNA repeat integrated into GRCh38.

### Unaligned BAMs

Input data are provided via a tab-separated sample sheet (samples.tsv). Each row corresponds to 1 unaligned BAM file, so samples can be processed in parallel.

Example `samples.tsv`:
```bash
file_path
/g/data/xc17/zaka/nextflow/rDNA-CN-pipeline/unaligned_bams/test.bam
```

## Run the pipeline

From the pipeline directory:
```
nextflow run map.nf -config nextflow.config --samplesheet samples.tsv
```

## Outputs

For each input BAM, the pipeline produces:

**1. Alignment outputs**
* `<sample>_GRCh38_rDNAx5.bam`
* `<sample>_GRCh38_rDNAx5.bam.bai`

**2. Copy number estimation**
* Coverage outputs for rDNA and SCGs
* rDNA CN estimates

**3. Methylation outputs**
* BedMethyl file
* Per-read base modification output table for rDNA





