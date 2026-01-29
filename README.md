# rDNA-CN-pipeline

A long-read sequencing–based pipeline to estimate **ribosomal DNA (rDNA) copy number (CN)** using coverage normalization against single-copy genes (SCGs)

The pipeline is designed for **Oxford Nanopore Technologies (ONT)** whole-genome sequencing data and supports both **CN estimation** and **DNA methylation analysis**

![workflow](images/rDNA-CN-workflow.png)

## Overview

This pipeline processes alignment, coverage quantification, and rDNA CN estimation. rDNA CN is calculated by normalizing rDNA coverage against genome-wide SCG coverage.

## Assumptions and current scope

This pipeline is currently designed and tested only on **NCI Australia (Gadi)** using PBSPro as the job scheduler.

Absolute paths are used for the tools (Minimap2, Samtools) used in the pipeline to ensure reproducibility on HPC.

## Get the pipeline

Clone the repository on Gadi:
```
git clone https://github.com/comprna/rDNA-CN-pipeline.git
cd rDNA-CN-pipeline
```

## Set up

Load Nextflow on Gadi:
```
module load nextflow/25.04.6
```

Check Nextflow is loaded:
```
nextflow -version
```

Configure project and storage

Edit `nextflow.config` and update the following fields to match your own NCI project:
```bash
process {
  executor = 'pbspro'
  project  = 'jd21' # replace with your own project
  storage  = 'gdata/jd21+scratch/jd21+gdata/xc17' # replace with your own project
}
```
* project: your NCI project code
* storage: filesystems where you saved the pipeline and access the files from

## Pipeline inputs

### Reference genome

Only GRCh38_rDNAx5.fa is available:
```
params.ref = "/g/data/xc17/zaka/res/human/rDNAx5/GRCh38_rDNAx5.fa"
```
***A portable reference distribution will be provided in a future update so users can download***

### Unaligned BAM

Provide the file path of your input BAM(s) in the `sample.tsv`. Each row corresponds to one unaligned BAM, and samples will be processed in parallel.

## Run the pipeline

From the pipeline directory, run:
```
nextflow run map.nf -config nextflow.config --samplesheet samples.tsv
```

## Outputs

For each input BAM, the pipeline produces:
* `<inputBasename>`_GRCh38_rDNAx5.bam
* `<inputBasename>`_GRCh38_rDNAx5.bam.bai




