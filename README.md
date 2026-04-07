# RICO

RICO (**Ri**bosomal DNA **Co**py number) is a long-read sequencing–based pipeline for estimating **ribosomal DNA (rDNA) copy number (CN)** using coverage normalization against **single-copy genes (SCGs)**, with **DNA methylation analysis**.

The pipeline is designed for **Oxford Nanopore Technologies (ONT)** whole-genome sequencing data and supports both **rDNA CN estimation** and **CpG methylation profiling**.

![workflow](images/rico_workflow.png)


## Table of contents
* [Overview](#overview)
* [Installation](#installation)
* [How to Use](#how-to-use)

---

## Overview

RICO performs the following major steps:

1. Alignment of unaligned ONT BAMs to a custom reference genome containing rDNA repeats  
2. Coverage quantification over rDNA loci and curated SCG panel  
3. Estimation of rDNA CN by normalizing rDNA coverage against SCG coverage  
4. Extraction of CpG methylation calls within rDNA regions

The workflow is implemented in *Nextflow (DSL2)* and is optimised for *HPC* environments.

### Execution environment

RICO is implemented in Nextflow (DSL2).
It is developed and tested on a PBSPro-based HPC system ([NCI Australia](https://nci.org.au/)).

**UPDATE (April 2026):** 
- **RICO now supports execution using Singularity containers**
- **This means all required software dependencies are provided**

The pipeline currently assumes:
- A PBSPro scheduler (executor = 'pbspro')
- HPC-style job submission
- Singularity available on the system (`module load singularity` on Gadi)

Other PBSPro-based HPC systems may work with minimal configuration changes.
Support for `Slurm` or other schedulers has not yet been tested.

### Tools used in the pipeline

RICO was developed and tested with the following versions (bundled in the container image):
* Nextflow (25.04.6)
* Minimap2 (2.30)
* Samtools (1.23)
* Bedtools (2.31.1)
* Modkit (0.6.1)

> When using the Singularity container, manual installation of these tools is no longer required.

---

## Installation

### 1. Clone the repository on Gadi:

```
git clone https://github.com/comprna/RICO.git
cd RICO
```

### 2. Load Nextflow and singularity on Gadi

```bash
module load nextflow/25.04.6
module load singularity
nextflow -version #check
```

### 3. Pull the Singularity container

```bash
mkdir -p /path/to/containers
singularity pull /path/to/containers/rico_2026.03.sif docker://zakayuen/rico:2026.03
```
> Replace `/path/to/containers` with your preferred location (e.g. /g/data/.../containers)

### 4. Download reference files

Human and mouse reference genomes and annotation files are available on Zenodo:

DOI: https://doi.org/10.5281/zenodo.18500657

Download the archive and extract:
```
tar -xzf RICO_ref_v1.tar.gz
```

> Place the extracted `ref/` directory inside the cloned RICO repository

RICO currently supports
- Human (GRCh38 + rDNAx5)
- Mouse (GRCm39 + rDNAx5)

### 5. Configure NCI project, storage, and mirror path

Edit `nextflow.config` and update the following fields:
```bash
params {
  project   = 'jd21'
  storage   = 'gdata/jd21+scratch/jd21+gdata/xc17+gdata/qq78'
  container = '/g/data/xc17/zaka/containers/rico_2026.03.sif'
...
```
* project: your NCI project code
* storage: file systems used by the pipeline (where input data and results are stored)
* container: /path/to/containers/rico_2026.03.sif

---

## How to use

### Input - unaligned BAM

Path to your input BAM needs to be provided in `samples.tsv`:

```bash
file_path
/path/to/your.bam 
```

### Minimum command (default: human, SCG-2)

From the RICO directory:
```
nextflow run rico.nf \
 --samplesheet samples.tsv
```

The results are written to the `results` folder by default.

But you can specify the results directory using `--out_dir`:
```
nextflow run rico.nf \
 --samplesheet samples.tsv \
 --out_dir /path/to/output
```

### Other SCGs for human samples

For human samples, three curated SCG panels are provided: `SCG-1`, `SCG-2` (default), `SCG-3`
Specify using `--scg`, for example:
```
nextflow run rico.nf \
 --samplesheet samples.tsv \
 --scg 3
```
If not specified, `SCG-2` is used.

### Mouse samples

```
nextflow run rico.nf \
 --samplesheet samples.tsv \
 --species mouse
```

For mouse samples, it uses a single curated SCG panel only. The `--scg` parameter is ignored when `--species mouse` is selected.

## Outputs

For each input BAM, RICO produces:

**1. Alignment outputs**
* `<sample>_<reference>.bam`
* `<sample>_<reference>.bam.bai`

**2. Copy number estimation**
* Coverage outputs for rDNA and SCGs
* rDNA CN estimates

**3. Methylation outputs**
* BedMethyl file
* Per-read base modification output table for rDNA

## Citation

> Yuen, Leeder, Hannan, Eyras & Hein. Accurate estimation of ribosomal DNA copy number using nanopore long-read sequencing. (Manuscript in preparation)



