# rDNA-CN-pipeline

A long-read sequencing–based pipeline to estimate **ribosomal DNA (rDNA) copy number (CN)** using coverage normalization against single-copy genes (SCGs)

The pipeline is designed for **Oxford Nanopore Technologies (ONT)** whole-genome sequencing data and supports both **CN estimation** and **DNA methylation analysis**

## Overview

This pipeline processes alignment, coverage quantification, and rDNA CN estimation. rDNA CN is calculated by normalizing rDNA coverage against genome-wide SCG coverage.
