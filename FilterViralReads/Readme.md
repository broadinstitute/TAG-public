# Viral BAM Filter Pipeline

This repository contains a WDL (Workflow Description Language) pipeline designed to extract 
viral-specific reads from a sample BAM file using `BBDuk`. It identifies reads (and their mates) matching a 
provided list of k-mers unique to the viral reference and outputs a subsetted BAM file 
containing only those reads.

The pipeline is optimized for use on the **Terra** platform and can be found in the workspace:
[tag_2460_lentiviral_insertion](https://app.terra.bio/#workspaces/broadtagteam/tag_2460_lentiviral_insertion)

---

## Overview

The workflow filters bam files to isolate viral sequences (e.g., Lentiviral) by converting 
alignments to FASTQ, performing k-mer matching against a reference, and re-subsetting the original BAM 
to contain only reads that have k-mers from the virus. 

### Workflow Steps

1.  **Count Input Reads**: Calculates the initial number of aligned reads in the source BAM. 

2.  **Conversion**: Converts the input BAM to FASTQ format using `samtools`. 

3.  **Viral Filtering**: Uses `bbduk.sh` to match reads against the provided viral reference FASTA using k-mer matching (k = 19, hdist = 0). 

4.  **Name Extraction**: Extracts and cleans unique read names that matched the viral reference. 

5.  **BAM Subsetting**: Uses `samtools view -N` to create a new BAM containing only the identified viral reads. 

---

## File Structure

| File | Description |
| --- | --- |
| `FilterViralReads.wdl` | The main WDL script defining the workflow and task logic. 
| `Dockerfile` | Defines the environment (Java 11, Samtools, BBTools) for the pipeline. 
| `docker_build.sh` | Utility script to build and push the image to Docker Hub (`fleharty/viral-bam-filter:v2`). |
| `FilterViralReads.inputs.json` | Example JSON template for workflow parameters. |
| `FilterViralReads.outputs.json` | Example JSON template for workflow parameters. |

---

## Inputs and Outputs

### Key Inputs

* **input_bam**: The original BAM file to be filtered.
* **viral_reference**: A FASTA file containing the viral k-mers to filter against.  It's important that these k-mers are unique to the viral sequence, and not found in the reference of the sample.
* **sample_name**: String used for naming output files. 
* **threads/memory_gb**: Resource allocations (default 8 threads, 16 GB RAM). 

### Primary Outputs

* **viral_bam**: A BAM file containing only the reads matching the viral reference.
* **viral_aligned_read_count**: Integer count of aligned reads in the final viral BAM. 
* **bbduk_stats**: A diagnostic text file showing the filtering statistics from BBMap. 

---

## Setup and Execution

### Docker Image

The pipeline relies on a custom Docker image. If you need to rebuild it:

1. Ensure you have Docker installed and authenticated.
2. Run the build script:

```bash
./docker_build.sh

```

### Running on Terra

The Terra Workspace is available at:
https://app.terra.bio/#workspaces/broadtagteam/tag_2460_lentiviral_insertion

