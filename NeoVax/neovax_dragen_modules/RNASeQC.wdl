version 1.0

struct Runtime {
  Int max_retries
  Int preemptible
  Int mem
  Int cpu
  Int boot_disk_size
  Int initial_disk_size
  String docker
}

workflow RNASeQC {
  input {
    String prefix
    File input_bam
    File genes_gtf

    # RNASeQC parameters
    # Strand-specific metrics: 'FR' or 'RF' (default: off)
    String? strandedness
    # BED file containing non-overlapping exons used for fragment size calculations
    File? intervals_bed
    # Number of counts on a gene to consider the gene is detected
    Int detection_threshold = 5
    # Maximum number of allowed mismatches between a read and the reference sequence
    Int base_mismatch = 6
    
    String? rnaseqc_extra_args

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 20
    Int boot_disk_size = 15
    Int mem = 10
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-rnaseq:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }

  call RunRNASeQC {
    input:
      prefix = prefix,
      input_bam = input_bam,
      genes_gtf = genes_gtf,
      strandedness = strandedness,
      intervals_bed = intervals_bed,
      base_mismatch = base_mismatch,
      detection_threshold = detection_threshold,
      rnaseqc_extra_args = rnaseqc_extra_args,
      runtime_params = standard_runtime
  }

  output {
    File gene_tpm = RunRNASeQC.gene_tpm
    File gene_counts = RunRNASeQC.gene_counts
    File exon_counts = RunRNASeQC.exon_counts
    File metrics = RunRNASeQC.metrics
    File fragment_size = RunRNASeQC.fragment_size
    File gene_duplicates = RunRNASeQC.gene_duplicates
  }
}

task RunRNASeQC {
  input {
    String prefix
    File input_bam
    File genes_gtf
    Int base_mismatch
    Int detection_threshold
    String? strandedness
    File? intervals_bed
    String? rnaseqc_extra_args
    Runtime runtime_params
  }
  
  Int disk_gb = runtime_params.initial_disk_size + ceil(size(input_bam,"GB")+size(genes_gtf,"GB"))
                   
  command {
    set -euo pipefail

    touch ~{prefix}.fragmentSizes.txt
    touch ~{prefix}.gene_duplicates.gct
    
    rnaseqc ~{genes_gtf} ~{input_bam} . \
      --sample ~{prefix} \
      --detection-threshold ~{detection_threshold} \
      --base-mismatch ~{base_mismatch} \
      ~{"--bed " + intervals_bed} \
      ~{"--stranded " + strandedness} \
      ~{rnaseqc_extra_args}

    gzip *.gct
  }

  output {
    File gene_tpm = "~{prefix}.gene_tpm.gct.gz"
    File gene_counts = "~{prefix}.gene_reads.gct.gz"
    File exon_counts = "~{prefix}.exon_reads.gct.gz"
    File metrics = "~{prefix}.metrics.tsv"
    File fragment_size = "~{prefix}.fragmentSizes.txt"
    File gene_duplicates = "~{prefix}.gene_duplicates.gct.gz"
  }

  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}
