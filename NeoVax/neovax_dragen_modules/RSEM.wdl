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

workflow RSEM {
  input {
    String prefix  
    File transcriptome_bam
    File rsem_reference

    # RSEM options
    # maximum fragment length
    Int max_frag_len = 1000
    # input reads are paired
    Boolean is_paired_end = true
    # estimate the read start position distribution from data
    Boolean estimate_rspd = true
    # calculate 95% credibility intervals and posterior mean estimates
    Boolean calculate_ci = false
    # choices for the strandedness option
    #   'none'   : Non-strand-specific protocols (default)
    #   'forward': All (upstream) reads are derived from the forward strand
    #   'reverse': All (upstream) reads are derived from the reverse strand
    #              For Illumina TruSeq Stranded protocols, use 'reverse'.
    String strandedness = "none"
    String? rsem_extra_args

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 50
    Int boot_disk_size = 15
    Int mem = 24
    Int cpu = 4
  }
  
  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-rnaseq:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }
  call RunRSEM {
    input:
      prefix = prefix,
      transcriptome_bam = transcriptome_bam,
      rsem_reference = rsem_reference,
      max_frag_len = max_frag_len,
      estimate_rspd = estimate_rspd,
      strandedness = strandedness,
      is_paired_end = is_paired_end,
      calculate_ci = calculate_ci,      
      rsem_extra_args = rsem_extra_args,
      runtime_params = standard_runtime
  }
  output {
    File genes_expr = RunRSEM.genes
    File isoforms_expr = RunRSEM.isoforms
  }
}

task RunRSEM {
  input {
    String prefix
    File transcriptome_bam
    File rsem_reference
    
    Int max_frag_len
    Boolean estimate_rspd
    Boolean is_paired_end
    String strandedness
    Boolean calculate_ci
    String? rsem_extra_args

    Runtime runtime_params
  }
  
  Int disk_gb = runtime_params.initial_disk_size +
                ceil((size(transcriptome_bam,"GB")+size(rsem_reference,"GB")) * 2.5)

  command {
    set -euo pipefail
    mkdir rsem_reference
    tar -xvvf ~{rsem_reference} -C rsem_reference --strip-components=1
    
    REF_PREFIX=`basename rsem_reference/*.grp .grp`
    if [[ -z "$REF_PREFIX" ]] ; then
      echo "Can't find RSEM reference in ~{rsem_reference}"
      exit 1
    fi

    rsem-calculate-expression \
      --no-bam-output \
      --num-threads ~{runtime_params.cpu} \
      --fragment-length-max ~{max_frag_len} \
      --strandedness ~{strandedness} \
      ~{true="--estimate-rspd" false="" estimate_rspd} \
      ~{true="--paired-end" false="" is_paired_end} \
      ~{true="--calc-ci" false="" calculate_ci} \
      ~{rsem_extra_args} \
      --bam ~{transcriptome_bam} \
      rsem_reference/$REF_PREFIX \
      "~{prefix}.rsem"

    gzip *.results
  }

  output {
    File genes = "~{prefix}.rsem.genes.results.gz"
    File isoforms = "~{prefix}.rsem.isoforms.results.gz"
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
