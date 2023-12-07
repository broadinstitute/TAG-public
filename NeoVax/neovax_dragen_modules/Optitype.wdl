version 1.0

import "./FormatUtil.wdl" as FormatUtil

struct Runtime {
  Int max_retries
  Int preemptible
  Int mem
  Int cpu
  Int boot_disk_size
  Int initial_disk_size
  String docker
}

workflow Optitype {
  input {
    String prefix
    File normal_bam
    File normal_bam_index
    File? hla_fasta_override
    
    # BWA parameters
    Int bwa_cores = 16
    Int bwa_extra_mem = 8
    Int bwa_mismatch_penalty = 4
    
    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 50
    Int boot_disk_size = 15
    Int mem = 16
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-optitype:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }

  call FormatUtil.SamToFastq as InputFastq {
    input:
      prefix = prefix,
      input_bam = normal_bam,
      max_retries = standard_runtime.max_retries,
      preemptible = standard_runtime.preemptible,
      mem = standard_runtime.mem,
      cpu = standard_runtime.cpu,
      boot_disk_size = standard_runtime.boot_disk_size,
      disk_pad = standard_runtime.initial_disk_size
  }
  
  call AlignReadsToHLA {
    input:
      prefix = prefix,
      hla_fasta_override = hla_fasta_override,
      fastq1 = InputFastq.fastq1,
      fastq2 = InputFastq.fastq2,
      mem = standard_runtime.mem + bwa_extra_mem,
      cpu = bwa_cores,
      mismatch_penalty = bwa_mismatch_penalty,
      docker = standard_runtime.docker,
      preemptible = 1,
      max_retries = standard_runtime.max_retries,
      additional_disk = standard_runtime.initial_disk_size * 3,
      boot_disk_size = standard_runtime.boot_disk_size
  }
  
  call FormatUtil.SamToFastq as HlaFastq {
    input:
      prefix = prefix,
      input_bam = AlignReadsToHLA.alignment,
      max_retries = standard_runtime.max_retries,
      preemptible = standard_runtime.preemptible,
      mem = standard_runtime.mem,
      cpu = standard_runtime.cpu,
      boot_disk_size = standard_runtime.boot_disk_size,
      disk_pad = standard_runtime.initial_disk_size
  }

  call RunOptitype {
    input:
      prefix = prefix,
      fastq1 = HlaFastq.fastq1,
      fastq2 = HlaFastq.fastq2,
      runtime_params = standard_runtime
  }

  output {
    File output_hla = RunOptitype.output_hla
    File output_raw = RunOptitype.output_raw
    File coverage_plot = RunOptitype.coverage_plot
  }
}


task AlignReadsToHLA {
  input {
    String prefix
    File? hla_fasta_override
    File fastq1
    File fastq2
    Int mismatch_penalty
    
    String docker
    Int preemptible
    Int max_retries
    Int additional_disk
    Int boot_disk_size
    Int mem
    Int cpu
  }

  Int fastq_size = ceil(size(fastq1, "GB") + size(fastq2, "GB"))
  Int disk_gb = fastq_size * 5 + additional_disk

  command <<<
    set -euxo pipefail
    
    HLA_REFERENCE=~{default="/root/OptiType-1.3.2/data/hla_reference_dna.fasta" hla_fasta_override}
    
    FASTQ1_INPUT=$(echo ~{fastq1} | rev | cut -f 2- -d '.' | rev)
    FASTQ2_INPUT=$(echo ~{fastq2} | rev | cut -f 2- -d '.' | rev)

    gunzip ~{fastq1}
    gunzip ~{fastq2}

    /root/bwa-0.7.17/bwa index $HLA_REFERENCE
    /root/bwa-0.7.17/bwa mem \
      -t ~{cpu} \
      -B ~{mismatch_penalty} \
      $HLA_REFERENCE $FASTQ1_INPUT $FASTQ2_INPUT | \
    samtools sort -@~{cpu} -O BAM -o ~{prefix}.aligned.sorted.bam

    # Only select mapped reads
    samtools view -F 4 -b ~{prefix}.aligned.sorted.bam > ~{prefix}.aligned.filtered.bam
  >>>
  
  output {    
    File alignment = "~{prefix}.aligned.filtered.bam"
  }
  
  runtime {
    docker        : docker
    bootDiskSizeGb: boot_disk_size
    preemptible   : preemptible
    cpu           : cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : mem + "GB"
    maxRetries    : max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}

task RunOptitype {
  input {
    String prefix
    File fastq1
    File fastq2
    Runtime runtime_params
  }
  
  Int disk_gb = ceil((size(fastq1,"GB") + size(fastq2,"GB")) * 2.5) + runtime_params.initial_disk_size

  command <<<
    set -euxo pipefail
    
    python /root/OptiType-1.3.2/OptiTypePipeline.py \
      -i ~{fastq1} ~{fastq2} --prefix ~{prefix} --dna -v -o .

    # Change output names
    mv ~{prefix}_result.tsv ~{prefix}.optitype.output.txt
    mv ~{prefix}_coverage_plot.pdf ~{prefix}.optitype.coverage_plot.pdf
    
    # Format predicted HLA alleles
    python <<CODE
    fo = open('~{prefix}.optitype.hla.txt','w')
    for line in open('~{prefix}.optitype.output.txt'):
        if line[0] != '0':
            continue
        for allele in line.replace('*','').split('\t'):
            if allele[0] in ['A', 'B', 'C']:
                fo.write('HLA-'+allele+'\n')
    fo.close()
    CODE
  >>>

  output {
    File output_hla = "~{prefix}.optitype.hla.txt"
    File output_raw = "~{prefix}.optitype.output.txt"
    File coverage_plot = "~{prefix}.optitype.coverage_plot.pdf"
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
