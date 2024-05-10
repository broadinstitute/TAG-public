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

workflow Polysolver {
  input {
    String prefix
    File normal_bam
    File normal_bam_index
    String? race

    # Reference genome version
    # option: 'hg18', 'hg19' or 'hg38'
    String reference_version

    # FASTQ quality score encoding format option for Novoalign
    # options:
    #   STDFQ  = Sanger coding: -10log10(Perr) + '!'
    #   ILMFQ  = Illumina coding: -10log10(Perr) + '@'
    #   ILM1.8 = Illumina Casava V1.8 with Sanger coding: -10log10(Perr) + '!'
    #   SLXFQ  = Solexa coding: 10log10(P/(1-P)) + '@'
    String fastq_format = "STDFQ"
    
    # whether population-level allele frequencies should be used as priors
    Boolean include_freq = true
    
    # whether empirical insert size distribution should be used in the model
    Boolean insert_calc = false

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 20
    Int boot_disk_size = 15
    Int mem = 8
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-polysolver:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }
  
  call RunPolysolver {
    input:
      prefix = prefix,
      normal_bam = normal_bam,
      normal_bam_index = normal_bam_index,
      race = race,
      reference_version = reference_version,
      fastq_format = fastq_format,
      include_freq = include_freq,
      insert_calc = insert_calc,
      runtime_params = standard_runtime
  }

  output {
    File output_hla = RunPolysolver.output_hla
    File output_raw = RunPolysolver.output_raw
  }
}


task RunPolysolver {
  input {
    String prefix
    File normal_bam
    File normal_bam_index
    String reference_version
    String? race
    Boolean include_freq
    Boolean insert_calc
    String fastq_format

    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size + ceil(size(normal_bam,"GB") + size(normal_bam_index,"GB"))

  command <<<
    set -e

    EFFECTIVE_POLYSOLVER_RACE="Unknown"
    if [[ -n "~{race}" ]] ; then
       EFFECTIVE_POLYSOLVER_RACE="~{race}"
    fi
    echo "set EFFECTIVE_POLYSOLVER_RACE = $EFFECTIVE_POLYSOLVER_RACE"

    # create output dir
    mkdir -pv hla_out

    python /home/process_monitor.py \
      /bin/bash /home/polysolver/scripts/shell_call_hla_type \
      ~{normal_bam} \
      $EFFECTIVE_POLYSOLVER_RACE \
      ~{true="1" false="0" include_freq} \
      ~{reference_version} \
      ~{fastq_format} \
      ~{true="1" false="0" insert_calc} \
      $PWD/hla_out

    # rename output file
    mv hla_out/winners.hla.txt ~{prefix}.polysolver.output.txt
    
    # format output file
    python <<CODE
    fo = open('~{prefix}.polysolver.hla.txt','w')
    for line in open('~{prefix}.polysolver.output.txt'):
        line = line.rstrip().split('\t')
        for allele in line[1:]:
            fo.write(line[0]+':'.join(allele.split('_')[2:4])+'\n')
    fo.close()
    CODE
  >>>

  output {
    File output_hla = "~{prefix}.polysolver.hla.txt"
    File output_raw = "~{prefix}.polysolver.output.txt"
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
}
