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

workflow IndexReferenceRNA {
  input {
    String star_prefix
    String rsem_prefix
    File ref_fasta
    File annotation_gtf

    # Parameter (sequencing read length - 1) for STAR index.
    # The index should be built to match the sequencing read length,
    # specified by the sjdbOverhang parameter. For example,
    # 2x76 bp paired-end sequencing protocol, this will be 75.
    Int overhang

    # workflow level runtime parameters
    Int preemptible = 1
    Int max_retries = 1
    Int additional_disk = 80
    Int boot_disk_size = 15
    Int mem = 52
    Int cpu = 4
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-rnaseq:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }
  call StarIndex {
    input:
      prefix = star_prefix + "_oh" + overhang,
      ref_fasta = ref_fasta,
      annotation_gtf = annotation_gtf,
      overhang = overhang,
      runtime_params = standard_runtime
  }
  call RsemReference {
    input:
      prefix = rsem_prefix,
      ref_fasta = ref_fasta,
      annotation_gtf = annotation_gtf,
      runtime_params = standard_runtime
  }
  
  output {
    File star_index = StarIndex.star_index
    File rsem_reference = RsemReference.rsem_reference
  }
}

task StarIndex {
  input {
    String prefix
    File ref_fasta
    File annotation_gtf
    Int overhang

    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size +
                ceil(size(ref_fasta,"GB") + size(annotation_gtf,"GB")) * 2

  command {
    set -e
    
    mkdir ~{prefix}
    
    STAR --runMode genomeGenerate \
         --genomeDir ~{prefix} \
         --genomeFastaFiles ~{ref_fasta} \
         --sjdbGTFfile ~{annotation_gtf} \
         --sjdbOverhang ~{overhang} \
         --runThreadN ~{runtime_params.cpu}
          
    tar -cvzf ~{prefix}.tar.gz ~{prefix}
  }

  output {
    File star_index = "~{prefix}.tar.gz"
  }

  runtime {
    docker     : runtime_params.docker
    memory     : runtime_params.mem + "GB"
    disks      : "local-disk " + disk_gb + " HDD"
    cpu        : runtime_params.cpu 
    preemptible: runtime_params.preemptible
    maxRetries : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}

task RsemReference {
  input {
    String prefix
    File ref_fasta
    File annotation_gtf
    
    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size +
                ceil(size(ref_fasta,"GB") + size(annotation_gtf,"GB")) * 2  

  command {
    set -e
    mkdir ~{prefix} && cd ~{prefix}

    # set prefix of the reference as 'rsem_reference'
    rsem-prepare-reference ~{ref_fasta} rsem_reference \
      --gtf ~{annotation_gtf} --num-threads ~{runtime_params.cpu}

    cd .. && tar -cvzf ~{prefix}.tar.gz ~{prefix}
  }

  output {
    File rsem_reference = "~{prefix}.tar.gz"
  }
  
  runtime {
    docker     : runtime_params.docker
    memory     : runtime_params.mem + "GB"
    disks      : "local-disk " + disk_gb + " HDD"
    cpu        : runtime_params.cpu 
    preemptible: runtime_params.preemptible
    maxRetries : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}
