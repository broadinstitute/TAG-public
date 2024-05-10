version 1.0

# formatting related tasks with hard-coded public docker images
workflow FormatUtil {}

task GetSampleID {
  input {
    File bam
    File bam_index

    Int max_retries = 1
    Int preemptible = 2
    Int mem = 4
    Int cpu = 1
    Int boot_disk_size = 15
    Int disk_pad = 50
  }

  Int disk_gb = disk_pad + ceil(size(bam,"GB") + size(bam_index,"GB"))
  
  command <<<
    set -e

    samtools view -H ~{bam} | \
    grep ^@RG | tr "\t" "\n" | \
    grep ^SM: | sort -u | cut -f2 -d":" > sample_id.txt

    # exit with error if there are more than one SM tag
    if [[ `cat sample_id.txt | wc -l` -gt 1 ]] ; then
        echo "error: more than one SM tag!" 1>&2
        exit 1
    fi
  >>>

  output {
    String sample_id = read_string("sample_id.txt")
  }

  runtime {
    docker        : "us.gcr.io/tag-public/neovax-tag-gatk:v1"
    bootDiskSizeGb: boot_disk_size
    preemptible   : preemptible
    cpu           : cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : mem + "GB"
    maxRetries    : max_retries
  }

  meta {
    author: "Junko Tsuji"
    email: "jtsuji@broadinstitute.org"
  }
}

task SamToFastq {
  input {
    File input_bam
    String prefix

    Int max_retries = 1
    Int preemptible = 2
    Int mem = 16
    Int cpu = 2
    Int boot_disk_size = 15
    Int disk_pad = 50
  }
  
  Int disk_gb = disk_pad + ceil(size(input_bam,"GB")*2.5)
  Int compute_mem = mem - 1

  command <<<
    set -euxo pipefail
    
    gatk --java-options "-Xmx~{compute_mem}g" SamToFastq \
        --INPUT ~{input_bam} \
        --FASTQ ~{prefix}_1.fastq.gz \
        --SECOND_END_FASTQ ~{prefix}_2.fastq.gz \
        --UNPAIRED_FASTQ ~{prefix}_unpaired.fastq.gz \
        --INTERLEAVE false \
        --INCLUDE_NON_PF_READS true \
        --INCLUDE_NON_PRIMARY_ALIGNMENTS false \
        --VALIDATION_STRINGENCY SILENT
  >>>

  output {
    File fastq1 = "~{prefix}_1.fastq.gz"
    File fastq2 = "~{prefix}_2.fastq.gz"
  }

  runtime {
    docker        : "us.gcr.io/tag-public/neovax-tag-gatk:v1"
    bootDiskSizeGb: boot_disk_size
    preemptible   : preemptible
    cpu           : cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : mem + "GB"
    maxRetries    : max_retries
  }

  meta {
    author: "Junko Tsuji"
    email: "jtsuji@broadinstitute.org"
  }
}

task FastqToSam {
  input {
    File fastq1
    File fastq2
    String prefix

    Int max_retries = 1
    Int preemptible = 2
    Int mem = 16
    Int cpu = 2
    Int boot_disk_size = 15
    Int disk_pad = 50
  }
  
  Int disk_gb = disk_pad + ceil((size(fastq1,"GB")+size(fastq2,"GB"))*2.5)
  Int compute_mem = mem - 1

  command <<<
    set -euxo pipefail
    
    gatk --java-options "-Xmx~{compute_mem}g" FastqToSam \
        --FASTQ ~{fastq1} \
        --FASTQ2 ~{fastq2} \
        --SAMPLE_NAME ~{prefix} \
        --LIBRARY_NAME ~{prefix} \
        --PLATFORM "ILLUMINA" \
        --READ_GROUP_NAME "RG_~{prefix}" \
        --PLATFORM_UNIT "barcode_~{prefix}" \
        --OUTPUT ~{prefix}.u.bam \
        --VALIDATION_STRINGENCY SILENT
  >>>

  output {
    File unmapped_bam = "~{prefix}.u.bam"
  }

  runtime {
    docker        : "us.gcr.io/tag-public/neovax-tag-gatk:v1"
    bootDiskSizeGb: boot_disk_size
    preemptible   : preemptible
    cpu           : cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : mem + "GB"
    maxRetries    : max_retries
  }

  meta {
    author: "Junko Tsuji"
    email: "jtsuji@broadinstitute.org"
  }
}

task AggregateMAF {
  input {
    String output_prefix
    Array[File] input_mafs
    
    Boolean strip_comments = false
    Boolean sort_by_chromosome = true
    String comment_char = "#"

    Int max_retries = 1
    Int preemptible = 2
    Int mem = 4
    Int cpu = 1
    Int boot_disk_size = 15
    Int disk_pad = 20
  }

  Int disk_gb = disk_pad + ceil(size(input_mafs,"GB"))

  command {
    set -euxo pipefail

    python /root/scripts/tsvConcatFiles.py \
      --outputFilename ~{output_prefix}.maf \
      ~{sep=" " input_mafs}
    
    # exclude any comments from aggregated MAF
    if [[ ~{strip_comments} = true ]] ; then
      grep -v ^"~{comment_char}" ~{output_prefix}.maf > ~{output_prefix}.stripped.maf
      mv ~{output_prefix}.stripped.maf ~{output_prefix}.maf
    fi

    # sort mutations by genomic coordinates
    if [[ ~{sort_by_chromosome} = true ]] ; then
      cat <(grep ^"~{comment_char}" ~{output_prefix}.maf) \
          <(grep ^Hugo_Symbol ~{output_prefix}.maf) \
          <(grep -v ^"~{comment_char}" ~{output_prefix}.maf | grep -v ^Hugo_Symbol | sort -k5,7 -V) > ~{output_prefix}.sorted.maf
      mv ~{output_prefix}.sorted.maf ~{output_prefix}.maf
    fi
  }

  output {
    File output_maf = "~{output_prefix}.maf"
  }

  runtime {
    docker        : "us.gcr.io/tag-public/neovax-tag-tumorqc:v1"
    bootDiskSizeGb: boot_disk_size
    preemptible   : preemptible
    cpu           : cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : mem + "GB"
    maxRetries    : max_retries
  }

  meta {
    author: "Junko Tsuji"
    email: "jtsuji@broadinstitute.org"
  }
}
