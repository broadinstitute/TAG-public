version 1.0

##First version?
##Copyright Broad Institute, 2019
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One VCF file and its index
##
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.


# WORKFLOW DEFINITION
workflow HaplotypeCaller_vcf_gatk4 {
  input {
    File input_bam
    File input_bam_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    File intervals

    Int scatter_count
    String? split_intervals_extra_args

    Boolean make_vcf = true
    Boolean make_bamout = false
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    String gatk_path = "/gatk/gatk"
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String samtools_path = "samtools"


    Int? additional_disk
    Int disk_pad = select_first([additional_disk, 10])
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB"))

    Int? preemptible_tries
    Int preemptible_attempts = select_first([preemptible_tries, 3])

    Int boot_disk_size = 12
    Int small_task_mem = 4
    Int small_task_disk = 100
    Int? max_retries
    Int small_task_cpu = 2

  }

    Array[File] scattered_calling_intervals = read_lines(intervals)

    #is the input a cram file?
    Boolean is_cram = sub(basename(input_bam), ".*\\.", "") == "cram"

    String sample_basename = if is_cram then  basename(input_bam, ".cram") else basename(input_bam, ".bam")
    String vcf_basename = sample_basename
    String output_suffix = if make_vcf then ".vcf.gz" else ".g.vcf.gz"
    String output_filename = vcf_basename + output_suffix

    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
    # If we take the number we are scattering by and reduce by 20 we will have enough disk space
    # to account for the fact that the data is quite uneven across the shards.
    Int potential_hc_divisor = length(scattered_calling_intervals) - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

    Int machine_mem = small_task_mem * 1000
    Int disk = small_task_disk + disk_pad
    Int max_retries_or_default = select_first([max_retries, 2])
    Int cpu = small_task_cpu
    Int command_mem = small_task_mem * 1000 - 500





    call SplitIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fasta_index,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            #runtime_params = standard_runtime
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts,
            gatk_path = gatk_path,
            gatk_docker = gatk_docker,
            machine_mem = machine_mem,
            disk = disk,
            max_retries = max_retries_or_default,
            cpu = cpu,
            command_mem = command_mem,
            boot_disk_size = boot_disk_size
    }


  if ( is_cram ) {
    call CramToBamTask {
      input:
        input_cram = input_bam,
        sample_name = sample_basename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = gitc_docker,
        samtools_path = samtools_path
    }
  }

  # Call variants in parallel over grouped calling intervals
  scatter (interval_file in SplitIntervals.interval_files) {

    # Generate VCF by interval
    call HaplotypeCaller {
      input:
        input_bam = select_first([CramToBamTask.output_bam, input_bam]),
        input_bam_index = select_first([CramToBamTask.output_bai, input_bam_index]),
        interval_list = interval_file,
        output_filename = output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        make_vcf = make_vcf,
        make_bamout = make_bamout,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }

  # Merge per-interval VCFs
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      output_filename = output_filename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
  }
}


# TASK DEFINITIONS


task SplitIntervals {
    input {
      File intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      Int scatter_count
      String? split_intervals_extra_args
      String gatk_path
      Int boot_disk_size
      String gatk_docker
      Int machine_mem
      Int disk
      Int cpu
      Int command_mem
      Int preemptible_attempts
      Int max_retries



      # runtime
      #Runtime runtime_params
    }

    command {
        set -e

        mkdir interval-files
        ~{gatk_path} --java-options "-Xmx~{command_mem}m" SplitIntervals \
            -R ~{ref_fasta} \
            ~{"-L " + intervals} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{split_intervals_extra_args}
        cp interval-files/*.interval_list .
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: boot_disk_size
        memory: machine_mem + " MB"
        disks: "local-disk " + disk + " HDD"
        preemptible: preemptible_attempts
        maxRetries: max_retries
        cpu: cpu
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}




task CramToBamTask {
  input {
    # Command parameters
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
    String samtools_path
  }
    Float output_bam_size = size(input_cram, "GB") / 0.60
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20

  command {
    set -e
    set -o pipefail

    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
    ~{samtools_path} view -b -o ~{sample_name}.bam -
    ~{samtools_path} index -b ~{sample_name}.bam
    mv ~{sample_name}.bam.bai ~{sample_name}.bai
  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
 }


  output {
    File output_bam = "~{sample_name}.bam"
    File output_bai = "~{sample_name}.bai"
  }
}


# HaplotypeCaller per-sample in GVCF mode turned off
task HaplotypeCaller {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    File interval_list
    String output_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_vcf
    Boolean make_bamout
    Int hc_scatter

    String? gcs_project_for_requester_pays

    String gatk_path
    String? java_options

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
  }

  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

  String vcf_basename = if make_vcf then basename(output_filename, ".vcf") else basename(output_filename, ".gvcf")
  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      description: "a bam file",
      localization_optional: true
    }
    input_bam_index: {
      description: "an index file for the bam input",
      localization_optional: true
    }
  }
  command {
    set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_filename} \
      -contamination ~{default="0" contamination} \
      ~{false="-G AS_StandardAnnotation" true="-G StandardAnnotation -G StandardHCAnnotation" make_vcf} \
      ~{true="-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90" false= "" make_vcf} \
      ~{false="-ERC GVCF" true="" make_vcf} \
      ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
    File bamout = "~{vcf_basename}.bamout.bam"
  }
}

# Merge VCFs generated per-interval for the same sample
task MergeVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename

    String gatk_path

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1

  command {
  set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}