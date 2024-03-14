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

workflow Funcotator {
  input {
    String output_basename
    String? case_name
    String? control_name

    File input_vcf
    File? input_vcf_index

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    # if gnomAD needs to be added, include it into the database beforehand
    File data_sources_tar_gz
    String reference_version

    String output_format = "MAF"
    Boolean compress = false
    String transcript_selection_mode = "BEST_EFFECT"

    # only report passing variants
    Boolean filter_funcotations = true

    File? interval_list
    File? transcript_selection_list
    String? sequencing_center
    String? sequence_source
    String? funcotator_extra_args

    # workflow level runtime parameters
    String? gatk_docker
    File? gatk_jar_override
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 50
    Int boot_disk_size = 15
    Int mem = 8
    Int cpu = 1
  }

  String docker = select_first([gatk_docker, "us.gcr.io/tag-public/neovax-tag-gatk:v1"])
  Runtime standard_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": docker,
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }

  call Funcotate {
    input:
      output_basename = output_basename,
      case_name = case_name,
      control_name = control_name,
      sequencing_center = sequencing_center,
      sequence_source = sequence_source,
      filter_funcotations = filter_funcotations,
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      reference_version = reference_version,
      output_format = output_format,
      compress = compress,
      interval_list = interval_list,
      data_sources_tar_gz = data_sources_tar_gz,
      transcript_selection_mode = transcript_selection_mode,
      transcript_selection_list = transcript_selection_list,
      funcotator_extra_args = funcotator_extra_args,
      gatk_jar_override = gatk_jar_override,
      runtime_params = standard_runtime
  }

  output {
    File funcotated_output = Funcotate.funcotated_output
    File funcotated_output_index = Funcotate.funcotated_output_index
  }
}


task Funcotate {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String output_basename
    String? case_name
    String? control_name
    File input_vcf
    File? input_vcf_index

    String reference_version
    File data_sources_tar_gz

    String output_format
    Boolean compress

    String? sequencing_center
    String? sequence_source
    String transcript_selection_mode
    File? transcript_selection_list
    Boolean filter_funcotations
    File? interval_list
    String? funcotator_extra_args

    File? gatk_jar_override
    Runtime runtime_params
  }

  String output_maf = output_basename + ".maf"
  String output_maf_index = output_maf + ".idx"
  String output_vcf = output_basename + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"

  String output_file = if output_format == "MAF" then output_maf else output_vcf
  String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx

  Int ref_size = ceil(size(ref_fasta,"GB") + size(ref_fasta_index,"GB") + size(ref_dict,"GB"))
  Int vcf_size = ceil(size(input_vcf,"GB") + size(input_vcf_index,"GB")) * 4
  
  Int compute_mem = runtime_params.mem * 1000 - 500
  Int disk_gb = runtime_params.initial_disk_size + ref_size + vcf_size + ceil(size(data_sources_tar_gz, "GB"))

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar_override}

    # check output format and exit if it's not MAF or VCF
    if [[ "~{output_format}" != "MAF" ]] && [[ "~{output_format}" != "VCF" ]] ; then
      echo "Error: output format must be MAF or VCF."
      exit 1
    fi

    # validate input VCF and modify it for exposing tumor and normal sample names to Funcotator
    INPUT_VCF_PATH="~{input_vcf}"
    if [[ -n "~{case_name}" ]] && [[ -n "~{control_name}" ]] ; then
      echo "Validating input VCF..."
      INPUT_VCF_PATH="~{output_basename}.validated.vcf"
      python /gatk/scripts/validate_funcovcf.py \
        --tumor "~{case_name}" \
        --normal "~{control_name}" \
        --output $INPUT_VCF_PATH \
        "~{input_vcf}"
      gatk --java-options "-Xmx~{compute_mem}m" IndexFeatureFile -I $INPUT_VCF_PATH
    fi

    # extract funcotator database
    echo "Extracting data sources tar/gzip file..."
    mkdir datasources_dir
    tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
    DATA_SOURCES_FOLDER="$PWD/datasources_dir"

    gatk --java-options "-Xmx~{compute_mem}m" Funcotator \
      --data-sources-path $DATA_SOURCES_FOLDER \
      --ref-version ~{reference_version} \
      --output-file-format ~{output_format} \
      -R ~{ref_fasta} \
      -V $INPUT_VCF_PATH \
      -O ~{output_file} \
      ~{"-L " + interval_list} \
      --annotation-default tumor_barcode:~{default="Unknown" case_name} \
      --annotation-default normal_barcode:~{default="Unknown" control_name} \
      --annotation-default Center:~{default="Unknown" sequencing_center} \
      --annotation-default source:~{default="Unknown" sequence_source} \
      --transcript-selection-mode ~{transcript_selection_mode} \
      ~{"--transcript-list " + transcript_selection_list} \
      ~{true="--remove-filtered-variants" false="" filter_funcotations} \
      ~{funcotator_extra_args}

      # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
      if [[ "~{output_format}" == "MAF" ]] ; then
        touch ~{output_maf_index}
      fi
  >>>

  output {
    File funcotated_output = "~{output_file}"
    File funcotated_output_index = "~{output_file_index}"
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
    email: "jtsuji@broadinstitute.org"
  }
}
