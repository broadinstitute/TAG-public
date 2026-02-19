version 1.0

import "./FormatUtil.wdl" as FormatUtil
import "./Funcotator.wdl" as Funcotator

struct Runtime {
   Int max_retries
   Int preemptible
   Int mem
   Int cpu
   Int boot_disk_size
   Int initial_disk_size
   String docker
}

workflow CollectMNP {
  input {
    String output_prefix
    File input_maf

    # VCF input is optional, but recommended to supply to the wdl.
    # SNV merging tool uses MAF for determining MNPs, but utilizes
    # the corresponding VCF for preserving the flags and tags in a
    # newly generated MNP VCF
    File? input_vcf
    File? input_vcf_index
      
    File tumor_bam
    File tumor_bam_index
    File? normal_bam
    File? normal_bam_index
      
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # SNV merging criteria
    # Minimum overlap ratio among reads haboring neighboring SNPs
    Float min_overlap_ratio = 0.7
    # Maximum alt count difference among neighboring SNPs
    Int max_alt_count_diff = 10
    # Maximum MNP length
    Int max_mnp_len = 6

    # Funcotator parameters
    String reference_version
    File data_sources_tar_gz
    String transcript_selection_mode = "BEST_EFFECT"
    File? transcript_selection_list
    String? gatk_docker
    File? gatk_jar_override    
    File? interval_list
    String? sequencing_center
    String? sequence_source
    String? funcotator_extra_args

    String merge_snv_docker = "us.gcr.io/tag-team-160914/snv_merger:v1"
    String merge_vcf_docker = "us.gcr.io/tag-public/neovax-tag-gatk:v1"

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 50
    Int boot_disk_size = 15
    Int mem = 8
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-snvmerge:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }

  Int ref_size = ceil(size(ref_fasta,"GB") + size(ref_fasta_index,"GB") + size(ref_dict,"GB"))
  Int tumor_size = ceil(size(tumor_bam,"GB") + size(tumor_bam_index,"GB"))
  Int normal_size = ceil(size(normal_bam,"GB") + size(normal_bam_index,"GB"))

  call MergeSNV {
    input:
      output_prefix = output_prefix,
      input_maf = input_maf,
      input_vcf = input_vcf,
      tumor_bam = tumor_bam,
      tumor_bam_index = tumor_bam_index,
      normal_bam = normal_bam,
      normal_bam_index = normal_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      min_overlap_ratio = min_overlap_ratio,
      max_alt_count_diff = max_alt_count_diff,
      max_mnp_len = max_mnp_len,
      disk_size = ref_size + tumor_size + normal_size,
      docker = merge_snv_docker,
      runtime_params = standard_runtime
  }
  
  if(MergeSNV.mnp_counts > 0){
    call Funcotator.Funcotator as FuncotateMNP {
      input:
        output_basename = output_prefix,
        case_name = MergeSNV.tumor_name,
        control_name = MergeSNV.normal_name,
        input_vcf = MergeSNV.output_mnp_vcf,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        data_sources_tar_gz = data_sources_tar_gz,
        reference_version = reference_version,        
        gatk_docker = gatk_docker,
        gatk_jar_override = gatk_jar_override,
        interval_list = interval_list,
        transcript_selection_list = transcript_selection_list,
        transcript_selection_mode = transcript_selection_mode,
        sequencing_center = sequencing_center,
        sequence_source = sequence_source,
        funcotator_extra_args = funcotator_extra_args,
        output_format = "MAF",
        filter_funcotations = true,
        compress = false,
        preemptible = standard_runtime.preemptible,
        max_retries = standard_runtime.max_retries,
        additional_disk = standard_runtime.initial_disk_size,
        boot_disk_size = standard_runtime.boot_disk_size,
        mem = standard_runtime.mem,
        cpu = standard_runtime.cpu
    }
    
    call FormatUtil.AggregateMAF as AggregateMNP {
      input:
        output_prefix = output_prefix,
        input_mafs = [MergeSNV.output_snv_maf, FuncotateMNP.funcotated_output],
        strip_comments = true,
        comment_char = "#",
        sort_by_chromosome = true,
        max_retries = standard_runtime.max_retries,
        preemptible = standard_runtime.preemptible,
        mem = standard_runtime.mem,
        cpu = standard_runtime.cpu,
        boot_disk_size = standard_runtime.boot_disk_size,
        disk_pad = standard_runtime.initial_disk_size
    }
    
    if(defined(input_vcf) && defined(input_vcf_index)){
      call MergeMNPVCF {
        input: 
          input_vcf = select_first([input_vcf]),
          mnp_vcf = MergeSNV.output_mnp_vcf,
          input_vcf_index = select_first([input_vcf_index]),
          output_prefix = output_prefix,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size = ref_size + 50,
          merge_vcf_docker = merge_vcf_docker,
          runtime_params = standard_runtime
      }
    }
  }
  
  

  output {
    # if MNPs are not found, return the original MAF
    File output_maf = select_first([AggregateMNP.output_maf, input_maf])
    
    # SNV merging tool outputs
    Int mnp_counts = MergeSNV.mnp_counts
    File output_mnp_maflite = MergeSNV.output_mnp_maflite
    File output_mnp_vcf = MergeSNV.output_mnp_vcf
    File output_mnp_tool_log = MergeSNV.output_log
    File? mnp_merged_vcf = MergeMNPVCF.output_merged_vcf
    File? mnp_merged_vcf_index = MergeMNPVCF.output_merged_vcf_index
  }
}

task MergeSNV {
  input {
    String output_prefix  
    File input_maf
    File? input_vcf
    
    File tumor_bam
    File tumor_bam_index
    File? normal_bam
    File? normal_bam_index
    File ref_fasta
    File ref_fasta_index
    Int disk_size
    Runtime runtime_params
    String docker = "us.gcr.io/tag-team-160914/snv_merger:v1"
    # MNP parameters
    Float min_overlap_ratio
    Int max_alt_count_diff
    Int max_mnp_len
  }
   
  Int disk_gb = runtime_params.initial_disk_size + disk_size + ceil(size(input_maf,"GB") + size(input_vcf,"GB"))

  command <<<
    set -e
    # get start position column name - this should be able to handle both 'Start_position' and 'Start_Position'
    START_POS_NAME=`grep ^Hugo ~{input_maf} | tr "\t" "\n" | grep ^Start_[Pp]osition`
      
    python /scripts/snv_merger.py \
      --input-maf ~{input_maf} ~{"--input-vcf " + input_vcf} \
      --tumor-bam ~{tumor_bam} ~{"--normal-bam " + normal_bam} \
      --ref-fasta ~{ref_fasta} \
      --output-prefix ~{output_prefix} \
      --start-pos-name $START_POS_NAME \
      --min-overlap-ratio ~{min_overlap_ratio} \
      --max-altc-diff ~{max_alt_count_diff} \
      --max-mnp-len ~{max_mnp_len}

      # check whether MNPs are found by reading a log file
      grep "FOUND_MNP=" snv_merger.log | cut -f2 -d"=" > mnp_counts.txt
      grep "TUMOR_SAMPLE_NAME=" snv_merger.log | cut -f2 -d"=" > tumor_sample_name.txt
      grep "NORMAL_SAMPLE_NAME=" snv_merger.log | cut -f2 -d"=" | grep -v None > normal_sample_name.txt
  >>>
   
  output {
      File output_snv_maf = "~{output_prefix}.snv.maf"
      File output_mnp_maflite = "~{output_prefix}.mnp.maflite"
      File output_mnp_vcf = "~{output_prefix}.mnp.vcf"
      Int mnp_counts = read_int("mnp_counts.txt")
      String tumor_name = read_string("tumor_sample_name.txt")
      String normal_name = read_string("normal_sample_name.txt")
      File output_log = "snv_merger.log"
  }
  
  runtime {
    docker        : docker
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

task MergeMNPVCF {
  input {
    String output_prefix  
    File mnp_vcf
    File input_vcf
    File input_vcf_index
    
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int disk_size
    Runtime runtime_params
    
    String merge_vcf_docker = "us.gcr.io/tag-public/neovax-tag-gatk:v1"
  
  }
   
  Int disk_gb = runtime_params.initial_disk_size + disk_size + ceil(size(input_vcf,"GB"))
  Int compute_mem = runtime_params.mem * 1000 - 500

  command <<<
    set -e
    
    gatk --java-options "-Xmx~{compute_mem}m" SelectVariants \
    -R ~{ref_fasta} -V ~{input_vcf} \
    -XL ~{mnp_vcf} -O ~{output_prefix}.snp.vcf.gz
    
    gatk --java-options "-Xmx~{compute_mem}m" MergeVcfs \
    -I ~{mnp_vcf} -I ~{output_prefix}.snp.vcf.gz \
    -O ~{output_prefix}.merged.vcf.gz \
    -R ~{ref_fasta}
  >>>
   
  output {
      File output_merged_vcf = "~{output_prefix}.merged.vcf.gz"
      File output_merged_vcf_index = "~{output_prefix}.merged.vcf.gz.tbi"
  }
  
  runtime {
    docker        : merge_vcf_docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}