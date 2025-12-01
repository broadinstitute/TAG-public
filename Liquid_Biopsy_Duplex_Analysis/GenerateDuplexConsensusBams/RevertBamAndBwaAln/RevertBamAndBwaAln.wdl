import "RevertBamAndBwaAln/subworkflows/CopyUmiFromReadName.wdl" as CopyUmiFromReadName
import "RevertBamAndBwaAln/subworkflows/RevertSam.wdl" as RevertSam
import "RevertBamAndBwaAln/subworkflows/BwaAlignment.wdl" as bwa_aln
import "RevertBamAndBwaAln/subworkflows/MergeBamAlignment.wdl" as MergeBamAlignment
import "RevertBamAndBwaAln/subworkflows/SamToFastq.wdl" as samtofastq

workflow AlignRawReadsBwaAln {
	File input_bam
  File input_bam_index
    Boolean extract_umis
    String sample_name
    String bwa_path
    String? gitc_docker
    String gitc_docker_or_default = select_first([gitc_docker, "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"])
	File ref_fasta
	File ref_fai
	File ref_dict
    File ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    Int compression_level
    String? bwa_tool = "bwa"
    Int? bwa_cpu = 8
    Int? bwa_disk_gb = 500
    Int? bwa_mem_gb = 32

    ## CopyUmiTask extra args
    String? copy_umi_docker_override
    Int? copy_umi_extra_disk
    Float? copy_umi_extra_mem
    Boolean? remove_umi_from_read_name

    ## RevertSam extra args
    String? revertsam_extra_args
    Int? revertsam_extra_disk
    Float? revertsam_mem

    ## SamToFastq extra args
    Int? samtofastq_premptible = 0
    Int? samtofastq_disk
    Float? samtofastq_mem

    ## MergeBamAlignment extra args
    Int? mba_extra_mem
    Int? mba_extra_disk
    String? mba_extra_args

    ## SortBam extra args
    Int? sortbam_extra_disk
    Float? sortbam_extra_mem


	call GetBwaVersion {
    	input: gitc_docker = gitc_docker_or_default,
             bwa_path = bwa_path
    }
    
    if(extract_umis){
    	call CopyUmiFromReadName.CopyUmiTask as CopyUmiTask {
    		input: bam_file = input_bam,
        	   	   bam_index = input_bam_index,
               	   base_name = sample_name,
                   bloodbiopsydocker = copy_umi_docker_override,
                   disk_pad = copy_umi_extra_disk,
                   extra_mem = copy_umi_extra_mem,
                   remove_umi_from_read_name = remove_umi_from_read_name
    	}
    }
    
    call RevertSam.RevertSam as revertsam_task {
    	input: input_bam = select_first([CopyUmiTask.umi_extracted_bam, input_bam]),
        	   base_name = sample_name,
               ref_fasta = ref_fasta,
               ref_fasta_index = ref_fai,
               ref_fasta_dict = ref_dict,
               additional_args = revertsam_extra_args,
               disk_buffer = revertsam_extra_disk,
               mem = revertsam_mem
    }
    
    call samtofastq.samtofastq as samtofastq_task {
    	input: input_bam = revertsam_task.output_bam,
             num_preempt = samtofastq_premptible,
             disk_space = samtofastq_disk,
             memory = samtofastq_mem
    }
    
    scatter(i in range(length(samtofastq_task.firstEndFastqs))){
      call bwa_aln.BwaAlignment as bwa_alignment {
         input: refFasta = ref_fasta,
                refFastaIndex = ref_fai,
                refFastaDict = ref_dict,
                ref_alt = ref_alt,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_bwt = ref_bwt,
                ref_pac = ref_pac,
                ref_sa = ref_sa,
                firstEndFastq = samtofastq_task.firstEndFastqs[i],
                secondEndFastq = samtofastq_task.secondEndFastqs[i],
                sampleName = sample_name,
                gitc_docker = gitc_docker_or_default,
                cpu = bwa_cpu,
                diskSpaceGb = bwa_disk_gb,
                memoryGb = bwa_mem_gb
      }
    }

    call MergeBamAlignment.MergeBamAlignmentTask as MBATask {
      input: mapped_bam = bwa_alignment.raw_aligned_bam,
             unmapped_bam = revertsam_task.output_bam,
             bwa_version = GetBwaVersion.version,
             bwa_tool = bwa_tool,
             compression_level = compression_level,
             bwa_commandline = bwa_alignment.bwa_command,
             ref_fasta = ref_fasta,
             ref_fasta_index = ref_fai,
             ref_dict = ref_dict,
             output_bam_basename = sample_name,
             extra_mem = mba_extra_mem,
             diskgb_buffer = mba_extra_disk,
             mba_extra_args = mba_extra_args

    }
    
    call sortbam {
    	input: input_bam = MBATask.output_bam,
        output_bam_basename = sample_name,
        diskgb_buffer = sortbam_extra_disk,
        extra_mem = sortbam_extra_mem
    }

  output {
    File bwa_aln_output_bam = sortbam.output_bam
    File bwa_aln_output_bam_index = sortbam.output_bam_index
  }
}

task GetBwaVersion {
   String gitc_docker
   String bwa_path
   Int? preemptible_attempts

   command {
      ${bwa_path} 2>&1 | \
      grep -e '^Version' | \
      sed 's/Version: //'
   }
   runtime {
      docker: gitc_docker
      memory: "1 GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 2])
   }
   output {
      String version = read_string(stdout())
   }
}

task sortbam {
  File input_bam
  String output_bam_basename
  Int? preemptible_tries = 1
  Int? compression_level = 2
  Int? diskgb_buffer
  Int diskSpaceGb = 50 + select_first([diskgb_buffer, 0])
  Float? extra_mem
  Float memory = 10 + select_first([extra_mem, 0])

  command <<<


        set -euxo pipefail


        java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
        SortSam \
        INPUT=${input_bam} \
        OUTPUT=${output_bam_basename}.bam \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        CREATE_MD5_FILE=true \
        MAX_RECORDS_IN_RAM=300000

  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk ${diskSpaceGb} HDD"
    bootDiskSizeGb: 12
    memory: memory + " GB"
    preemptible: select_first([preemptible_tries])
  }

  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }

}