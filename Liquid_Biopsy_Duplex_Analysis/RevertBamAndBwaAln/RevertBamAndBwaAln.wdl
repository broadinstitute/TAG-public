import "./subworkflows/CopyUmiFromReadName.wdl" as CopyUmiFromReadName
import "./subworkflows/RevertSam.wdl" as RevertSam
import "./subworkflows/BwaAlignment.wdl" as bwa_aln
import "./subworkflows/MergeBamAlignment.wdl" as MergeBamAlignment
import "./subworkflows/SamToFastq.wdl" as samtofastq

workflow AlignRawReadsBwaAln {
	File input_bam
    File input_bam_index
    Boolean extract_umis
    String sample_name
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

	call GetBwaVersion {
    	input: gitc_docker = gitc_docker_or_default
    }
    
    if(extract_umis){
    	call CopyUmiFromReadName.CopyUmiTask as CopyUmiTask {
    		input: bam_file = input_bam,
        	   	   bam_index = input_bam_index,
               	   base_name = sample_name
    	}
    }
    
    call RevertSam.RevertSam as revertsam_task {
    	input: input_bam = select_first([CopyUmiTask.umi_extracted_bam, input_bam]),
        	   base_name = sample_name,
               ref_fasta = ref_fasta,
               ref_fasta_index = ref_fai,
               ref_fasta_dict = ref_dict
    }
    
    call samtofastq.samtofastq as samtofastq_task {
    	input: input_bam = revertsam_task.output_bam
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
                gitc_docker = gitc_docker_or_default
      }
    }

    call MergeBamAlignment.MergeBamAlignmentTask as MBATask {
      input: mapped_bam = bwa_alignment.raw_aligned_bam,
             unmapped_bam = revertsam_task.output_bam,
             bwa_commandline = bwa_alignment.bwa_command,
             ref_fasta = ref_fasta,
             ref_fasta_index = ref_fai,
             ref_dict = ref_dict,
             output_bam_basename = sample_name
    }
    
    call sortbam {
    	input: input_bam = MBATask.output_bam,
        output_bam_basename = sample_name
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