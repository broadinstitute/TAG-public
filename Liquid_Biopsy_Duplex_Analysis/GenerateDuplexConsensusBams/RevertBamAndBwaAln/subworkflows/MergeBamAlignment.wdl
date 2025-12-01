workflow RunMBA{
	File sample_name
    
	call MergeBamAlignmentTask{
    	input: output_bam_basename = sample_name
    }
    
    call sortbam {
    	input: input_bam = MergeBamAlignmentTask.output_bam,
        output_bam_basename = sample_name
    }
}

task MergeBamAlignmentTask {
  Array[File] mapped_bam
  File unmapped_bam
  Array[String] bwa_commandline
  String bwa_version
  String bwa_tool
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Int? extra_mem
  String? mba_extra_args
  Int memGb = 64 + select_first([extra_mem,0])
  String? sort_order = "coordinate"

  Int? diskgb_buffer
  Float disk_size = 50 + size(mapped_bam[0], "GB")*length(mapped_bam)*3 + size(unmapped_bam, "GB")*3 + select_first([diskgb_buffer, 0])
  Int compression_level
  Int? preemptible_tries = 1
  String? gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
  Int? cpu = 16

  command <<<
    set -o pipefail
    set -e

    /gatk/gatk \
        MergeBamAlignment \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM ${sep=" --ALIGNED_BAM " mapped_bam} \
        --UNMAPPED_BAM ${unmapped_bam} \
        --OUTPUT ${output_bam_basename}.bam \
        --REFERENCE_SEQUENCE ${ref_fasta} \
        --PAIRED_RUN true \
        --SORT_ORDER ${sort_order} \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --CLIP_ADAPTERS false \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --PROGRAM_RECORD_ID "${bwa_tool}" \
        --PROGRAM_GROUP_VERSION "${bwa_version}" \
        --PROGRAM_GROUP_COMMAND_LINE "${sep=' / ' bwa_commandline}" \
        --PROGRAM_GROUP_NAME "${bwa_tool}" \
        --ADD_PG_TAG_TO_READS false \
        ${mba_extra_args}
     
     du --block-size=kB ${output_bam_basename}.bam | \
     awk -F "kB" '{print $1/1000000}' > output_bam_size.txt
  >>>
  runtime {
    preemptible: select_first([preemptible_tries])
    memory: memGb + " GB"
    bootDiskSizeGb: 12
    docker: select_first([gatk_docker])
    cpu: select_first([cpu])
    disks: "local-disk " + ceil(disk_size) + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    Float output_bam_size = read_float("output_bam_size.txt")
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