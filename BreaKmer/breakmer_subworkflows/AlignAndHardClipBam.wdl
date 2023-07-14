workflow AlignAndHardClipBam {
  File input_bam
  File? output_map
  String base_file_name
  Boolean? fix_mate = false

  Int? additional_disk
  Int compression_level = 2

  Int? preemptible_tries = 2
  Int? agg_preemptible_tries = 2

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  Int disk_pad = select_first([additional_disk, 20])

  # Get reference file sizes
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float bwa_ref_size = ref_size + size(ref_amb, "GB") + size(ref_ann, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_sa, "GB")

  # Disk multipliers
  # md_disk_multiplier: Mark Duplicates outputs a slightly smaller aggregated bam
  # sort_sam_disk_multiplier: SortSam spills to disk a lot more because we are only store 300000 records in RAM to make it faster
  #                           Also it spills to disk in an uncompressed format so we need a larger multiplier.
  #                           This multiplier is used in CramToBam as well.
  Float bwa_disk_multiplier = 2.5
  Float md_disk_multiplier = 2.25
  Float sort_sam_disk_multiplier = 3.25
  Int mba_extra_mem = 0
  
  Boolean? remove_duplicates = false

  Float input_size = size(input_bam, "GB")
  String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y " + ref_fasta

  String? gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
  String? picard_docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.5"
  String? bwa_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"

  call GetBwaVersion {
      input: bwa_docker = select_first([bwa_docker])
  }

  if(!defined(output_map)){
    call GenerateOutputMap {
      input:
        input_bam = input_bam,
        disk_size = ceil(input_size) + disk_pad
    }
  }
  File final_output_map = select_first([output_map, GenerateOutputMap.output_map])

  call RevertSam {
    input:
      input_bam = input_bam,
      output_map = final_output_map,
      gatk_docker = select_first([gatk_docker]),
      disk_size = ceil(input_size)  + disk_pad
  }

  scatter (rg_bam in RevertSam.unmapped_bams) {
    String rg_bam_basename = basename(rg_bam,".coord.sorted.unmapped.bam")
    Float rg_bam_size = size(rg_bam, "GB")
    call SortSam {
      input:
        input_bam = rg_bam,
      	fix_mate = select_first([fix_mate]),
        sorted_bam_name = rg_bam_basename + ".unmapped.bam",
        gatk_docker = select_first([gatk_docker]),
        disk_size = ceil(rg_bam_size * (sort_sam_disk_multiplier * 1.75)) + disk_pad
    }
    call CollectQualityYieldMetrics {
      input:
        input_bam = rg_bam,
        picard_docker = select_first([picard_docker]),
        metrics_filename = rg_bam_basename + ".unmapped.quality_yield_metrics",
        disk_size = ceil(rg_bam_size) + disk_pad,
        preemptible_tries = preemptible_tries
    }
  }

  scatter (unmapped_bam in SortSam.sorted_bam) {
    String output_basename = basename(unmapped_bam, ".unmapped.bam")
    Float unmapped_bam_size = size(unmapped_bam, "GB")

    # Map reads to reference
    call SamToFastq {
      input:
        input_bam = unmapped_bam,
        gatk_docker = select_first([gatk_docker]),
        output_basename = output_basename,
        disk_size = ceil(unmapped_bam_size + bwa_ref_size + (bwa_disk_multiplier * unmapped_bam_size)) + disk_pad,
        preemptible_tries = preemptible_tries
    }

  Float input_fastq_size = size(SamToFastq.output_fastq, "GB")

    call BwaMem {
      input:
        input_fastq = SamToFastq.output_fastq,
        output_bam_basename = output_basename + ".aligned.unsorted.mapped",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        bwa_docker = select_first([bwa_docker]),
        disk_size = ceil(input_fastq_size + bwa_ref_size + (bwa_disk_multiplier * input_fastq_size)) + disk_pad,
        preemptible_tries = preemptible_tries
    }

    call MergeBamAlignment{
      input:
        mapped_bam = BwaMem.output_bam,
        unmapped_bam = unmapped_bam,
        gatk_docker = select_first([gatk_docker]),
        bwa_commandline = bwa_commandline,
        output_bam_basename = output_basename + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        extra_mem = mba_extra_mem,
        bwa_version = GetBwaVersion.version,
        disk_size = ceil(unmapped_bam_size + bwa_ref_size + (bwa_disk_multiplier * unmapped_bam_size)) + disk_pad,
        compression_level = compression_level,
        preemptible_tries = preemptible_tries
    }
  }

  # Sum the read group bam sizes to approximate the aggregated bam size
  call SumFloats {
    input:
      sizes = MergeBamAlignment.output_bam_size,
      preemptible_tries = preemptible_tries
  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs
  # and the merged output.
  call MarkDuplicates {
    input:
      input_bams = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      disk_size = ceil(md_disk_multiplier * SumFloats.total_size) + disk_pad,
      compression_level = compression_level,
      remove_duplicates = remove_duplicates,
        gatk_docker = select_first([gatk_docker]),
      preemptible_tries = agg_preemptible_tries
  }

  Float agg_bam_size = size(MarkDuplicates.output_bam, "GB")

  # Sort aggregated+deduped BAM file and fix tags
  # This task spills to disk so we need space for the input bam, the output bam, and any spillage.
  call SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      picard_docker = select_first([picard_docker]),
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      disk_size = ceil(sort_sam_disk_multiplier * agg_bam_size) + disk_pad,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  # QC the final BAM some more (no such thing as too much QC)
  call CollectAggregationMetrics {
    input:
      input_bam = SortSampleBam.output_bam,
      picard_docker = select_first([picard_docker]),
      input_bam_index = SortSampleBam.output_bam_index,
      output_bam_prefix = base_file_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = ceil(agg_bam_size + ref_size) + disk_pad,
      preemptible_tries = agg_preemptible_tries
  }

  output {
    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.metrics
    File duplicate_metrics = MarkDuplicates.duplicate_metrics

    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
    File agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
    File agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
    File agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics

    # main output files
    File output_bam = SortSampleBam.output_bam
    File output_bam_index = SortSampleBam.output_bam_index
  }
}

task GenerateOutputMap {
  File input_bam
  Int disk_size

  command {
    set -e

    samtools view -H ${input_bam} | grep @RG | cut -f2 | sed s/ID:// > readgroups.txt

    echo -e "READ_GROUP_ID\tOUTPUT" > output_map.tsv

    for rg in `cat readgroups.txt`; do
      echo -e "$rg\t$rg.coord.sorted.unmapped.bam" >> output_map.tsv
    done

    ### edytam: adding echo to check if the task ends up correctly
    echo "generation of the output map finished"
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
  }
  output {
    File output_map = "output_map.tsv"
  }
}

task RevertSam {
  File input_bam
  String input_basename = basename(input_bam, ".bam")
  File output_map
  Int disk_size
  Float? extra_mem
  Float memory = 1.2 + select_first([extra_mem,0])
  Int java_mem = floor((memory - 0.2)*1000)

  String gatk_docker

  command <<<
  	set -o pipefail
    set -e
  
    java -Xmx${java_mem}m -jar /root/gatk.jar \
    RevertSam \
    --INPUT ${input_bam} \
    --OUTPUT_MAP ${output_map} \
    --OUTPUT_BY_READGROUP true \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER coordinate
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GB"
  }
  output {
    Array[File] unmapped_bams = glob("*.bam")
  }
}


task SortSam {
  File input_bam
  String sorted_bam_name
  Boolean fix_mate
  Int disk_size
  String gatk_docker

  command <<<
    
    java -Xmx3000m -jar /root/gatk.jar  \
    SortSam \
    --INPUT ${input_bam} \
    --OUTPUT sorted.bam \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 1000000
    
    if [ "${fix_mate}" == "true" ]; then
    	samtools fixmate sorted.bam ${sorted_bam_name}
    else
    	mv sorted.bam ${sorted_bam_name}
    fi
    
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 12
    memory: "3500 MB"
    preemptible: 3
  }
  output {
    File sorted_bam = "${sorted_bam_name}"
  }
}

task CollectQualityYieldMetrics {
  File input_bam
  String metrics_filename
  Float disk_size
  Int preemptible_tries

  String picard_docker

  command {
    java -Xms2000m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=${input_bam} \
      OQ=true \
      OUTPUT=${metrics_filename}
  }
  runtime {
    disks: "local-disk " + ceil(disk_size) + " HDD"
    memory: "3 GB"
    bootDiskSizeGb: 12
    docker: picard_docker
    preemptible: preemptible_tries
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

task SamToFastq {
  File input_bam
  String output_basename

  Float disk_size
  Int preemptible_tries
  String gatk_docker

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    java -Xms5000m -jar /root/gatk.jar \
        SamToFastq \
        --INPUT ${input_bam} \
        --FASTQ ${output_basename}.fastq \
        --INTERLEAVE true \
        -NON_PF true
  >>>
    runtime {
    preemptible: preemptible_tries
    memory: "14 GB"
    bootDiskSizeGb: 12
    docker: gatk_docker
    cpu: "16"
    disks: "local-disk " + ceil(disk_size) + " HDD"
  }
  output {
    File output_fastq = "${output_basename}.fastq"
  }
}

task BwaMem {
  File input_fastq
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  Float disk_size
  Int preemptible_tries
  String bwa_docker
  String? bwa_mem_command = "bwa mem -M -v3 -t 13 -p -K 100000000 -Y"

  command <<<
    set -o pipefail
    set -e

    /usr/gitc/${bwa_mem_command} ${ref_fasta} ${input_fastq} - 2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) > ${output_bam_basename}.bam
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "14 GB"
    bootDiskSizeGb: 12
    docker: bwa_docker
    cpu: "16"
    disks: "local-disk " + ceil(disk_size) + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task MergeBamAlignment {
  File mapped_bam
  File unmapped_bam
  String bwa_commandline
  String bwa_version
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Int? extra_mem
  String? mba_extra_args
  Int? memGb = 64 + select_first([extra_mem,0])
  Int command_mem = floor((memGb - 0.5) * 1000)
  String? sort_order = "coordinate"

  Float disk_size
  Int compression_level
  Int preemptible_tries
  String gatk_docker

  command <<<
    set -o pipefail
    set -e

    java -Dsamjdk.compression_level=${compression_level} -Xms${command_mem}m -jar /root/gatk.jar \
        MergeBamAlignment \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM ${mapped_bam} \
        --UNMAPPED_BAM ${unmapped_bam} \
        --OUTPUT ${output_bam_basename}.bam \
        --REFERENCE_SEQUENCE ${ref_fasta} \
        --PAIRED_RUN true \
        --SORT_ORDER ${sort_order} \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --HARD_CLIP_OVERLAPPING_READS true \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --PROGRAM_RECORD_ID "bwamem" \
        --PROGRAM_GROUP_VERSION "${bwa_version}" \
        --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
        --PROGRAM_GROUP_NAME "bwamem" \
        --ADD_PG_TAG_TO_READS false \
        ${mba_extra_args}
     
     du --block-size=kB ${output_bam_basename}.bam | \
     awk -F "kB" '{print $1/1000000}' > output_bam_size.txt
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: memGb + " GB"
    bootDiskSizeGb: 12
    docker: gatk_docker
    cpu: "16"
    disks: "local-disk " + ceil(disk_size) + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    Float output_bam_size = read_float("output_bam_size.txt")
  }
}

task GetBwaVersion {
    String bwa_docker
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    docker: bwa_docker
    memory: "1 GB"
  }
  output {
    String version = read_string(stdout())
  }
}

task SumFloats {
  Array[Float] sizes
  Int preemptible_tries

  command <<<
  python -c "print ${sep="+" sizes}"
  >>>
  output {
    Float total_size = read_float(stdout())
  }
  runtime {
    disks: "local-disk " + 10 + " HDD"
    memory: "5000 MB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    preemptible: preemptible_tries
  }
}

task MarkDuplicates {
  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  Float disk_size
  Int compression_level
  Int preemptible_tries
  String gatk_docker

  # The program default for READ_NAME_REGEX is appropriate in nearly every case.
  # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
  # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
  String? read_name_regex
  Boolean? remove_duplicates
  

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar /root/gatk.jar \
      MarkDuplicates \
      --INPUT ${sep=' --INPUT ' input_bams} \
      --OUTPUT ${output_bam_basename}.bam \
      --METRICS_FILE ${metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      ${"--READ_NAME_REGEX " + read_name_regex} \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --ADD_PG_TAG_TO_READS false \
      --REMOVE_DUPLICATES ${default="false" remove_duplicates}
  }
  runtime {
    preemptible: preemptible_tries
    docker: gatk_docker
    memory: "7 GB"
    bootDiskSizeGb: 12
    disks: "local-disk " + ceil(disk_size) + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortSampleBam {
  File input_bam
  String output_bam_basename
  Int preemptible_tries
  Int compression_level
  Float disk_size
  Float? extra_mem
  Float memory = 5 + select_first([extra_mem,0])
  Int java_mem = floor(memory-1)*1000
  String picard_docker

  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms${java_mem}m -jar /usr/picard/picard.jar \
      SortSam \
      INPUT=${input_bam} \
      OUTPUT=${output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    disks: "local-disk " + ceil(disk_size) + " HDD"
    cpu: "1"
    bootDiskSizeGb: 12
    docker: picard_docker
    memory: memory+ " GB"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}


task CollectAggregationMetrics {
  File input_bam
  File input_bam_index
  String output_bam_prefix
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int preemptible_tries
  Float disk_size
  String picard_docker

  command {
    java -Xms5000m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=${input_bam} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectAlignmentSummaryMetrics" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="CollectSequencingArtifactMetrics" \
      PROGRAM="CollectGcBiasMetrics" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="SAMPLE" \
      METRIC_ACCUMULATION_LEVEL="LIBRARY"

    touch ${output_bam_prefix}.insert_size_metrics
    touch ${output_bam_prefix}.insert_size_histogram.pdf
  }
  runtime {
    memory: "7 GB"
    docker: picard_docker
    disks: "local-disk " + ceil(disk_size) + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File alignment_summary_metrics = "${output_bam_prefix}.alignment_summary_metrics"
    File bait_bias_detail_metrics = "${output_bam_prefix}.bait_bias_detail_metrics"
    File bait_bias_summary_metrics = "${output_bam_prefix}.bait_bias_summary_metrics"
    File gc_bias_detail_metrics = "${output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "${output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_bam_prefix}.gc_bias.summary_metrics"
    File insert_size_histogram_pdf = "${output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_bam_prefix}.insert_size_metrics"
    File pre_adapter_detail_metrics = "${output_bam_prefix}.pre_adapter_detail_metrics"
    File pre_adapter_summary_metrics = "${output_bam_prefix}.pre_adapter_summary_metrics"
    File quality_distribution_pdf = "${output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "${output_bam_prefix}.quality_distribution_metrics"
  }
}
