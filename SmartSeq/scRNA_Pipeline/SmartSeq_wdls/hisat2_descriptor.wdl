## QC pipeline
## HISAT2 as aligner to align reads to genome reference
## output: genome reference aligned bam file 
## Picard will produce a set of QC metricss
## output: a set of metrics files. 
workflow RunHisat2Pipeline {
  File fastq_read1
  File fastq_read2
  File gtf
  File stranded
  File ref_fasta
  File rrna_interval
  File ref_flat
  File hisat2_ref
  String output_prefix
  String hisat2_ref_name
  String sample_name
  ## variables to estimate disk size
  Float hisat2_ref_size = size(hisat2_ref, "GB")
  Float fastq_size = size(fastq_read1, "GB") + size(fastq_read2, "GB")
  Float reference_bundle_size = size(ref_fasta, "GB") + size(ref_flat, "GB") + size(rrna_interval, "GB") + size(gtf, "GB")
  Float bam_disk_multiplier = 10.0
  Int? increase_disk_size
  Float? increase_mem
  Int? input_threads
  Int threads = select_first([input_threads,4])
  Int additional_disk = select_first([increase_disk_size, 10])
  Float memory = 5 + select_first([increase_mem,0])
  
  Boolean zip_fq1 = if basename(fastq_read1,".gz") == basename(fastq_read1) then true else false
  Boolean zip_fq2 = if basename(fastq_read2,".gz") == basename(fastq_read2) then true else false
 
  call HISAT2PE as Hisat2 {
    input:
      hisat2_ref = hisat2_ref,
      fq1 = fastq_read1,
      fq2 = fastq_read2,
      zip_fq1 = zip_fq1,
      zip_fq2 = zip_fq2,
      ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_name = output_prefix,
      disk_size = ceil(hisat2_ref_size + fastq_size * bam_disk_multiplier + additional_disk * 5.0),
      memory = memory,
      threads = threads
  }
  
  Float bam_size = size(Hisat2.output_bam, "GB")

  call CollectMultipleMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename = output_prefix,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)
	}
  call CollectRnaMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_flat = ref_flat,
      rrna_interval = rrna_interval,
      output_filename = output_prefix,
      stranded = stranded,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)
  }
  call CollectDuplicationMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      output_filename = output_prefix,
      disk_size = ceil((bam_disk_multiplier * bam_size) + additional_disk)
  }
  
  output {
    File aligned_bam = Hisat2.output_bam
    File metfile = Hisat2.metfile
    File logfile = Hisat2.logfile
    File? alignment_summary_metrics = CollectMultipleMetrics.alignment_summary_metrics
    File? base_call_dist_metrics = CollectMultipleMetrics.base_call_dist_metrics
    File? base_call_pdf = CollectMultipleMetrics.base_call_pdf
    File? gc_bias_detail_metrics = CollectMultipleMetrics.gc_bias_detail_metrics
    File? gc_bias_dist_pdf = CollectMultipleMetrics.gc_bias_dist_pdf
    File? gc_bias_summary_metrics = CollectMultipleMetrics.gc_bias_summary_metrics
    File? insert_size_hist = CollectMultipleMetrics.insert_size_hist
    File? insert_size_metrics = CollectMultipleMetrics.insert_size_metrics
    File? quality_distribution_metrics = CollectMultipleMetrics.quality_distribution_metrics
    File? quality_distribution_dist_pdf = CollectMultipleMetrics.quality_distribution_dist_pdf
    File? quality_by_cycle_metrics = CollectMultipleMetrics.quality_by_cycle_metrics
    File? quality_by_cycle_pdf = CollectMultipleMetrics.quality_by_cycle_pdf
    File? pre_adapter_details_metrics = CollectMultipleMetrics.pre_adapter_details_metrics
    File? bait_bias_detail_metrics = CollectMultipleMetrics.bait_bias_detail_metrics
    File? bait_bias_summary_metrics = CollectMultipleMetrics.bait_bias_summary_metrics
    File? error_summary_metrics = CollectMultipleMetrics.error_summary_metrics
    File? rna_metrics = CollectRnaMetrics.rna_metrics
    File? rna_coverage = CollectRnaMetrics.rna_coverage_pdf
    File? dedup_metrics = CollectDuplicationMetrics.dedup_metrics
  }
}

## paired-end alignment
## run HISAT2 to genome reference with dedault parameters
## --seed to fix pseudo-random number and in order to produce deterministics results
## -k --secondary to output multiple mapping reads. --keep 10 will output up to 10 multiple mapping reads, which is default in HISAT2
task HISAT2PE {
  File hisat2_ref
  File fq1  # gz forward fastq file
  File fq2  # gz reverse fastq file
  Boolean zip_fq1
  Boolean zip_fq2
  String ref_name
  String output_name
  String sample_name
  Int disk_size
  Float memory
  Int threads
  command {
    # Note that files MUST be gzipped or the module will not function properly
    # This will be addressed in the future either by a change in how Hisat2 functions or a more
    # robust test for compression type.

    set -e

    # zip fastq files if necessary.
    if [ "${zip_fq1}" == "true" ]; then
        FQ1=${fq1}.fastq.gz
        gzip ${fq1} -c > ${fq1}.fastq.gz
    else
        FQ1=${fq1}
    fi
    echo $FQ1

    if [ "${zip_fq2}" == "true" ]; then
        FQ2=${fq2}.fastq.gz
        gzip ${fq2} -c > ${fq2}.fastq.gz
    else
        FQ2=${fq2}
    fi
    echo $FQ2

   tar -zxf "${hisat2_ref}"
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 $FQ1 \
      -2 $FQ2 \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_name}.log \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      --seed 12345 \
      -k 10 \
      --secondary \
      -p ${threads} -S ${output_name}.sam 
    samtools sort -@ ${threads} -O bam -o "${output_name}.bam" "${output_name}.sam" 
		
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: memory + " GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: threads
    preemptible: 5
  }
  output {
	#Boolean reads_aligned = read_boolean(stdout())
    File logfile = "${output_name}.log"
    File metfile = "${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}


task CollectMultipleMetrics {
  File aligned_bam
  File ref_genome_fasta
  String output_filename
  Int disk_size 
  command {
    java -Xmx6g -jar /usr/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_filename}" \
      FILE_EXTENSION=".txt" \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectGcBiasMetrics \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=QualityScoreDistribution \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=CollectQualityYieldMetrics \
      REFERENCE_SEQUENCE="${ref_genome_fasta}" \
      ASSUME_SORTED=true
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    memory:"7.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }
  output {
    File alignment_summary_metrics = "${output_filename}.alignment_summary_metrics.txt"
    File base_call_dist_metrics = "${output_filename}.base_distribution_by_cycle_metrics.txt"
    File base_call_pdf = "${output_filename}.base_distribution_by_cycle.pdf"
    File gc_bias_detail_metrics = "${output_filename}.gc_bias.detail_metrics.txt"
    File gc_bias_dist_pdf = "${output_filename}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_filename}.gc_bias.summary_metrics.txt"
    File insert_size_hist = "${output_filename}.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_filename}.insert_size_metrics.txt"
    File quality_distribution_metrics = "${output_filename}.quality_distribution_metrics.txt"
    File quality_distribution_dist_pdf = "${output_filename}.quality_distribution.pdf"
    File quality_by_cycle_metrics = "${output_filename}.quality_by_cycle_metrics.txt"
    File quality_by_cycle_pdf = "${output_filename}.quality_by_cycle.pdf"
    File pre_adapter_details_metrics = "${output_filename}.pre_adapter_detail_metrics.txt"
    File pre_adapter_summary_metrics = "${output_filename}.pre_adapter_summary_metrics.txt"
    File bait_bias_detail_metrics = "${output_filename}.bait_bias_detail_metrics.txt"
    File bait_bias_summary_metrics = "${output_filename}.bait_bias_summary_metrics.txt"
    File error_summary_metrics = "${output_filename}.error_summary_metrics.txt"
  }
}

task CollectRnaMetrics {
  File aligned_bam
  File ref_flat
  File rrna_interval
  String output_filename
  String stranded
  Int disk_size
  command{
    set -e
    java -Xmx3g -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_filename}.rna_metrics.txt" \
      REF_FLAT="${ref_flat}" \
      RIBOSOMAL_INTERVALS="${rrna_interval}" \
      STRAND_SPECIFICITY=${stranded} \
      CHART_OUTPUT="${output_filename}.rna.coverage.pdf"
    touch "${output_filename}.rna.coverage.pdf"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    memory:"3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }
  output {
    File rna_metrics = "${output_filename}.rna_metrics.txt"
    File rna_coverage_pdf = "${output_filename}.rna.coverage.pdf"
  }
}

## Here are use  -XX:ParallelGCThreads=2 to run MarkDuplication on mutlple
## thread. 
task CollectDuplicationMetrics {
  File aligned_bam
  String output_filename
  Int disk_size
  command {
    java -Xmx6g -XX:ParallelGCThreads=2  -jar /usr/picard/picard.jar  MarkDuplicates \
       VALIDATION_STRINGENCY=SILENT  \
       INPUT=${aligned_bam} \
       OUTPUT="${output_filename}.MarkDuplicated.bam" \
       ASSUME_SORTED=true \
       METRICS_FILE="${output_filename}.duplicate_metrics.txt" \
       REMOVE_DUPLICATES=false
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    memory: "7.5 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }
  output {
    File dedup_metrics = "${output_filename}.duplicate_metrics.txt"
  }
}
