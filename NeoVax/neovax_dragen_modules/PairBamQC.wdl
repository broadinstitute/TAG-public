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

workflow PairBamQC {
  input {
    String pair_name

    String case_name
    File tumor_bam
    File tumor_bam_index

    String? control_name
    File? normal_bam
    File? normal_bam_index

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File target_intervals
    File bait_intervals

    File haplotype_database
    File dbSNP_vcf
    File dbSNP_vcf_index

    String validation_stringency = "LENIENT"
    String? validation_errors_to_ignore

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

  String docker_image = select_first([gatk_docker, "us.gcr.io/tag-public/neovax-tag-gatk:v1"])
  Runtime standard_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": docker_image,
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }

  call PicardWesMetrics as TumorPicardMetrics {
    input:
      bam = tumor_bam,
      bam_index = tumor_bam_index,
      sample_name = case_name,
      ref_fasta = ref_fasta,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      target_intervals = target_intervals,
      bait_intervals = bait_intervals,
      validation_stringency = validation_stringency,
      validation_errors_to_ignore = validation_errors_to_ignore,
      gatk_jar_override = gatk_jar_override,
      runtime_params = standard_runtime
  }

  if(defined(normal_bam)) {
    
    call PicardWesMetrics as NormalPicardMetrics {
      input:
        bam = select_first([normal_bam]),
        bam_index = select_first([normal_bam_index]),
        sample_name = select_first([control_name]),
        ref_fasta = ref_fasta,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        target_intervals = target_intervals,
        bait_intervals = bait_intervals,
        validation_stringency = validation_stringency,
        validation_errors_to_ignore = validation_errors_to_ignore,
        gatk_jar_override = gatk_jar_override,
        runtime_params = standard_runtime
    }

    call CrossCheckLaneFingerprints {
      input:
        pair_name = pair_name,
        tumor_bam = tumor_bam,
        tumor_bam_index = tumor_bam_index,
        normal_bam = select_first([normal_bam]),
        normal_bam_index = select_first([normal_bam_index]),
        gatk_jar_override = gatk_jar_override,
        haplotype_database = haplotype_database,
        validation_stringency = validation_stringency,
        runtime_params = standard_runtime
    }
  }

  output {
    # Picard WES metrics for tumor
    File tumor_picard_metrics = TumorPicardMetrics.picard_metrics
    File tumor_picard_plots = TumorPicardMetrics.picard_plots
    File tumor_pre_adapter_detail_metrics = TumorPicardMetrics.pre_adapter_detail_metrics

    # Picard WES metrics for normal
    File? normal_picard_metrics = NormalPicardMetrics.picard_metrics
    File? normal_picard_plots = NormalPicardMetrics.picard_plots

    # Cross Check Lane Fingerprints
    File? cross_check_fingprt_metrics = CrossCheckLaneFingerprints.crosscheck_metrics
    File? cross_check_fingprt_report = CrossCheckLaneFingerprints.crosscheck_report
    Float? cross_check_fingprt_min_lod_value = CrossCheckLaneFingerprints.crosscheck_min_lod
    String? cross_check_fingprt_min_lod_lanes = CrossCheckLaneFingerprints.crosscheck_min_lod_lanes
  }
}


task PicardWesMetrics {
  input {
    String sample_name
    File bam
    File bam_index

    File target_intervals
    File bait_intervals

    File ref_fasta
    File dbSNP_vcf
    File dbSNP_vcf_index

    String validation_stringency
    String? validation_errors_to_ignore

    File? gatk_jar_override
    Runtime runtime_params
  }

  Int compute_mem = runtime_params.mem - 1
  Int disk_gb = runtime_params.initial_disk_size +
                ceil(size(target_intervals,"GB") + size(bait_intervals,"GB")) +
                ceil(size(bam,"GB") + size(bam_index,"GB") + size(ref_fasta,"GB")) +
                ceil(size(dbSNP_vcf,"GB") + size(dbSNP_vcf_index,"GB") + size(gatk_jar_override,"GB"))

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar_override}

    gatk --java-options "-Xmx~{compute_mem}g" ValidateSamFile \
      --INPUT ~{bam} \
      --OUTPUT ~{sample_name}.bam_validation \
      --MODE SUMMARY \
      --IGNORE_WARNINGS true ~{"--IGNORE " + validation_errors_to_ignore} \
      --REFERENCE_SEQUENCE ~{ref_fasta} \
      --VALIDATION_STRINGENCY ~{validation_stringency}

    gatk --java-options "-Xmx~{compute_mem}g" CollectMultipleMetrics \
      --INPUT ~{bam} \
      --OUTPUT ~{sample_name}.multiple_metrics \
      --REFERENCE_SEQUENCE ~{ref_fasta} \
      --DB_SNP ~{dbSNP_vcf} \
      --PROGRAM CollectAlignmentSummaryMetrics \
      --PROGRAM CollectInsertSizeMetrics \
      --PROGRAM QualityScoreDistribution \
      --PROGRAM MeanQualityByCycle \
      --PROGRAM CollectBaseDistributionByCycle \
      --PROGRAM CollectSequencingArtifactMetrics \
      --PROGRAM CollectQualityYieldMetrics \
      --PROGRAM CollectGcBiasMetrics
        
    gatk --java-options "-Xmx~{compute_mem}g" ConvertSequencingArtifactToOxoG \
      --INPUT_BASE ~{sample_name}.multiple_metrics \
      --OUTPUT_BASE ~{sample_name}.multiple_metrics.converted \
      --REFERENCE_SEQUENCE ~{ref_fasta} \
      --VALIDATION_STRINGENCY ~{validation_stringency}

    gatk --java-options "-Xmx~{compute_mem}g" CollectHsMetrics \
      --INPUT ~{bam} \
      --BAIT_INTERVALS ~{target_intervals} \
      --TARGET_INTERVALS ~{bait_intervals} \
      --OUTPUT ~{sample_name}.hybrid_selection_metrics \
      --VALIDATION_STRINGENCY ~{validation_stringency}
      
    # create metrics and plot directories
    mkdir ~{sample_name}.picard_files ~{sample_name}.picard_plots
    mv ~{sample_name}.bam_validation ~{sample_name}.*_metrics ~{sample_name}.picard_files
    mv ~{sample_name}.*pdf ~{sample_name}.picard_plots
    cp ~{sample_name}.picard_files/~{sample_name}.multiple_metrics.pre_adapter_detail_metrics ~{sample_name}.pre_adapter_detail_metrics
    # archive picard metrics
    tar -cvzf ~{sample_name}.picard_files.tar.gz ~{sample_name}.picard_files
    tar -cvzf ~{sample_name}.picard_plots.tar.gz ~{sample_name}.picard_plots
  >>>

  output {
    File picard_metrics = "~{sample_name}.picard_files.tar.gz"
    File picard_plots = "~{sample_name}.picard_plots.tar.gz"
    File pre_adapter_detail_metrics = "~{sample_name}.pre_adapter_detail_metrics"
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
  }
}


task CrossCheckLaneFingerprints {
  input {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    File haplotype_database
    String validation_stringency

    File? gatk_jar_override
    Runtime runtime_params
  }

  Int compute_mem = runtime_params.mem - 1
  Int disk_gb = runtime_params.initial_disk_size + 
                ceil(size(tumor_bam,"GB") + size(tumor_bam_index,"GB")) +
                ceil(size(normal_bam,"GB") + size(normal_bam_index,"GB")) +
                ceil(size(haplotype_database,"GB") + size(gatk_jar_override,"GB"))

  command <<<
    set -euxo pipefail
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar_override}

    # remove chromosomes in the haplotype dabase if these do not exist in the BAM header
    samtools view -H ~{tumor_bam} | grep ^@SQ | cut -f2 | cut -f2 -d":" > chromosome_names.txt
    python <<CODE
    chroms = [c.rstrip() for c in open('chromosome_names.txt')]
    fo = open('~{pair_name}.database.txt', 'w')

    for hp in open('~{haplotype_database}'):
        hp = hp.rstrip('\n').split('\t')
        
        if hp[0].startswith('@SQ'):
            if not (hp[1].lstrip('SN:') in chroms):
                continue
        elif not hp[0].startswith('@') and not hp[0].startswith('#'):
            if not (hp[0] in chroms):
                continue

        fo.write('\t'.join(hp) + '\n')
    fo.close()
    CODE

    gatk --java-options "-Xmx~{compute_mem}g" CrosscheckFingerprints \
      -I ~{tumor_bam} \
      -I ~{normal_bam} \
      -H ~{pair_name}.database.txt \
      --TMP_DIR /tmp \
      --QUIET false \
      --EXIT_CODE_WHEN_MISMATCH 0 \
      --OUTPUT ~{pair_name}.crosscheck_metrics \
      --VALIDATION_STRINGENCY ~{validation_stringency}

    # generate the crosscheck report html
    grep -v "#" ~{pair_name}.crosscheck_metrics | sed 1d > crosscheck_metrics.input
    python /gatk/scripts/crosscheck_report.py crosscheck_metrics.input
    mv report.html ~{pair_name}.crosscheck_report.html
  >>>

  output {
    File crosscheck_metrics = "~{pair_name}.crosscheck_metrics"
    File crosscheck_report = "~{pair_name}.crosscheck_report.html"
    Float crosscheck_min_lod = read_float("crosscheck_min_lod_value.txt")
    String crosscheck_min_lod_lanes = read_string("crosscheck_min_lod_lanes.txt")
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
  }
}
