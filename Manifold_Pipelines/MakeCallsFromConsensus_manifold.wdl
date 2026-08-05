version 1.0

## ============================================================================
## Combined single-file WDL
##
## This file merges three previously separate WDLs into one document:
##   - MakeCallsFromConsensus.wdl  (top-level workflow)
##   - mutect2.wdl                 (imported as `m2`)
##   - checkBaitSetName.wdl        (imported as `checkBaitSetName`)
##
## WDL permits only ONE workflow per document, so the `Mutect2` sub-workflow
## (previously `call m2.Mutect2 as M2Duplex`) has been expanded inline inside
## MakeCallsFromConsensus. Its workflow-level declarations are prefixed with
## `m2_` to avoid name collisions, and its original runtime defaults are
## preserved (standard tasks use preemptible=2, the M2 scatter uses 10 --
## these are independent of the workflow-level `preemptible` input, exactly as
## in the original where M2Duplex did not forward `preemptible`).
##
## The imported tasks `Funcotate` and `compareBaitSetName` are inlined verbatim.
## The Mutect2 tasks that were never reached in this pipeline (MergeBamOuts and
## FilterAlignmentArtifacts -- guarded by make_bamout / realignment_index_bundle,
## which MakeCallsFromConsensus never sets) have been omitted.
## ============================================================================

struct Runtime {
    String gatk_docker
    File? gatk_override
    Int max_retries
    Int preemptible
    Int cpu
    Int machine_mem
    Int command_mem
    Int disk
    Int boot_disk_size
}

workflow MakeCallsFromConsensus {

   input {
      # dockers and override jars
      String bloodbiopsydocker
      String reference_version
      Int preemptible

      # Mutect2 specific arguments
      String basename
      File tumor_bam
      File tumor_bam_idx
      String tumor_sample_name
      File? normal_bam
      File? normal_bam_idx
      String? normal_sample_name
      File? raw_tumor_bam
      File? raw_tumor_bam_idx
      File? raw_normal_bam
      File? raw_normal_bam_idx

      File target_intervals
      Boolean fail_on_intervals_mismatch
      File reference
      File reference_idx
      File reference_dict
      File? pon
      File? pon_idx
      File? pon_for_benchmarking
      File? pon_for_benchmarking_idx
      File? variants_for_contamination
      File? variants_for_contamination_idx
      File? gnomad
      File? gnomad_idx
      File? gnomad_for_benchmarking
      File? gnomad_for_benchmarking_idx
      Boolean? run_orientation_bias_mixture_model_filter
      Boolean run_ob_filter = select_first([run_orientation_bias_mixture_model_filter, false])
      Int scatter_count
      String? m2_extra_args
      String? m2_extra_filtering_args
      String m2_extra_filtering_args_or_default = select_first([m2_extra_filtering_args, ""])
      Boolean? compress_vcfs
      String mapping_filter_python_script
      File blastdb_nhr
      File blastdb_nin
      File blastdb_nsq
      String blastn_path

      Boolean? is_benchmark

      File? gatk_override

      # Fingerprinting haplotype database
      File haplotype_map
      Boolean? fingerprint_tumor_normal

      ## Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk
      # This is added to every task as padding, should increase if systematically you need more disk for every call
      Int disk_pad = 10 + select_first([emergency_extra_disk,0])
      #if set to true, use the strand bias filter on the raw annotation
      #if set to false (default), uses the strand bias filter on the duplex annotations

      # Funcotator inputs
      Boolean filter_funcotations
      File? data_sources_tar_gz
      String? funcotator_extra_args

      # --- Inlined-Mutect2 runtime knobs (previously set via M2Duplex.* in the inputs JSON) ---
      # Standard M2 helper tasks (SplitIntervals, MergeVCFs, MergeStats, MergePileupSummaries,
      # CalculateContamination, Filter) run with these; the M2 scatter uses m2_task_mem +
      # m2_preemptible / m2_max_retries.
      Int m2_small_task_mem = 8
      Int m2_small_task_disk = 100
      Int m2_small_task_cpu = 2
      Int m2_boot_disk_size = 12
      Int? m2_preemptible
      Int? m2_max_retries
      Float m2_task_mem = 3.5
      Boolean? m2_make_bamout
      String? m2_getpileupsummaries_extra_args
   }

   Int small_task_mem = 4
   Int small_task_disk = 100
   Int boot_disk_size = 12
   Int small_task_cpu = 2
   Int max_retries_or_default = 2

   Runtime standard_runtime = {"gatk_docker": bloodbiopsydocker, "gatk_override": gatk_override,
      "max_retries": max_retries_or_default, "preemptible": preemptible, "cpu": small_task_cpu,
      "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
      "disk": small_task_disk + disk_pad, "boot_disk_size": boot_disk_size}

   if(select_first([fingerprint_tumor_normal, true])) {
      call CrosscheckFingerprints {
         input:
            bloodbiopsydocker = bloodbiopsydocker,
            preemptible = preemptible,
            disk_pad = disk_pad,
            input_bam = tumor_bam,
            second_input_bam = normal_bam,
            haplotype_map = haplotype_map
      }
      call FingerprintsResult{
         input:
            bloodbiopsydocker = bloodbiopsydocker,
            preemptible = preemptible,
            disk_pad = disk_pad,
            fingerprint_metrics = CrosscheckFingerprints.fingerprint_metrics,

      }
   }
   # (was checkBaitSetName.compareBaitSetName) -- `bait_set` remains a required
   # task input, surfaced as MakeCallsFromConsensus.checkBaitSetName.bait_set
   call compareBaitSetName as checkBaitSetName {
      input:
         target_intervals = target_intervals,
         fail_task = fail_on_intervals_mismatch
   }
   # Collect Sequencing Artifact Metrics after deduplication by start and stop
   # positions (but not including UMIs).
   call CollectSequencingArtifactMetrics as ConsensusArtifactMetricsTumor {
      input:
         bam_file = tumor_bam,
         bam_idx = tumor_bam_idx,
         basename = tumor_sample_name,
         reference = reference,
         reference_idx = reference_idx,
         bloodbiopsydocker = bloodbiopsydocker,
         preemptible = preemptible,
         disk_pad = disk_pad
   }

   if (defined(normal_bam)) {
      call CollectSequencingArtifactMetrics as ConsensusArtifactMetricsNormal {
         input:
            bam_file = select_first([normal_bam]),
            bam_idx =  select_first([normal_bam_idx]),
            basename = normal_sample_name,
            reference = reference,
            reference_idx = reference_idx,
            bloodbiopsydocker = bloodbiopsydocker,
            preemptible = preemptible,
            disk_pad = disk_pad
      }
   }

   Boolean is_benchmark_or_default = select_first([is_benchmark, false])
   File? pon_to_use = if is_benchmark_or_default then pon_for_benchmarking else pon
   File? pon_to_use_idx = if is_benchmark_or_default then pon_for_benchmarking_idx else pon_idx

   File? gnomad_to_use = if is_benchmark_or_default then gnomad_for_benchmarking else gnomad
   File? gnomad_to_use_idx = if is_benchmark_or_default then gnomad_for_benchmarking_idx else gnomad_idx

   # =========================================================================
   # Inlined Mutect2 pipeline  (previously: call m2.Mutect2 as M2Duplex)
   #
   # Only filtered_vcf / filtered_vcf_idx / contamination_table are consumed
   # downstream, so the Funcotate/FilterAlignmentArtifacts/MergeBamOuts branches
   # of the original Mutect2 workflow (all disabled for this pipeline) are not
   # reproduced. Runtime defaults are kept identical to the original sub-workflow.
   # =========================================================================

   Int m2_preemptible_or_default = select_first([m2_preemptible, 2])
   Int m2_max_retries_or_default = select_first([m2_max_retries, 2])
   Boolean m2_compress = select_first([compress_vcfs, false])
   Boolean m2_make_bamout_or_default = select_first([m2_make_bamout, false])
   Int learn_read_orientation_mem = 8000
   Float small_input_to_output_multiplier = 2.0
   Float cram_to_bam_multiplier = 6.0

   # The original Mutect2 sub-workflow computed its own disk padding as
   # 10 + gatk_override_size + emergency_extra_disk. M2Duplex forwarded neither
   # a gatk_override nor emergency_extra_disk, so that value was always 10.
   # Kept as a constant (separate from the workflow-level disk_pad) to reproduce
   # the original Mutect2 disk sizing exactly, even when emergency_extra_disk is set.
   Int m2_disk_pad = 10

   Int m2_ref_size = ceil(size(reference, "GB") + size(reference_dict, "GB") + size(reference_idx, "GB"))
   Int m2_tumor_reads_size = ceil(size(tumor_bam, "GB") + size(tumor_bam_idx, "GB"))
   Int m2_gnomad_vcf_size = if defined(gnomad_to_use) then ceil(size(gnomad_to_use, "GB")) else 0
   Int m2_normal_reads_size = if defined(normal_bam) then ceil(size(normal_bam, "GB") + size(normal_bam_idx, "GB")) else 0

   String m2_output_basename = basename(basename(tumor_bam, ".bam"), ".cram")  #hacky way to strip either .bam or .cram
   String m2_unfiltered_name = m2_output_basename + "-unfiltered"
   String m2_filtered_name = m2_output_basename + "-filtered"

   Int m2_tumor_cram_to_bam_disk = ceil(m2_tumor_reads_size * cram_to_bam_multiplier)
   Int m2_normal_cram_to_bam_disk = ceil(m2_normal_reads_size * cram_to_bam_multiplier)

   # NOTE: gatk_override is intentionally omitted (optional struct member) so the
   # inlined Mutect2 tasks fall back to /root/gatk.jar, matching the original
   # where M2Duplex did not forward a gatk_override.
   Runtime m2_standard_runtime = {"gatk_docker": bloodbiopsydocker,
      "max_retries": m2_max_retries_or_default, "preemptible": m2_preemptible_or_default, "cpu": m2_small_task_cpu,
      "machine_mem": m2_small_task_mem * 1000, "command_mem": m2_small_task_mem * 1000 - 500,
      "disk": m2_small_task_disk + m2_disk_pad, "boot_disk_size": m2_boot_disk_size}

   if (basename(tumor_bam) != basename(tumor_bam, ".cram")) {
      call CramToBam as TumorCramToBam {
         input:
            ref_fasta = reference,
            ref_fai = reference_idx,
            ref_dict = reference_dict,
            cram = tumor_bam,
            crai = tumor_bam_idx,
            name = m2_output_basename,
            disk_size = m2_tumor_cram_to_bam_disk
      }
   }

   String m2_normal_or_empty = select_first([normal_bam, ""])
   if (basename(m2_normal_or_empty) != basename(m2_normal_or_empty, ".cram")) {
      String m2_normal_basename = basename(basename(m2_normal_or_empty, ".bam"), ".cram")
      call CramToBam as NormalCramToBam {
         input:
            ref_fasta = reference,
            ref_fai = reference_idx,
            ref_dict = reference_dict,
            cram = normal_bam,
            crai = normal_bam_idx,
            name = m2_normal_basename,
            disk_size = m2_normal_cram_to_bam_disk
      }
   }

   File m2_tumor_bam = select_first([TumorCramToBam.output_bam, tumor_bam])
   File m2_tumor_bai = select_first([TumorCramToBam.output_bai, tumor_bam_idx])
   File? m2_normal_bam = if defined(normal_bam) then select_first([NormalCramToBam.output_bam, normal_bam]) else normal_bam
   File? m2_normal_bai = if defined(normal_bam) then select_first([NormalCramToBam.output_bai, normal_bam_idx]) else normal_bam_idx

   Int m2_tumor_bam_size = ceil(size(m2_tumor_bam, "GB") + size(m2_tumor_bai, "GB"))
   Int m2_normal_bam_size = if defined(m2_normal_bam) then ceil(size(m2_normal_bam, "GB") + size(m2_normal_bai, "GB")) else 0

   Int m2_output_size = m2_tumor_bam_size / scatter_count
   Int m2_per_scatter_size = (m2_tumor_bam_size + m2_normal_bam_size) + m2_ref_size + m2_gnomad_vcf_size + m2_output_size + m2_disk_pad

   call SplitIntervals {
      input:
         intervals = target_intervals,
         ref_fasta = reference,
         ref_fai = reference_idx,
         ref_dict = reference_dict,
         scatter_count = scatter_count,
         runtime_params = m2_standard_runtime
   }

   scatter (subintervals in SplitIntervals.interval_files) {
      call M2 {
         input:
            intervals = subintervals,
            ref_fasta = reference,
            ref_fai = reference_idx,
            ref_dict = reference_dict,
            tumor_bam = m2_tumor_bam,
            tumor_bai = m2_tumor_bai,
            normal_bam = m2_normal_bam,
            normal_bai = m2_normal_bai,
            pon = pon_to_use,
            pon_idx = pon_to_use_idx,
            gnomad = gnomad_to_use,
            gnomad_idx = gnomad_to_use_idx,
            m2_extra_args = m2_extra_args,
            variants_for_contamination = variants_for_contamination,
            make_bamout = m2_make_bamout_or_default,
            run_ob_filter = run_ob_filter,
            mem = m2_task_mem,
            compress = m2_compress,
            getpileupsummaries_extra_args = m2_getpileupsummaries_extra_args,
            variants_for_contamination_idx = variants_for_contamination_idx,
            preemptible = m2_preemptible,
            max_retries = m2_max_retries,
            gatk_docker = bloodbiopsydocker,
            disk_space = m2_per_scatter_size
      }
   }

   if (run_ob_filter) {
      call LearnReadOrientationModel {
         input:
            f1r2_tar_gz = M2.f1r2_counts,
            runtime_params = m2_standard_runtime,
            mem = learn_read_orientation_mem
      }
   }

   call MergeVCFs {
      input:
         input_vcfs = M2.unfiltered_vcf,
         input_vcf_indices = M2.unfiltered_vcf_idx,
         output_name = m2_unfiltered_name,
         compress = m2_compress,
         runtime_params = m2_standard_runtime
   }

   call MergeStats { input: stats = M2.stats, runtime_params = m2_standard_runtime }

   if (defined(variants_for_contamination)) {
      call MergePileupSummaries as MergeTumorPileups {
         input:
            input_tables = flatten(M2.tumor_pileups),
            output_name = m2_output_basename,
            ref_dict = reference_dict,
            runtime_params = m2_standard_runtime
      }

      if (defined(m2_normal_bam)){
         call MergePileupSummaries as MergeNormalPileups {
            input:
               input_tables = flatten(M2.normal_pileups),
               output_name = m2_output_basename,
               ref_dict = reference_dict,
               runtime_params = m2_standard_runtime
         }
      }

      call CalculateContamination {
         input:
            tumor_pileups = MergeTumorPileups.merged_table,
            normal_pileups = MergeNormalPileups.merged_table,
            runtime_params = m2_standard_runtime
      }
   }

   call Filter {
      input:
         ref_fasta = reference,
         ref_fai = reference_idx,
         ref_dict = reference_dict,
         intervals = target_intervals,
         unfiltered_vcf = MergeVCFs.merged_vcf,
         unfiltered_vcf_idx = MergeVCFs.merged_vcf_idx,
         output_name = m2_filtered_name,
         compress = m2_compress,
         mutect_stats = MergeStats.merged_stats,
         contamination_table = CalculateContamination.contamination_table,
         maf_segments = CalculateContamination.maf_segments,
         artifact_priors_tar_gz = LearnReadOrientationModel.artifact_prior_table,
         m2_extra_filtering_args = m2_extra_filtering_args_or_default,
         runtime_params = m2_standard_runtime,
         disk_space = ceil(size(MergeVCFs.merged_vcf, "GB") * small_input_to_output_multiplier) + m2_disk_pad
   }

   # Resolved outputs of the inlined Mutect2 pipeline (formerly M2Duplex.*)
   File m2_filtered_vcf = Filter.filtered_vcf
   File m2_filtered_vcf_idx = Filter.filtered_vcf_idx
   File? m2_contamination_table = CalculateContamination.contamination_table
   # =========================================================================
   # End inlined Mutect2 pipeline
   # =========================================================================

   if(defined(m2_contamination_table)){
       call ExtractContam {
           input: bloodbiopsydocker = bloodbiopsydocker,
                  preemptible = preemptible,
                  contamination_table = select_first([m2_contamination_table])
       }
   }

   # Apply MRD mapping filter to duplex
   call RunMappingFilter {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         blastn_path = blastn_path,
         mapping_filter_python_script = mapping_filter_python_script,
         vcf_file = m2_filtered_vcf,
         vcf_idx = m2_filtered_vcf_idx,
         reference = reference,
         reference_idx = reference_idx,
         reference_dict = reference_dict,
         blastdb_nhr = blastdb_nhr,
         blastdb_nin = blastdb_nin,
         blastdb_nsq = blastdb_nsq,
         preemptible = preemptible,
         disk_pad = disk_pad
   }

   # Run additional filtering tasks
   call VariantFiltration {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         mapping_filter_name = "mapping_filter",
         duplex_vcf_file = m2_filtered_vcf,
         duplex_vcf_idx = m2_filtered_vcf_idx,
         filter_vcf_file = RunMappingFilter.map_filtered_vcf,
         filter_vcf_idx = RunMappingFilter.map_filtered_vcf_idx,
         reference_fasta = reference,
         reference_fasta_idx = reference_idx,
         reference_dict = reference_dict,
         preemptible = preemptible,
         filter_not_in_mask = true,
         disk_pad = disk_pad
   }

   # Split VCFs by snps and indels
   call SplitVCFs {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         vcf_file = VariantFiltration.output_vcf,
         vcf_idx = VariantFiltration.output_vcf_idx,
         reference_fasta = reference,
         reference_fasta_idx = reference_idx,
         reference_dict = reference_dict,
         preemptible = preemptible,
         disk_pad = disk_pad
   }

   # Create text file of variants
   call VariantsToTable {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         vcf = SplitVCFs.output_snp_vcf,
         vcf_idx = SplitVCFs.output_snp_vcf_idx,
         output_name = basename,
         preemptible = preemptible,
         disk_pad = disk_pad
   }

   # Create annotated maf file containing SNPs  (was m2.Funcotate)
   # NOTE: ref_fai and ref_dict are intentionally left unbound to match the
   # original interface -- they surface as required workflow inputs
   # (MakeCallsFromConsensus.FuncotateMafSnps.ref_fai / .ref_dict).
   call Funcotate as FuncotateMafSnps {
      input:
         ref_fasta = reference,
         input_vcf = SplitVCFs.output_snp_vcf,
         input_vcf_idx = SplitVCFs.output_snp_vcf_idx,
         case_id = tumor_sample_name,
         control_id = normal_sample_name,
         reference_version = reference_version,
         data_sources_tar_gz = data_sources_tar_gz,
         filter_funcotations = filter_funcotations,
         extra_args = funcotator_extra_args,
         use_gnomad = false,
         output_format = "MAF",
         output_file_base_name = basename + "_snp",
         compress = true,
         runtime_params = standard_runtime
   }

   # Create annotated maf file containing indels  (was m2.Funcotate)
   # NOTE: ref_fai and ref_dict are intentionally left unbound to match the
   # original interface -- they surface as required workflow inputs
   # (MakeCallsFromConsensus.FuncotateMafIndels.ref_fai / .ref_dict).
   call Funcotate as FuncotateMafIndels {
      input:
         ref_fasta = reference,
         input_vcf = SplitVCFs.output_indel_vcf,
         input_vcf_idx = SplitVCFs.output_indel_vcf_idx,
         case_id = tumor_sample_name,
         control_id = normal_sample_name,
         reference_version = reference_version,
         data_sources_tar_gz = data_sources_tar_gz,
         filter_funcotations = filter_funcotations,
         extra_args = funcotator_extra_args,
         use_gnomad = false,
         output_format = "MAF",
         output_file_base_name = basename + "_indel",
         compress = true,
         runtime_params = standard_runtime
   }

   Array[String] input_files = select_all([SplitVCFs.output_snp_vcf, SplitVCFs.output_indel_vcf, tumor_bam, normal_bam])

   # Create IGV session so that it is easy to examine results of pipeline
   call GenerateIGVSession {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         preemptible = preemptible,
         input_files = input_files,
         file_name =  basename,
         reference_version = reference_version
   }

   output {

      File igv_session = GenerateIGVSession.igv_session

      File filtered_vcf = VariantFiltration.output_vcf
      File filtered_vcf_idx = VariantFiltration.output_vcf_idx

      File filtered_snp_vcf = SplitVCFs.output_snp_vcf
      File filtered_snp_vcf_idx = SplitVCFs.output_snp_vcf_idx

      File filtered_indel_vcf = SplitVCFs.output_indel_vcf
      File filtered_indel_vcf_idx = SplitVCFs.output_indel_vcf_idx

      Int n_passing_snps = SplitVCFs.passing_SNP
      Int n_filtered_snps = SplitVCFs.filtered_SNP

      Int n_passing_mnps = SplitVCFs.passing_MNP
      Int n_filtered_mnps = SplitVCFs.filtered_MNP

      Int n_passing_indels = SplitVCFs.passing_INDEL
      Int n_filtered_indels = SplitVCFs.filtered_INDEL

      File? contamination_table = m2_contamination_table
      Float? contamination_fraction = ExtractContam.fracContam

      File funcotated_snp_maf = FuncotateMafSnps.funcotated_output_file
      File funcotated_indel_maf = FuncotateMafIndels.funcotated_output_file

      File variant_table = VariantsToTable.output_table
      File? fingerprint_metrics = CrosscheckFingerprints.fingerprint_metrics
      Int? expected_match = FingerprintsResult.expected_match

   }
}

## ============================================================================
## Tasks from checkBaitSetName.wdl
## ============================================================================

task compareBaitSetName {
  input {
    String bait_set
    Boolean fail_task
    File target_intervals
    File? bait_intervals
  }
  String target_intervals_name = basename(target_intervals)
  String bait_intervals_name = if defined(bait_intervals) then basename(select_first([bait_intervals])) else bait_set + "."

  command <<<
python <<CODE
import sys
target_correct = "~{target_intervals_name}".startswith("~{bait_set}.")
bait_correct = "~{bait_intervals_name}".startswith("~{bait_set}.")
if not target_correct and not bait_correct:
  print("Bait and target intervals do not match the bait_set.")
  open("set_mismatch.txt", 'w').write("1")
  sys.stderr.write("1")
elif not target_correct:
  print("Target intervals do not match the bait_set.")
  open("set_mismatch.txt", 'w').write("1")
  sys.stderr.write("1")
elif not bait_correct:
  print("Bait intervals do not match the bait_set.")
  open("set_mismatch.txt", 'w').write("1")
  sys.stderr.write("1")
else:
  print("bait_set matches the provided intervals.")
  open("set_mismatch.txt", 'w').write("0")
CODE
  >>>

  output {
    String mismatch_message = read_string(stdout())
    Int bait_mismatch = read_int("set_mismatch.txt")
  }

  runtime {
    docker: "python:3"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 1
    maxRetries: 0
    failOnStderr: fail_task
  }
}

## ============================================================================
## Tasks from MakeCallsFromConsensus.wdl
## ============================================================================

# If pipeline is run in tumor/normal mode, ensure that they share the same fingerprint
task CrosscheckFingerprints {
   input {
      String bloodbiopsydocker
      Int preemptible
      Int disk_pad
      File input_bam
      File? second_input_bam
      File haplotype_map
      Int? memory

      Int disk_size = 20 + disk_pad
      Int mem = select_first([memory, 5])
   }
   command {
      set -e

      gatk CrosscheckFingerprints \
         --I ~{input_bam} \
         ~{"--I " + second_input_bam} \
         --EXIT_CODE_WHEN_MISMATCH 1 \
         --CROSSCHECK_BY SAMPLE \
         --EXPECT_ALL_GROUPS_TO_MATCH true \
         --HAPLOTYPE_MAP ~{haplotype_map} \
         --OUTPUT fingerprint_metrics.txt
    }
    runtime {
       docker: bloodbiopsydocker
       memory: mem + " GB"
       disks: "local-disk " + disk_size + " HDD"
       preemptible: preemptible
    }
    output {
       File fingerprint_metrics = "fingerprint_metrics.txt"
    }
}

task FingerprintsResult{
   input{
      File fingerprint_metrics
      String bloodbiopsydocker
      Int preemptible
      Int disk_pad
      Int disk_size = 20 + disk_pad
      Int? memory
      Int mem = select_first([memory, 5])
   }
   command{
      # Extract output from the fingerprint matrix

      python3 <<CODE
      import pandas as pd

      fps_mtx = pd.read_csv("${fingerprint_metrics}",sep='\t',skiprows=6)
      # Get a brief output from crosscheck-result.txt
      # if all samples matched with each other
      # it means they are very likely from the same individual
      with open('crosscheck-result.txt', 'w') as f:
         if 'EXPECTED_MATCH' in fps_mtx['RESULT'].unique():
               f.write('1')
         else:
               f.write('0')
      CODE
   }
   runtime {
       docker: bloodbiopsydocker
       memory: mem + " GB"
       disks: "local-disk " + disk_size + " HDD"
       preemptible: preemptible
    }
   output{
       Int expected_match = read_int('crosscheck-result.txt')
   }
}


# Collect Sequencing Artifact Metrics
# This is generally a good set of metrics to collect, but
# we will also need this later for the variant filtration
# step in FilterByOrientationBias
task CollectSequencingArtifactMetrics {
   input {
      String bloodbiopsydocker
      File bam_file
      File bam_idx
      File reference
      File reference_idx
      String? basename
      Int disk_pad
      Int preemptible
      Int ref_size = ceil(size(reference, "GB") + size(reference_idx, "GB"))
      Int? memory
      Int disk_size = 500
   }
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e

      gatk --java-options "-Xmx4G" \
         CollectSequencingArtifactMetrics \
         --I ${bam_file} \
         --O "${basename}.gatk" \
         --R ${reference} \
         --VALIDATION_STRINGENCY LENIENT
   }
   runtime {
      docker: bloodbiopsydocker
      memory: mem + " GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: preemptible
   }
   output {
      File pre_adapter_metrics = "${basename}.gatk.pre_adapter_detail_metrics"
   }
}


task RunMappingFilter {

   input {
      String bloodbiopsydocker
      String basename
      File vcf_file
      File vcf_idx
      File reference
      File reference_idx
      File reference_dict
      String mapping_filter_python_script
      String blastn_path
      File blastdb_nhr
      File blastdb_nin
      File blastdb_nsq
      Int preemptible
      Int disk_pad
      Int? memory
   }

   Int ref_size = ceil(size(reference, "GB") + size(reference_idx, "GB") + size(reference_dict, "GB"))
   Int blast_size = ceil(size(blastdb_nhr, "GB") + size(blastdb_nin, "GB") + size(blastdb_nsq, "GB"))
   Int mem = select_first([memory, 10])
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ref_size + blast_size + disk_pad

   command {

      set -e

      python2.7 ${mapping_filter_python_script} \
         --vcf ${vcf_file} \
         --outfile ${basename}.filtered.vcf \
         --reference_fasta ${reference} \
         --blastn ${blastn_path}

      bgzip "${basename}.filtered.vcf"
      tabix "${basename}.filtered.vcf.gz"

   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: preemptible
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
   }
   output {
      File map_filtered_vcf = "${basename}.filtered.vcf.gz"
      File map_filtered_vcf_idx = "${basename}.filtered.vcf.gz.tbi"
   }
}


# This task filters out anything not included in the filter_vcf_file
task VariantFiltration {

   input {
      String bloodbiopsydocker
      String basename
      File? gatk_override
      String mapping_filter_name
      File duplex_vcf_file
      File duplex_vcf_idx
      File filter_vcf_file
      File filter_vcf_idx
      File reference_fasta
      File reference_fasta_idx
      File reference_dict
      Int preemptible
      Boolean filter_not_in_mask
      Int? memory
      Int disk_pad
   }

   Int mem = select_first([memory, 5])
   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(duplex_vcf_file, "GB") * 1.25) + ceil(size(filter_vcf_file, "GB") * 1.25) + ref_size + disk_pad

   command {
      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      # Filter out variant with POPAF < 3.0 (this is phred scaled), which is
      # equivalent to filtering population allele frequency > 0.001.
      gatk VariantFiltration \
         -R ${reference_fasta} \
         -V ${duplex_vcf_file} \
         -O popaf.filtered.vcf.gz \
         --filter-expression "POPAF < 3.0" \
         --filter-name germline

      # Apply mapping filter
      gatk VariantFiltration \
         -R ${reference_fasta} \
         -V popaf.filtered.vcf.gz \
         -O ${basename}.${mapping_filter_name}.vcf.gz \
         --mask ${filter_vcf_file} \
         --filter-not-in-mask ${filter_not_in_mask} \
         --mask-name ${mapping_filter_name}

   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: preemptible
   }

   output {
      File output_vcf = "${basename}.${mapping_filter_name}.vcf.gz"
      File output_vcf_idx = "${basename}.${mapping_filter_name}.vcf.gz.tbi"
   }

}


task SplitVCFs {

   input {
      String bloodbiopsydocker
      String basename
      File? gatk_override
      File vcf_file
      File vcf_idx
      File reference_fasta
      File reference_fasta_idx
      File reference_dict
      Int preemptible
      Int? memory
      Int disk_pad
   }

   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ceil(size(vcf_idx)) + ref_size + disk_pad
   Int mem = select_first([memory, 16])

   command {

      set -e
      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # Include all snp and mnp calls in the snp vcf
      gatk --java-options "-Xmx15G" SelectVariants \
         -V ${vcf_file} \
         -O "${basename}.snp.vcf.gz" \
         -R ${reference_fasta} \
         -select-type SNP \
         -select-type MNP

      #passing only SNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -O "${basename}.snp.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered \
         -select-type SNP

      passing_SNP="$(gatk CountVariants -V ${basename}.snp.passing.vcf.gz | tail -1)"
      echo "$passing_SNP" > passing_snp.txt

      #passing only MNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -O "${basename}.mnp.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered \
         -select-type MNP

      passing_MNP="$(gatk CountVariants -V ${basename}.mnp.passing.vcf.gz | tail -1)"
      echo "$passing_MNP" > passing_mnp.txt

      #filtered only SNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -XL "${basename}.snp.passing.vcf.gz" \
         -O "${basename}.snp.filtered.vcf.gz" \
         -select-type SNP \
         -R ${reference_fasta}

      filtered_SNP="$(gatk CountVariants -V ${basename}.snp.filtered.vcf.gz | tail -1)"
      echo "$filtered_SNP" > filtered_snp.txt

      #filtered only MNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -XL "${basename}.mnp.passing.vcf.gz" \
         -O "${basename}.mnp.filtered.vcf.gz" \
         -select-type MNP \
         -R ${reference_fasta}

      filtered_MNP="$(gatk CountVariants -V ${basename}.mnp.filtered.vcf.gz | tail -1)"
      echo "$filtered_MNP" > filtered_mnp.txt

      gatk --java-options "-Xmx15G" SelectVariants \
         -V ${vcf_file} \
         -O "${basename}.indel.vcf.gz" \
         -R ${reference_fasta} \
         -select-type INDEL \
         -select-type MIXED

      #passing only
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.indel.vcf.gz" \
         -O "${basename}.indel.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered

      passing_INDEL="$(gatk CountVariants -V ${basename}.indel.passing.vcf.gz | tail -1)"
      echo "$passing_INDEL" > passing_indel.txt

      #filtered only
      bcftools isec -C -O v -o ${basename}.indel.filtered.vcf -w1 ${basename}.indel.vcf.gz ${basename}.indel.passing.vcf.gz
      filtered_INDEL="$(gatk CountVariants -V ${basename}.indel.filtered.vcf | tail -1)"
      echo "$filtered_INDEL" > filtered_indel.txt

   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: preemptible
   }
   output {

      File output_snp_vcf = "${basename}.snp.vcf.gz"
      File output_snp_vcf_idx = "${basename}.snp.vcf.gz.tbi"

      Int passing_SNP = read_int("passing_snp.txt")
      Int filtered_SNP = read_int("filtered_snp.txt")

      Int passing_MNP = read_int("passing_mnp.txt")
      Int filtered_MNP = read_int("filtered_mnp.txt")

      File output_indel_vcf = "${basename}.indel.vcf.gz"
      File output_indel_vcf_idx = "${basename}.indel.vcf.gz.tbi"

      Int passing_INDEL = read_int("passing_indel.txt")
      Int filtered_INDEL = read_int("filtered_indel.txt")

   }
}

#creates an IGV session
#given a list of IGV compatible files (as strings)
#reference is either "hg19" or "hg38"
task GenerateIGVSession {

   input {
      String bloodbiopsydocker
      Array[String] input_files
      String reference_version
      String file_name
      Array[String]? input_names
      Int preemptible
   }

   Array[String] input_names_prefix = if defined(input_names) then prefix('-n ', select_first([input_names])) else []

   command {
      bash /usr/writeIGV.sh ${reference_version} ${sep=" " input_files} ${sep=" " input_names_prefix}  > "${file_name}.xml"
   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: preemptible
   }
   output {
      File igv_session = "${file_name}.xml"
   }
}

task ExtractContam {

   input {
      String bloodbiopsydocker
      File contamination_table
      Int preemptible
   }

   command {
      tail -n1 ~{contamination_table} | cut -f2 > fraction_contamination.txt
   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: preemptible
   }
   output {
      Float fracContam=read_float("fraction_contamination.txt")
   }
}

task VariantsToTable {
   input {
      String bloodbiopsydocker
      File vcf
      File vcf_idx
      String output_name
      File? gatk_override
      Int? memory
      Int disk_pad
      Int preemptible
   }

   Int mem = if defined(memory) then memory * 1000 else 3500
   Int compute_mem = mem - 500
   Int disk_size = ceil(size(vcf, "GB") * 1.5 ) + ceil(size(vcf_idx, "GB")) + disk_pad

   command {

      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx${compute_mem}m" VariantsToTable \
         -V ${vcf} \
         -F CHROM -F POS -F ID -F REF -F ALT -F QUAL \
         -F FILTER -ASGF AD -ASGF AF -GF NCount -ASGF SB \
         -SMA \
         -O "${output_name}.table"

   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + 500 + " HDD"
      memory: mem + " MB"
      preemptible: preemptible
   }
   output {
      File output_table = "${output_name}.table"
   }
}

## ============================================================================
## Tasks from mutect2.wdl
## ============================================================================

task CramToBam {
    input {
      File ref_fasta
      File ref_fai
      File ref_dict
      #cram and crai must be optional since Normal cram is optional
      File? cram
      File? crai
      String name
      Int disk_size
      Int? mem
    }

    Int machine_mem = if defined(mem) then mem * 1000 else 6000

    #Calls samtools view to do the conversion
    command {
        #Set -e and -o says if any command I run fails in this script, make sure to return a failure
        set -e
        set -o pipefail

        samtools view -h -T ~{ref_fasta} ~{cram} |
            samtools view -b -o ~{name}.bam -
        samtools index -b ~{name}.bam
        mv ~{name}.bam.bai ~{name}.bai
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_bam = "~{name}.bam"
        File output_bai = "~{name}.bai"
    }
}

task SplitIntervals {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      Int scatter_count
      String? split_intervals_extra_args

      # runtime
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        mkdir interval-files
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" SplitIntervals \
            -R ~{ref_fasta} \
            ~{"-L " + intervals} \
            -scatter ~{scatter_count} \
            -O interval-files \
            ~{split_intervals_extra_args}
        cp interval-files/*.interval_list .
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task M2 {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File tumor_bam
      File tumor_bai
      File? normal_bam
      File? normal_bai
      File? pon
      File? pon_idx
      File? gnomad
      File? gnomad_idx
      String? m2_extra_args
      String? getpileupsummaries_extra_args
      Boolean? make_bamout
      Boolean? run_ob_filter
      Boolean compress
      File? gga_vcf
      File? gga_vcf_idx
      File? variants_for_contamination
      File? variants_for_contamination_idx

      File? gatk_override

      # runtime
      String gatk_docker
      Float? mem
      Int? preemptible
      Int? max_retries
      Int? disk_space
      Int? cpu
      Boolean use_ssd = false
    }

    String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then ceil(mem * 1000) else 3500
    Int command_mem = machine_mem - 500

    # NOTE: localization_optional removed for the Manifold/AWS (S3) backend.
    # GATK NIO can stream from gs:// but not s3://, so these inputs must be
    # localized to disk; otherwise GATK receives the raw s3:// URI and fails
    # with "fasta file ... does not exist".

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        touch f1r2.tar.gz
        echo "" > normal_name.txt

        gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{tumor_bam} -O tumor_name.txt -encode
        tumor_command_line="-I ~{tumor_bam} -tumor `cat tumor_name.txt`"

        if [[ ! -z "~{normal_bam}" ]]; then
            gatk --java-options "-Xmx~{command_mem}m" GetSampleName -R ~{ref_fasta} -I ~{normal_bam} -O normal_name.txt -encode
            normal_command_line="-I ~{normal_bam} -normal `cat normal_name.txt`"
        fi

        gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
            -R ~{ref_fasta} \
            $tumor_command_line \
            $normal_command_line \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            ~{"-L " + intervals} \
            ~{"--alleles " + gga_vcf} \
            -O "~{output_vcf}" \
            ~{true='--bam-output bamout.bam' false='' make_bamout} \
            ~{true='--f1r2-tar-gz f1r2.tar.gz' false='' run_ob_filter} \
            ~{m2_extra_args}

        m2_exit_code=$?

        ### GetPileupSummaries

        # If the variants for contamination and the intervals for this scatter don't intersect, GetPileupSummaries
        # throws an error.  However, there is nothing wrong with an empty intersection for our purposes; it simply doesn't
        # contribute to the merged pileup summaries that we create downstream.  We implement this by with array outputs.
        # If the tool errors, no table is created and the glob yields an empty array.
        set +e

        if [[ ! -z "~{variants_for_contamination}" ]]; then
            gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries -R ~{ref_fasta} -I ~{tumor_bam} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O tumor-pileups.table ~{getpileupsummaries_extra_args}


            if [[ ! -z "~{normal_bam}" ]]; then
                gatk --java-options "-Xmx~{command_mem}m" GetPileupSummaries -R ~{ref_fasta} -I ~{normal_bam} ~{"--interval-set-rule INTERSECTION -L " + intervals} \
                    -V ~{variants_for_contamination} -L ~{variants_for_contamination} -O normal-pileups.table ~{getpileupsummaries_extra_args}
            fi
        fi

        # the script only fails if Mutect2 itself fails
        exit $m2_exit_code
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible, 10])
        maxRetries: select_first([max_retries, 0])
        cpu: select_first([cpu, 1])
    }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        File output_bamOut = "bamout.bam"
        String tumor_sample = read_string("tumor_name.txt")
        String normal_sample = read_string("normal_name.txt")
        File stats = "~{output_stats}"
        File f1r2_counts = "f1r2.tar.gz"
        Array[File] tumor_pileups = glob("*tumor-pileups.table")
        Array[File] normal_pileups = glob("*normal-pileups.table")
    }
}

task MergeVCFs {
    input {
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      String output_name
      Boolean compress
      Runtime runtime_params
    }

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    # using MergeVcfs instead of GatherVcfs so we can create indices
    # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}
        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task MergeStats {
    input {
      Array[File]+ stats
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}


        gatk --java-options "-Xmx~{runtime_params.command_mem}m" MergeMutectStats \
            -stats ~{sep=" -stats " stats} -O merged.stats
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_stats = "merged.stats"
    }
}

task MergePileupSummaries {
    input {
      Array[File] input_tables
      String output_name
      File ref_dict
      Runtime runtime_params
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" GatherPileupSummaries \
        --sequence-dictionary ~{ref_dict} \
        -I ~{sep=' -I ' input_tables} \
        -O ~{output_name}.tsv
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File merged_table = "~{output_name}.tsv"
    }
}

# Learning step of the orientation bias mixture model, which is the recommended orientation bias filter as of September 2018
task LearnReadOrientationModel {
    input {
      Array[File] f1r2_tar_gz
      Runtime runtime_params
      Int? mem  #override memory
    }

    Int machine_mem = select_first([mem, runtime_params.machine_mem])
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{command_mem}m" LearnReadOrientationModel \
            -I ~{sep=" -I " f1r2_tar_gz} \
            -O "artifact-priors.tar.gz"
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File artifact_prior_table = "artifact-priors.tar.gz"
    }

}

task CalculateContamination {
    input {
      String? intervals
      File tumor_pileups
      File? normal_pileups
      Runtime runtime_params
    }

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" CalculateContamination -I ~{tumor_pileups} \
        -O contamination.table --tumor-segmentation segments.table ~{"-matched " + normal_pileups}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + runtime_params.disk + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File contamination_table = "contamination.table"
        File maf_segments = "segments.table"
    }
}

task Filter {
    input {
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File unfiltered_vcf
      File unfiltered_vcf_idx
      String output_name
      Boolean compress
      File? mutect_stats
      File? artifact_priors_tar_gz
      File? contamination_table
      File? maf_segments
      String? m2_extra_filtering_args

      Runtime runtime_params
      Int? disk_space
    }

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    # NOTE: localization_optional removed for the Manifold/AWS (S3) backend
    # (GATK cannot stream from s3://; inputs must be localized to disk).

    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

        gatk --java-options "-Xmx~{runtime_params.command_mem}m" FilterMutectCalls -V ~{unfiltered_vcf} \
            -R ~{ref_fasta} \
            -O ~{output_vcf} \
            ~{"--contamination-table " + contamination_table} \
            ~{"--tumor-segmentation " + maf_segments} \
            ~{"--ob-priors " + artifact_priors_tar_gz} \
            ~{"-stats " + mutect_stats} \
            --filtering-stats filtering.stats \
            ~{m2_extra_filtering_args}
    }

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "filtering.stats"
    }
}

task Funcotate {
     input {
       File ref_fasta
       File ref_fai
       File ref_dict
       File input_vcf
       File input_vcf_idx
       String reference_version
       String output_file_base_name
       String output_format
       Boolean compress
       Boolean use_gnomad
       # This should be updated when a new version of the data sources is released
       # TODO: Make this dynamically chosen in the command.
       File? data_sources_tar_gz = "gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz"
       String? control_id
       String? case_id
       String? sequencing_center
       String? sequence_source
       String? transcript_selection_mode
       File? transcript_selection_list
       Array[String]? annotation_defaults
       Array[String]? annotation_overrides
       Array[String]? funcotator_excluded_fields
       Boolean? filter_funcotations
       File? interval_list

       String? extra_args

       # ==============
       Runtime runtime_params
       Int? disk_space   #override to request more disk than default small task params

       # You may have to change the following two parameter values depending on the task requirements
       Int default_ram_mb = 3000
       # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
       Int default_disk_space_gb = 100
     }

     # ==============
     # Process input args:
     String output_maf = output_file_base_name + ".maf"
     String output_maf_index = output_maf + ".idx"
     String output_vcf = output_file_base_name + if compress then ".vcf.gz" else ".vcf"
     String output_vcf_idx = output_vcf +  if compress then ".tbi" else ".idx"
     String output_file = if output_format == "MAF" then output_maf else output_vcf
     String output_file_index = if output_format == "MAF" then output_maf_index else output_vcf_idx
     String transcript_selection_arg = if defined(transcript_selection_list) then " --transcript-list " else ""
     String annotation_def_arg = if defined(annotation_defaults) then " --annotation-default " else ""
     String annotation_over_arg = if defined(annotation_overrides) then " --annotation-override " else ""
     String filter_funcotations_args = if defined(filter_funcotations) && (filter_funcotations) then " --remove-filtered-variants " else ""
     String excluded_fields_args = if defined(funcotator_excluded_fields) then " --exclude-field " else ""
     String interval_list_arg = if defined(interval_list) then " -L " else ""
     String extra_args_arg = select_first([extra_args, ""])

     String dollar = "$"

     # NOTE: localization_optional removed for the Manifold/AWS (S3) backend
     # (GATK cannot stream from s3://; inputs must be localized to disk).

     command <<<
         set -e
         export GATK_LOCAL_JAR=~{default="/root/gatk.jar" runtime_params.gatk_override}

         # Extract our data sources:
         echo "Extracting data sources zip file..."
         mkdir datasources_dir
         tar zxvf ~{data_sources_tar_gz} -C datasources_dir --strip-components 1
         DATA_SOURCES_FOLDER="$PWD/datasources_dir"

         # Handle gnomAD:
         if ~{use_gnomad} ; then
             echo "Enabling gnomAD..."
             for potential_gnomad_gz in gnomAD_exome.tar.gz gnomAD_genome.tar.gz ; do
                 if [[ -f ~{dollar}{DATA_SOURCES_FOLDER}/~{dollar}{potential_gnomad_gz} ]] ; then
                     cd ~{dollar}{DATA_SOURCES_FOLDER}
                     tar -zvxf ~{dollar}{potential_gnomad_gz}
                     cd -
                 else
                     echo "ERROR: Cannot find gnomAD folder: ~{dollar}{potential_gnomad_gz}" 1>&2
                     false
                 fi
             done
         fi

         # Run Funcotator:
         gatk --java-options "-Xmx~{runtime_params.command_mem}m" Funcotator \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ~{reference_version} \
             --output-file-format ~{output_format} \
             -R ~{ref_fasta} \
             -V ~{input_vcf} \
             -O ~{output_file} \
             ~{interval_list_arg} ~{default="" interval_list} \
             --annotation-default normal_barcode:~{default="Unknown" control_id} \
             --annotation-default tumor_barcode:~{default="Unknown" case_id} \
             --annotation-default Center:~{default="Unknown" sequencing_center} \
             --annotation-default source:~{default="Unknown" sequence_source} \
             ~{"--transcript-selection-mode " + transcript_selection_mode} \
             ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \
             ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \
             ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \
             ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \
             ~{filter_funcotations_args} \
             ~{extra_args_arg}
         # Make sure we have a placeholder index for MAF files so this workflow doesn't fail:
         if [[ "~{output_format}" == "MAF" ]] ; then
            touch ~{output_maf_index}
         fi
     >>>

    runtime {
        docker: runtime_params.gatk_docker
        bootDiskSizeGb: runtime_params.boot_disk_size
        memory: runtime_params.machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, runtime_params.disk]) + " HDD"
        preemptible: runtime_params.preemptible
        maxRetries: runtime_params.max_retries
        cpu: runtime_params.cpu
    }

     output {
         File funcotated_output_file = "~{output_file}"
         File funcotated_output_file_index = "~{output_file_index}"
     }
}
