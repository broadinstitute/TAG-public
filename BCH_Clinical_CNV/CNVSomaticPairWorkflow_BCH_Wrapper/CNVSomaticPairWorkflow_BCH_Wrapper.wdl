version 1.0

import "./subworkflows/cnv_somatic_pair_workflow.wdl" as CNV_WDL
import "/BCH_Clinical_CNV/PlotCNVQC/plotCNVQC.wdl" as CNV_QC

workflow CNVSomaticPairWorkflow_BCH_Wrapper {
	input {
      ##################################
      #### required basic arguments ####
      ##################################
      File common_sites
      File intervals
      File? blacklist_intervals
      File tumor_bam
      File tumor_bam_idx
      File? normal_bam
      File? normal_bam_idx
      File read_count_pon
      File ref_fasta_dict
      File ref_fasta_fai
      File ref_fasta
      String gatk_docker

      ##################################
      #### optional basic arguments ####
      ##################################
       # For running oncotator
      Boolean? is_run_oncotator
       # For running funcotator
      Boolean? is_run_funcotator

      File? gatk4_jar_override
      Int? preemptible_attempts
      # Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk

      # Required if BAM/CRAM is in a requester pays bucket
      String? gcs_project_for_requester_pays

      ####################################################
      #### optional arguments for PreprocessIntervals ####
      ####################################################
      Int? padding
      Int? bin_length
      Int? mem_gb_for_preprocess_intervals

      ##############################################
      #### optional arguments for CollectCounts ####
      ##############################################
      String? collect_counts_format
      Int? mem_gb_for_collect_counts

      #####################################################
      #### optional arguments for CollectAllelicCounts ####
      #####################################################
      String? minimum_base_quality
      Int? mem_gb_for_collect_allelic_counts

      ##################################################
      #### optional arguments for DenoiseReadCounts ####
      ##################################################
      Int? number_of_eigensamples
      Int? mem_gb_for_denoise_read_counts

      ##############################################
      #### optional arguments for ModelSegments ####
      ##############################################
      Int? max_num_segments_per_chromosome
      Int? min_total_allele_count
      Int? min_total_allele_count_normal
      Float? genotyping_homozygous_log_ratio_threshold
      Float? genotyping_base_error_rate
      Float? kernel_variance_copy_ratio
      Float? kernel_variance_allele_fraction
      Float? kernel_scaling_allele_fraction
      Int? kernel_approximation_dimension
      Array[Int]+? window_sizes = [8, 16, 32, 64, 128, 256]
      Float? num_changepoints_penalty_factor
      Float? minor_allele_fraction_prior_alpha
      Int? num_samples_copy_ratio
      Int? num_burn_in_copy_ratio
      Int? num_samples_allele_fraction
      Int? num_burn_in_allele_fraction
      Float? smoothing_threshold_copy_ratio
      Float? smoothing_threshold_allele_fraction
      Int? max_num_smoothing_iterations
      Int? num_smoothing_iterations_per_fit
      Int? mem_gb_for_model_segments

      ######################################################
      #### optional arguments for CallCopyRatioSegments ####
      ######################################################
      Float? neutral_segment_copy_ratio_lower_bound
      Float? neutral_segment_copy_ratio_upper_bound
      Float? outlier_neutral_segment_copy_ratio_z_score_threshold
      Float? calling_copy_ratio_z_score_threshold
      Int? mem_gb_for_call_copy_ratio_segments

      #########################################
      #### optional arguments for plotting ####
      #########################################
      Int? minimum_contig_length
      # If maximum_copy_ratio = Infinity, the maximum copy ratio will be automatically determined
      String? maximum_copy_ratio
      Float? point_size_copy_ratio
      Float? point_size_allele_fraction
      Int? mem_gb_for_plotting

      ##########################################
      #### optional arguments for Oncotator ####
      ##########################################
      String? additional_args_for_oncotator
      String? oncotator_docker
      Int? mem_gb_for_oncotator
      Int? boot_disk_space_gb_for_oncotator

      ##################################################
      #### optional arguments for FuncotateSegments ####
      ##################################################
      String? additional_args_for_funcotator
      String? funcotator_ref_version
      Int? mem_gb_for_funcotator
      File? funcotator_transcript_selection_list
      File? funcotator_data_sources_tar_gz
      String? funcotator_transcript_selection_mode
      Array[String]? funcotator_annotation_defaults
      Array[String]? funcotator_annotation_overrides
      Array[String]? funcotator_excluded_fields
      Boolean? funcotator_is_removing_untared_datasources
      Int? funcotator_disk_space_gb
      Boolean? funcotator_use_ssd
      Int? funcotator_cpu


      #############################################
      #### optional arguments for clinical CNV ####
      #############################################
      File evaluate_cnv_script
      String python_docker = "python:3"
      String pairID
      String quicvizDocker = "us.gcr.io/tag-team-160914/cmi_quicviz:0.4.3"
	}

	call CNV_WDL.CNVSomaticPairWorkflow as CNVSomaticPairWorkflow{
		input :
			common_sites = common_sites,
			intervals = intervals,
			blacklist_intervals = blacklist_intervals,
            tumor_bam = tumor_bam,
            tumor_bam_idx = tumor_bam_idx,
            normal_bam = normal_bam,
            normal_bam_idx = normal_bam_idx,
            read_count_pon = read_count_pon,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_fai = ref_fasta_fai,
            ref_fasta = ref_fasta,
            gatk_docker = gatk_docker,
            is_run_oncotator = is_run_oncotator,
            is_run_funcotator = is_run_funcotator,
            gatk4_jar_override = gatk4_jar_override,
            preemptible_attempts = preemptible_attempts,
            emergency_extra_disk = emergency_extra_disk,
            gcs_project_for_requester_pays = gcs_project_for_requester_pays,
            padding = padding,
            bin_length = bin_length,
            mem_gb_for_preprocess_intervals = mem_gb_for_preprocess_intervals,
            collect_counts_format = collect_counts_format,
            mem_gb_for_collect_counts = mem_gb_for_collect_counts,
            minimum_base_quality = minimum_base_quality,
            mem_gb_for_collect_allelic_counts = mem_gb_for_collect_allelic_counts,
            number_of_eigensamples = number_of_eigensamples,
            mem_gb_for_denoise_read_counts = mem_gb_for_denoise_read_counts,
            max_num_segments_per_chromosome = max_num_segments_per_chromosome,
            min_total_allele_count = min_total_allele_count,
            min_total_allele_count_normal = min_total_allele_count_normal,
            genotyping_homozygous_log_ratio_threshold = genotyping_homozygous_log_ratio_threshold,
            genotyping_base_error_rate = genotyping_base_error_rate,
            kernel_variance_copy_ratio = kernel_variance_copy_ratio,
            kernel_variance_allele_fraction = kernel_variance_allele_fraction,
            kernel_scaling_allele_fraction = kernel_scaling_allele_fraction,
            kernel_approximation_dimension = kernel_approximation_dimension,
            window_sizes = window_sizes ,
            num_changepoints_penalty_factor = num_changepoints_penalty_factor,
            minor_allele_fraction_prior_alpha = minor_allele_fraction_prior_alpha,
            num_samples_copy_ratio = num_samples_copy_ratio,
            num_burn_in_copy_ratio = num_burn_in_copy_ratio,
            num_samples_allele_fraction = num_samples_allele_fraction,
            num_burn_in_allele_fraction = num_burn_in_allele_fraction,
            smoothing_threshold_copy_ratio = smoothing_threshold_copy_ratio,
            smoothing_threshold_allele_fraction = smoothing_threshold_allele_fraction,
            max_num_smoothing_iterations = max_num_smoothing_iterations,
            num_smoothing_iterations_per_fit = num_smoothing_iterations_per_fit,
            mem_gb_for_model_segments = mem_gb_for_model_segments,
            neutral_segment_copy_ratio_lower_bound = neutral_segment_copy_ratio_lower_bound,
            neutral_segment_copy_ratio_upper_bound = neutral_segment_copy_ratio_upper_bound,
            outlier_neutral_segment_copy_ratio_z_score_threshold = outlier_neutral_segment_copy_ratio_z_score_threshold,
            calling_copy_ratio_z_score_threshold = calling_copy_ratio_z_score_threshold,
            mem_gb_for_call_copy_ratio_segments = mem_gb_for_call_copy_ratio_segments,
            minimum_contig_length = minimum_contig_length,
            maximum_copy_ratio = maximum_copy_ratio,
            point_size_copy_ratio = point_size_copy_ratio,
            point_size_allele_fraction = point_size_allele_fraction,
            mem_gb_for_plotting = mem_gb_for_plotting,
            additional_args_for_oncotator = additional_args_for_oncotator,
            oncotator_docker = oncotator_docker,
            mem_gb_for_oncotator = mem_gb_for_oncotator,
            boot_disk_space_gb_for_oncotator = boot_disk_space_gb_for_oncotator,
            additional_args_for_funcotator = additional_args_for_funcotator,
            funcotator_ref_version = funcotator_ref_version,
            mem_gb_for_funcotator = mem_gb_for_funcotator,
            funcotator_transcript_selection_list = funcotator_transcript_selection_list,
            funcotator_data_sources_tar_gz = funcotator_data_sources_tar_gz,
            funcotator_transcript_selection_mode = funcotator_transcript_selection_mode,
            funcotator_annotation_defaults = funcotator_annotation_defaults,
            funcotator_annotation_overrides = funcotator_annotation_overrides,
            funcotator_excluded_fields = funcotator_excluded_fields,
            funcotator_is_removing_untared_datasources = funcotator_is_removing_untared_datasources,
            funcotator_disk_space_gb = funcotator_disk_space_gb,
            funcotator_use_ssd = funcotator_use_ssd,
            funcotator_cpu = funcotator_cpu
	}
	call CallCNVPassFail {
        input: 
            tumor_seg_file = CNVSomaticPairWorkflow.called_copy_ratio_segments_tumor,
            normal_seg_file = select_first([CNVSomaticPairWorkflow.called_copy_ratio_segments_normal]),
            tumor_MAD_value = CNVSomaticPairWorkflow.denoised_MAD_value_tumor,
            normal_MAD_value = select_first([CNVSomaticPairWorkflow.denoised_MAD_value_normal]),
            cnv_eval_script = evaluate_cnv_script,
            python_docker = python_docker
    }

	call QUICviz {
        input:
            pairID = pairID,
            tumor_sample_id = CNVSomaticPairWorkflow.read_counts_entity_id_tumor,
            normal_sample_id = select_first([CNVSomaticPairWorkflow.read_counts_entity_id_normal]),
            quicvizDocker = quicvizDocker,
            allelicCountsNormal = select_first([CNVSomaticPairWorkflow.allelic_counts_normal]),
            allelicCountsTumor = CNVSomaticPairWorkflow.allelic_counts_tumor,
            denoisedCopyRatiosNormal = select_first([CNVSomaticPairWorkflow.denoised_copy_ratios_normal]),
            denoisedCopyRatiosTumor = CNVSomaticPairWorkflow.denoised_copy_ratios_tumor,
            calledCopyRatioSegTumor = CNVSomaticPairWorkflow.called_copy_ratio_segments_tumor,
            oncotatedCalledTumor = CNVSomaticPairWorkflow.oncotated_called_file_tumor
    }

    call CNV_QC.plotQCTask as plotQCTask {
        input:
            called_copy_ratio_segments_tumor = CNVSomaticPairWorkflow.called_copy_ratio_segments_tumor,
            called_copy_ratio_segments_normal = select_first([CNVSomaticPairWorkflow.called_copy_ratio_segments_normal]),
            ref_fasta_index = ref_fasta_fai,
            output_basename = pairID,
            denoised_MAD_value_normal = select_first([CNVSomaticPairWorkflow.denoised_MAD_value_normal]),
            denoised_MAD_value_tumor = CNVSomaticPairWorkflow.denoised_MAD_value_tumor
    }

	output {
		File preprocessed_intervals = CNVSomaticPairWorkflow.preprocessed_intervals
        File read_counts_entity_id_tumor = CNVSomaticPairWorkflow.read_counts_entity_id_tumor
        File read_counts_tumor = CNVSomaticPairWorkflow.read_counts_tumor
        File allelic_counts_entity_id_tumor = CNVSomaticPairWorkflow.allelic_counts_entity_id_tumor
        File allelic_counts_tumor = CNVSomaticPairWorkflow.allelic_counts_tumor
        File denoised_copy_ratios_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_tumor
        File standardized_copy_ratios_tumor = CNVSomaticPairWorkflow.standardized_copy_ratios_tumor
        File het_allelic_counts_tumor = CNVSomaticPairWorkflow.het_allelic_counts_tumor
        File normal_het_allelic_counts_tumor = CNVSomaticPairWorkflow.normal_het_allelic_counts_tumor
        File copy_ratio_only_segments_tumor = CNVSomaticPairWorkflow.copy_ratio_only_segments_tumor
        File copy_ratio_legacy_segments_tumor = CNVSomaticPairWorkflow.copy_ratio_legacy_segments_tumor
        File allele_fraction_legacy_segments_tumor = CNVSomaticPairWorkflow.allele_fraction_legacy_segments_tumor
        File modeled_segments_begin_tumor = CNVSomaticPairWorkflow.modeled_segments_begin_tumor
        File copy_ratio_parameters_begin_tumor = CNVSomaticPairWorkflow.copy_ratio_parameters_begin_tumor
        File allele_fraction_parameters_begin_tumor = CNVSomaticPairWorkflow.allele_fraction_parameters_begin_tumor
        File modeled_segments_tumor = CNVSomaticPairWorkflow.modeled_segments_tumor
        File copy_ratio_parameters_tumor = CNVSomaticPairWorkflow.copy_ratio_parameters_tumor
        File allele_fraction_parameters_tumor = CNVSomaticPairWorkflow.allele_fraction_parameters_tumor
        File called_copy_ratio_segments_tumor = CNVSomaticPairWorkflow.called_copy_ratio_segments_tumor
        File called_copy_ratio_legacy_segments_tumor = CNVSomaticPairWorkflow.called_copy_ratio_legacy_segments_tumor
        File denoised_copy_ratios_plot_tumor = CNVSomaticPairWorkflow.denoised_copy_ratios_plot_tumor
        File standardized_MAD_tumor = CNVSomaticPairWorkflow.standardized_MAD_tumor
        Float standardized_MAD_value_tumor = CNVSomaticPairWorkflow.standardized_MAD_value_tumor
        File denoised_MAD_tumor = CNVSomaticPairWorkflow.denoised_MAD_tumor
        Float denoised_MAD_value_tumor = CNVSomaticPairWorkflow.denoised_MAD_value_tumor
        File delta_MAD_tumor = CNVSomaticPairWorkflow.delta_MAD_tumor
        Float delta_MAD_value_tumor = CNVSomaticPairWorkflow.delta_MAD_value_tumor
        File scaled_delta_MAD_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_tumor
        Float scaled_delta_MAD_value_tumor = CNVSomaticPairWorkflow.scaled_delta_MAD_value_tumor
        File modeled_segments_plot_tumor = CNVSomaticPairWorkflow.modeled_segments_plot_tumor

        File? read_counts_entity_id_normal = CNVSomaticPairWorkflow.read_counts_entity_id_normal
        File? read_counts_normal = CNVSomaticPairWorkflow.read_counts_normal
        File? allelic_counts_entity_id_normal = CNVSomaticPairWorkflow.allelic_counts_entity_id_normal
        File? allelic_counts_normal = CNVSomaticPairWorkflow.allelic_counts_normal
        File? denoised_copy_ratios_normal = CNVSomaticPairWorkflow.denoised_copy_ratios_normal
        File? standardized_copy_ratios_normal = CNVSomaticPairWorkflow.standardized_copy_ratios_normal
        File? het_allelic_counts_normal = CNVSomaticPairWorkflow.het_allelic_counts_normal
        File? normal_het_allelic_counts_normal = CNVSomaticPairWorkflow.normal_het_allelic_counts_normal
        File? copy_ratio_only_segments_normal = CNVSomaticPairWorkflow.copy_ratio_only_segments_normal
        File? copy_ratio_legacy_segments_normal = CNVSomaticPairWorkflow.copy_ratio_legacy_segments_normal
        File? allele_fraction_legacy_segments_normal = CNVSomaticPairWorkflow.allele_fraction_legacy_segments_normal
        File? modeled_segments_begin_normal = CNVSomaticPairWorkflow.modeled_segments_begin_normal
        File? copy_ratio_parameters_begin_normal = CNVSomaticPairWorkflow.copy_ratio_parameters_begin_normal
        File? allele_fraction_parameters_begin_normal = CNVSomaticPairWorkflow.allele_fraction_parameters_begin_normal
        File? modeled_segments_normal = CNVSomaticPairWorkflow.modeled_segments_normal
        File? copy_ratio_parameters_normal = CNVSomaticPairWorkflow.copy_ratio_parameters_normal
        File? allele_fraction_parameters_normal = CNVSomaticPairWorkflow.allele_fraction_parameters_normal
        File? called_copy_ratio_segments_normal = CNVSomaticPairWorkflow.called_copy_ratio_segments_normal
        File? called_copy_ratio_legacy_segments_normal = CNVSomaticPairWorkflow.called_copy_ratio_legacy_segments_normal
        File? denoised_copy_ratios_plot_normal = CNVSomaticPairWorkflow.denoised_copy_ratios_plot_normal
        File? standardized_MAD_normal = CNVSomaticPairWorkflow.standardized_MAD_normal
        Float? standardized_MAD_value_normal = CNVSomaticPairWorkflow.standardized_MAD_value_normal
        File? denoised_MAD_normal = CNVSomaticPairWorkflow.denoised_MAD_normal
        Float? denoised_MAD_value_normal = CNVSomaticPairWorkflow.denoised_MAD_value_normal
        File? delta_MAD_normal = CNVSomaticPairWorkflow.delta_MAD_normal
        Float? delta_MAD_value_normal = CNVSomaticPairWorkflow.delta_MAD_value_normal
        File? scaled_delta_MAD_normal = CNVSomaticPairWorkflow.scaled_delta_MAD_normal
        Float? scaled_delta_MAD_value_normal = CNVSomaticPairWorkflow.scaled_delta_MAD_value_normal
        File? modeled_segments_plot_normal = CNVSomaticPairWorkflow.modeled_segments_plot_normal

        File oncotated_called_file_tumor = CNVSomaticPairWorkflow.oncotated_called_file_tumor
        File oncotated_called_gene_list_file_tumor = CNVSomaticPairWorkflow.oncotated_called_gene_list_file_tumor
        File funcotated_called_file_tumor = CNVSomaticPairWorkflow.funcotated_called_file_tumor
        File funcotated_called_gene_list_file_tumor = CNVSomaticPairWorkflow.funcotated_called_gene_list_file_tumor

		String tumor_cnv_pass_fail = CallCNVPassFail.tumor_cnv_pass_fail
        Int tumor_called_segments_count = CallCNVPassFail.tumor_called_segments_count
        String normal_cnv_pass_fail = CallCNVPassFail.normal_cnv_pass_fail
        Int normal_called_segments_count = CallCNVPassFail.normal_called_segments_count
        String cnv_pass_fail = CallCNVPassFail.cnv_pass_fail

        File QUICvizPDF = QUICviz.QUICvizPDF
        File GeneLevelCNV = QUICviz.GeneLevelCNV
        File AllChrPlot = QUICviz.AllChrPlot

        File CNV_QC_plots = plotQCTask.CNV_QC_plots
	}
}

task CallCNVPassFail {
    input {
      File tumor_seg_file
      File normal_seg_file
      Float tumor_MAD_value
      Float normal_MAD_value
      File cnv_eval_script

      Float normal_MAD_threshold = 0.15
      Int normal_seg_threshold = 200
      Float tumor_MAD_threshold = 0.2
      Int tumor_seg_threshold = 1000

      # Runtime parameters
      String python_docker
      Int mem_gb = 5
      Int disk_space_gb = 50
      Int? preemptible_attempts
    }

    command <<<
        python ~{cnv_eval_script} --tumor-seg ~{tumor_seg_file} --normal-seg ~{normal_seg_file} \
        --tumor-mad ~{tumor_MAD_value} --normal-mad ~{normal_MAD_value} \
        --thresholds ~{normal_MAD_threshold} ~{normal_seg_threshold} ~{tumor_MAD_threshold} ~{tumor_seg_threshold}
    >>>

    runtime {
        docker: "~{python_docker}"
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_space_gb + " HDD"
        cpu: 1
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        Int tumor_called_segments_count = read_int("tumor_seg_count.txt")
        String tumor_cnv_pass_fail = read_string("tumor_res.txt")
        Int normal_called_segments_count = read_int("normal_seg_count.txt")
        String normal_cnv_pass_fail = read_string("normal_res.txt")
        String cnv_pass_fail = read_string("cnv_res.txt")
    }
}


task QUICviz {
    input {
        String tumor_sample_id
        String normal_sample_id
        String pairID
        String tumorType
        String quicvizDocker
        File allelicCountsNormal
        File allelicCountsTumor
        File denoisedCopyRatiosNormal
        File denoisedCopyRatiosTumor
        File calledCopyRatioSegTumor
        File oncotatedCalledTumor
        File? gene_list_override
        Int memory = 16
        Int cpu = 4
        Int maxRetries = 3
    }
    command <<<
        set -e
        mkdir outputs

        echo "Input Tumor Sample: "~{tumor_sample_id}
        echo "Input Normal Sample: "~{normal_sample_id}


        Rscript /BaseImage/CMI_QUICviz/scripts/CMI_QUICviz.R \
            --sample ~{tumor_sample_id} \
            --tumor_type ~{tumorType} \
            --normal_acf ~{allelicCountsNormal} \
            --normal_cr ~{denoisedCopyRatiosNormal} \
            --tumor_acf ~{allelicCountsTumor} \
            --tumor_cr ~{denoisedCopyRatiosTumor} \
            --tumor_cr_seg ~{calledCopyRatioSegTumor} \
            --tumor_seg_oncotated ~{oncotatedCalledTumor} \
            ~{'--gene_list '+ gene_list_override} \
            --output_dir outputs/

        mv outputs/*chromosome_plots.pdf outputs/~{pairID}_chromosome_plots.pdf
        mv outputs/*gene_level_calls.csv outputs/~{pairID}_gene_level_calls.csv
        mv outputs/*_all_chr.png outputs/~{pairID}_All_chr.png

    >>>
    output {
        File QUICvizPDF = "outputs/~{pairID}_chromosome_plots.pdf"
        File GeneLevelCNV = "outputs/~{pairID}_gene_level_calls.csv"
        File AllChrPlot = "outputs/~{pairID}_All_chr.png"
    }
    runtime {
        docker: quicvizDocker
        memory: memory + " GB"
        cpu: cpu
        disks: "local-disk 100 HDD"
        maxRetries: maxRetries
    }
}
