version 1.0

import "../NeoVax/neovax_dragen_modules/AllelicCNV.wdl" as AllelicCNV
import "../NeoVax/neovax_dragen_modules/Absolute.wdl" as Absolute
import "../NeoVax/neovax_dragen_modules/PairBamQC.wdl" as PairBamQC
import "../NeoVax/neovax_dragen_modules/SomaticQC.wdl" as SomaticQC
import "../NeoVax/neovax_dragen_modules/LegoPlotter.wdl" as LegoPlotter

struct Runtime {
  Int max_retries
  Int preemptible
  Int mem
  Int cpu
  Int boot_disk_size
  Int initial_disk_size
  String docker
}

workflow WESOutofSpecAnalysis {
    input {
    String pair_name
	String case_name
    File maf_file
    File somatic_vcf
    File mut_categs
    
    File tumor_bam
    File tumor_bam_index
    String cov_string = "exome"

    File normal_bam
    File normal_bam_index

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File dbSNP_vcf
    File dbSNP_vcf_index

	File bait_intervals
	File target_intervals

    String validation_stringency
    String? validation_errors_to_ignore

    # deTiN parameters
    File exac_pickle
    Float TiN_prior = 0.2
    Float mutation_prior = 0.05

    # ContEst parameters
    File snp6_bed
    File hapmap_vcf

    # Orientation bias Q-value
    File pre_adapter_detail_metrics

    File? gatk_jar_override

    File cnv_pon
    File common_snp_list

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 50
    Int boot_disk_size = 15
    Int mem = 16
    Int cpu = 1

    String cga_tools_docker = "us.gcr.io/tag-public/neovax-tag-cga:v1"
    String gatk_docker = "us.gcr.io/tag-public/neovax-tag-gatk:v1"
  }
  Runtime standard_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": cga_tools_docker,
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }
  Runtime gatk_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": gatk_docker,
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }
    
    call LegoPlotter.CreateLegoPlot as LegoPlotter {
    	input:
      pair_name = pair_name,
      maf_file = maf_file,
      mut_categs = mut_categs,
      cov_string = cov_string,
      runtime_params = standard_runtime
    }

    call AllelicCNV.AllelicCNV as AllelicCNV {
    	input:
    	pair_name = pair_name,
    	tumor_bam = tumor_bam,
    	tumor_bam_index = tumor_bam_index,
    	normal_bam = normal_bam,
    	normal_bam_index = normal_bam_index,
    	ref_fasta = ref_fasta,
    	ref_fasta_index = ref_fasta_index,
    	ref_dict = ref_dict,
    	cnv_pon = cnv_pon,
    	common_snp_list = common_snp_list
    }
    call PairBamQC.PicardWesMetrics as TumorPicardMetrics {
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
      runtime_params = gatk_runtime
  }

  call Absolute.RunAbsolute as Absolute {
    input:
      analysis_id = pair_name,
      seg_file = AllelicCNV.alleliccapseg_tsv,
      maf_file = maf_file,
      skew = AllelicCNV.alleliccapseg_skew,
      runtime_params = standard_runtime
  }

  call SomaticQC.SomaticQC as SomaticQC {
    input:
    	pair_name = pair_name,
    	tumor_bam = tumor_bam,
    	tumor_bam_index = tumor_bam_index,
    	normal_bam = normal_bam,
    	normal_bam_index = normal_bam_index,
    	ref_fasta = ref_fasta,
    	ref_fasta_index = ref_fasta_index,
    	ref_dict = ref_dict,
    	target_intervals = target_intervals,
    	input_vcf = somatic_vcf,
    	input_seg = AllelicCNV.alleliccapseg_tsv,
    	tumor_hets = AllelicCNV.tumor_hets,
    	normal_hets = AllelicCNV.normal_hets,
    	exac_pickle = exac_pickle,
    	TiN_prior = TiN_prior,
    	mutation_prior = mutation_prior,
    	snp6_bed = snp6_bed,
    	hapmap_vcf = hapmap_vcf,
    	pre_adapter_detail_metrics = TumorPicardMetrics.pre_adapter_detail_metrics
  }

  output {
    ## LegoPlotter Outputs
    Array[File] lego_plotter_ais = LegoPlotter.ais
    Array[File] lego_plotter_pngs = LegoPlotter.pngs
    Array[File] lego_plotter_figs = LegoPlotter.figs
    Array[File] lego_plotter_pss = LegoPlotter.pss
    File mut_legos_html = LegoPlotter.mut_legos_html

    ## Absolute Outputs
    File absolute_highres_plot = Absolute.absolute_highres_plot

    ## SomaticQC Outputs
    File contamination_data = SomaticQC.contamination_data
    File contest_base_report = SomaticQC.contest_base_report
    Float frac_contam = SomaticQC.frac_contam
    String frac_contam_CI = SomaticQC.frac_contam_CI
    Float TiN = SomaticQC.TiN
    String TiN_CI = SomaticQC.TiN_CI
    Float oxog_q_value = SomaticQC.oxog_q_value
    Float ffpe_q_value = SomaticQC.ffpe_q_value
  }
}