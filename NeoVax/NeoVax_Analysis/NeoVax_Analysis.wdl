version 1.0

import "../neovax_dragen_modules/reformat_hla.wdl" as reformat_hla
import "../neovax_dragen_modules/Optitype.wdl" as Optitype
import "../neovax_dragen_modules/Polysolver.wdl" as Polysolver
import "../neovax_dragen_modules/AlignTranscripts.wdl" as AlignTranscripts
import "../neovax_dragen_modules/RNASeQC.wdl" as RNASeQC
import "../neovax_dragen_modules/RSEM.wdl" as RSEM
import "../neovax_dragen_modules/LegoPlotter.wdl" as LegoPlotter
import "../neovax_dragen_modules/PairBamQC.wdl" as PairBamQC
import "../neovax_dragen_modules/AllelicCNV.wdl" as AllelicCNV
import "../neovax_dragen_modules/SomaticQC.wdl" as SomaticQC
import "../neovax_dragen_modules/Absolute.wdl" as Absolute
import "../neovax_dragen_modules/MutationValidator.wdl" as MutationValidator
import "../neovax_dragen_modules/NeoAntigen.wdl" as neoantigen
import "../neovax_dragen_modules/HLAthena.wdl" as HLAthena
import "../neovax_dragen_modules/CollectMNP.wdl" as CollectMNP



workflow NeoVax_Analysis {
	input {
		## Terra Table and CollectMNP Inputs ##
		String patient_id
		String normal_id
		String normal_crsp_id
		File normal_bam
		File normal_bam_index
		String tumor_id
		String tumor_crsp_id
		File tumor_bam
		File tumor_bam_index
		File? tumor_rna_bam
		File hla_file
		File somatic_maf
		File somatic_vcf
        File somatic_vcf_index

		## Polysolver Inputs ##
		String reference_version = "hg19"
		String fastq_format = "STDFQ"
		Boolean include_freq = true
		Boolean insert_calc = false

		## RNA Workflow Inputs ##
		File star_index
		Boolean extract_umi = false
		Boolean umi_aware_mode = false
		File genes_gtf
		File? rna_genome_bam_input
		File? rna_genome_bam_index_input
		File? transcriptome_bam_input
		File rsem_reference
		Boolean estimate_rspd = true
		Boolean is_paired_end = true
		Int max_frag_len = 1000
		String strandedness = "none"

		## WES QC Workflow Inputs ##
		File ref_fasta
		File ref_fasta_index
		File ref_fasta_dict
		File target_intervals
		File bait_intervals
		File dbSNP_vcf
		File dbSNP_vcf_index
		File haplotype_database
		File exac_pickle
		File hapmap_vcf
		File snp6_bed
		File mut_categs
		String validation_errors_to_ignore = "INVALID_TAG_NM"
		String cov_string = "exome"
		Boolean add_chr_prefix_to_rna = false
		String maf_type = "wex"
                ## WES CNV Workflow Inputs ##
		File cnv_pon
		File common_snp_list

		## Neoantigen Workflow Inputs ##
		File duplicate_genes
		File gene_id_name_mapping
		File hg19DBTarBall
		File wild_type_peptides_tar_gz
		Boolean aggregate_pep = false
		String assign_by_ranks_or_scores = "ranks"
		String assign_colors = "#8c0a82, #ab0c9f, #048edb, #00a0fa, #5F36A0, #ff6000, darkorange, #e05400, grey"
		String assign_threshold = "0.1"
		Boolean exists_ctex = true
		Boolean exists_expr = true
		String expr_col_name = "TPM"
		Boolean logtransform_expr = true
		String peptide_col_name = "pep"
		# run NetMHCPan if executable tarball is supplied
		File? netmhc_tar_gz

	}

    call CollectMNP.CollectMNP as CollectMNP {
		input:
			output_prefix = patient_id,
			input_maf = somatic_maf,
			input_vcf = somatic_vcf,
			input_vcf_index = somatic_vcf_index,
			tumor_bam = tumor_bam,
			tumor_bam_index = tumor_bam_index,
			normal_bam = normal_bam,
			normal_bam_index = normal_bam_index,
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			ref_dict = ref_fasta_dict
    }
    File maf_file = select_first([CollectMNP.output_maf, somatic_maf])
    File vcf_file = select_first([CollectMNP.mnp_merged_vcf, somatic_vcf])
    File vcf_index_file = select_first([CollectMNP.mnp_merged_vcf_index, somatic_vcf_index])

    call reformat_hla.reformat_hla as reformat_hla {
		input:
			hla_file = hla_file,
			patient_id = patient_id
	}

	call Optitype.Optitype as Optitype {
		input:
			normal_bam = normal_bam,
			normal_bam_index = normal_bam_index,
			prefix = patient_id
	}

	call Polysolver.Polysolver as Polysolver {
		input:
			normal_bam = normal_bam,
			normal_bam_index = normal_bam_index,
			prefix = patient_id,
			reference_version = reference_version,
			fastq_format = fastq_format,
			include_freq = include_freq,
			insert_calc = insert_calc
	}
	if (defined(tumor_rna_bam) && !defined(rna_genome_bam_input)){
		call AlignTranscripts.AlignTranscripts as AlignTranscripts {
			input:
				prefix = tumor_id,
				input_bam = tumor_rna_bam,
				star_index = star_index,
				extract_umi = extract_umi,
				umi_aware_mode = umi_aware_mode
		}
	
		call RNASeQC.RNASeQC as RNASeQC {
			input:
				genes_gtf = genes_gtf,
				prefix = patient_id,
				input_bam = AlignTranscripts.output_bam
		}
	}
	call RSEM.RSEM as RSEM {
		input:
			prefix = patient_id,
			rsem_reference = rsem_reference,
			transcriptome_bam = select_first([AlignTranscripts.transcriptome_bam, transcriptome_bam_input]),
			estimate_rspd = estimate_rspd,
			is_paired_end = is_paired_end,
			max_frag_len = max_frag_len,
			strandedness = strandedness
	}

	call LegoPlotter.LegoPlotter as LegoPlotter {
		input:
			maf_file = maf_file,
			mut_categs = mut_categs,
			pair_name = patient_id,
			cov_string = cov_string
	}

	call PairBamQC.PairBamQC as PairBamQC {
		input:
			pair_name = patient_id,
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			ref_dict = ref_fasta_dict,
			bait_intervals = bait_intervals,
			target_intervals = target_intervals,
			case_name = tumor_id,
			control_name = normal_id,
			dbSNP_vcf = dbSNP_vcf,
			dbSNP_vcf_index = dbSNP_vcf_index,
			normal_bam = normal_bam,
			normal_bam_index = normal_bam_index,
			tumor_bam = tumor_bam,
			tumor_bam_index = tumor_bam_index,
			haplotype_database = haplotype_database,
			validation_errors_to_ignore = validation_errors_to_ignore
	}

	call AllelicCNV.AllelicCNV as AllelicCNV {
		input:
			cnv_pon = cnv_pon,
			common_snp_list = common_snp_list,
			pair_name = patient_id,
			ref_dict = ref_fasta_dict,
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			tumor_bam = tumor_bam,
			tumor_bam_index = tumor_bam_index,
			normal_bam = normal_bam,
			normal_bam_index = normal_bam_index
	}

	call SomaticQC.SomaticQC as SomaticQC {
		input:
			exac_pickle = exac_pickle,
			hapmap_vcf = hapmap_vcf,
			input_seg = AllelicCNV.alleliccapseg_tsv,
			input_vcf = vcf_file,
			normal_bam = normal_bam,
			normal_bam_index = normal_bam_index,
			normal_hets = AllelicCNV.normal_hets,
			pair_name = patient_id,
			pre_adapter_detail_metrics = PairBamQC.tumor_pre_adapter_detail_metrics,
			ref_dict = ref_fasta_dict,
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			snp6_bed = snp6_bed,
			tumor_bam = tumor_bam,
			tumor_bam_index = tumor_bam_index,
			tumor_hets = AllelicCNV.tumor_hets
	}

	call MutationValidator.MutationValidator as MutationValidator {
		input:
			maf_file = maf_file,
			normal_bam = normal_bam,
			normal_bam_index = normal_bam_index,
			pair_name = patient_id,
			tumor_bam = tumor_bam,
			tumor_bam_index = tumor_bam_index,
			tumor_rna_bam = select_first([AlignTranscripts.output_bam, rna_genome_bam_input]),
			tumor_rna_bam_index = select_first([AlignTranscripts.output_bam_index, rna_genome_bam_index_input]),
			add_chr_prefix_to_rna = add_chr_prefix_to_rna,
			maf_type = maf_type
	}

	call Absolute.Absolute as Absolute {
		input:
			maf_file = maf_file,
			pair_name = patient_id,
			seg_file = AllelicCNV.alleliccapseg_tsv,
			skew = AllelicCNV.alleliccapseg_skew
	}

	call neoantigen.Neoantigen_Pipeline as neoantigen {
		input:
			caseName = tumor_crsp_id,
			duplicate_genes = duplicate_genes,
			gene_id_name_mapping = gene_id_name_mapping,
			geneExpr = RSEM.genes_expr,
			hg19DBTarBall = hg19DBTarBall,
			id = patient_id,
			mafFile = maf_file
	}
	
	call HLAthena.HLAthena as HLAthena {
		input:
			aggregate_pep = aggregate_pep,
			alleles_file = reformat_hla.reformatted_hla,
			assign_by_ranks_or_scores = assign_by_ranks_or_scores,
			assign_colors = assign_colors,
			assign_threshold = assign_threshold,
			exists_ctex = exists_ctex,
			exists_expr = exists_expr,
			expr_col_name = expr_col_name,
			logtransform_expr = logtransform_expr,
			peptide_col_name = peptide_col_name,
			peptide_list = neoantigen.snp_indel_neoantigens,
            netmhc_tar_gz = netmhc_tar_gz,
            analysis_name = patient_id,
			wild_type_peptides_tar_gz = wild_type_peptides_tar_gz,
			second_peptide_list = neoantigen.neoorf_neoantigens
	}

  output {
    ## CollectMNP Output Files ##
    File output_mnp_maf = maf_file
    Int mnp_counts = CollectMNP.mnp_counts
    File output_mnp_maflite = CollectMNP.output_mnp_maflite
    File output_mnp_vcf = CollectMNP.output_mnp_vcf
    File output_mnp_tool_log = CollectMNP.output_mnp_tool_log
    File? mnp_merged_vcf = CollectMNP.mnp_merged_vcf
    File? mnp_merged_vcf_index = CollectMNP.mnp_merged_vcf_index

    ## HLA Output Files ##
    File reformatted_hla = reformat_hla.reformatted_hla
    File optitype_coverage_plot = Optitype.coverage_plot
    File optitype_hla_output = Optitype.output_hla
    File polysolver_hla_output = Polysolver.output_hla

    ## RNA Pipeline Output Files ##
    File? rna_duplicate_metrics = AlignTranscripts.duplicate_metrics
    File? rna_genome_bam = AlignTranscripts.output_bam
    File? rna_genome_bam_index = AlignTranscripts.output_bam_index
    File? transcriptome_bam = AlignTranscripts.transcriptome_bam
    File? exon_counts = RNASeQC.exon_counts
    File? fragment_size = RNASeQC.fragment_size
    File? gene_counts = RNASeQC.gene_counts
    File? gene_duplicates = RNASeQC.gene_duplicates
    File? gene_tpm = RNASeQC.gene_tpm
    File? rna_metrics = RNASeQC.metrics
    File genes_expr = RSEM.genes_expr
    File isoforms_expr = RSEM.isoforms_expr

    ## WES QC Task Outputs ##
    File mut_legos_html = LegoPlotter.mut_legos_html
    Array[File] lego_plotter_pngs = LegoPlotter.lego_plotter_pngs
    File? cross_check_fingerprint_metrics = PairBamQC.cross_check_fingprt_metrics
    String? cross_check_fingerprint_min_lod_lanes = PairBamQC.cross_check_fingprt_min_lod_lanes
    Float? cross_check_fingerprint_min_lod_value = PairBamQC.cross_check_fingprt_min_lod_value
    File? cross_check_fingerprint_report = PairBamQC.cross_check_fingprt_report
    File? normal_picard_metrics = PairBamQC.normal_picard_metrics
    File? normal_picard_plots = PairBamQC.normal_picard_plots
    File tumor_picard_metrics = PairBamQC.tumor_picard_metrics
    File tumor_picard_plots = PairBamQC.tumor_picard_plots
    File tumor_pre_adapter_detail_metrics = PairBamQC.tumor_pre_adapter_detail_metrics
    Float ffpe_q_value = SomaticQC.ffpe_q_value
    Float oxog_q_value = SomaticQC.oxog_q_value
    Float frac_contam = SomaticQC.frac_contam
    String frac_contam_CI = SomaticQC.frac_contam_CI
    Float TiN = SomaticQC.TiN
    String TiN_CI = SomaticQC.TiN_CI
    File validated_maf = MutationValidator.validated_maf

    ## WES CNV Task Outputs ##
    File alleliccapseg_plot = AllelicCNV.alleliccapseg_plot
    Float alleliccapseg_skew = AllelicCNV.alleliccapseg_skew
    File alleliccapseg_tsv = AllelicCNV.alleliccapseg_tsv
    File copy_ratio_plot = AllelicCNV.copy_ratio_plot
    File copy_ratio_plot_lim4 = AllelicCNV.copy_ratio_plot_lim4
    File acnv_seg = AllelicCNV.output_seg
    Float pre_tn_mad = AllelicCNV.pre_tn_mad
    Float post_tn_mad = AllelicCNV.post_tn_mad
    File absolute_highres_plot = Absolute.absolute_highres_plot
    File absolute_rdata = Absolute.absolute_rdata

    ## Neoantigen Task Outputs ##
    File outINDEL_for_ms_predictions = neoantigen.outINDEL_for_ms_predictions
    File outSNP_for_ms_predictions = neoantigen.outSNP_for_ms_predictions
    File combinedFullpeptideFile = neoantigen.combinedFullpeptideFile
    File neoorf_neoantigens = neoantigen.neoorf_neoantigens
    File neoorfFile = neoantigen.neoorfFile
    File snp_indel_neoantigens = neoantigen.snp_indel_neoantigens
	File nonsense_vars_file = neoantigen.nonsense_vars_file
    File hlathena_all_binders = HLAthena.all_binders
    File hlathena_allele_assignment_counts = HLAthena.allele_assignment_counts
    File hlathena_plots = HLAthena.allele_assignment_plots
    File hlathena_mut_binders = HLAthena.mut_binders
  }
}
