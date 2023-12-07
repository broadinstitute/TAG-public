version 1.0

import "../neovax_dragen_modules/runTPMCutoff.wdl" as run_tpm_cutoff

workflow NeoVax_PackageData {
	input{
	String pair_id
	String patient_id
	}
	call absolute_extract {
		input: pair_id = pair_id
	}

	call pre_variant_review {
		input: patient_id = pair_id
	}

	call run_tpm_cutoff.findTPMCutoff as findTPMCutoff {
		input: patient_id = patient_id
	}

	call extract_metric_task as extract_depth_tumor
	call extract_plots_task as extract_plots_tumor {
		input: output_basename = pair_id
	}
	call extract_metric_task as extract_depth_normal
	call extract_plots_task as extract_plots_normal {
		input: output_basename = pair_id
	}

	output{
		## TPM cutoff plot outputs
		File cutoff_plot_pdf = findTPMCutoff.cutoff_plot_pdf
		File cutoff_plot_png = findTPMCutoff.cutoff_plot_png
		Float cutoff_tpm = findTPMCutoff.first_cutoff_tpm
		Int cutoff_genes = findTPMCutoff.first_cutoff_genes

		## Absolute extracted outputs
		File absolute_annotated_maf_capture = absolute_extract.absolute_annotated_maf_capture
		File absolute_seg_file = absolute_extract.absolute_seg_file
		File absolute_segdat_file = absolute_extract.absolute_segdat_file
		File absolute_table = absolute_extract.absolute_table
		String purity = absolute_extract.purity
		String ploidy = absolute_extract.ploidy

		## Automatic variant pre-review
		File coding_only_reviewed_maf = pre_variant_review.coding_only_reviewed_maf
		File coding_only_reviewed_maf_tsv = pre_variant_review.coding_only_reviewed_maf_tsv
		File missing_genes_list = pre_variant_review.missing_genes_list
		File variant_review_log = pre_variant_review.variant_review_log

		## Extracted HybSel outputs
		Float normal_pct_above_100x = extract_depth_normal.output_metric
		Float tumor_pct_above_100x = extract_depth_tumor.output_metric
		File normal_insert_size_histogram = extract_plots_normal.insert_size_histogram
		File normal_quality_by_cycle_plot = extract_plots_normal.quality_by_cycle_plot
		File normal_quality_distribution_plot = extract_plots_normal.quality_distribution_plot
		File tumor_insert_size_histogram = extract_plots_tumor.insert_size_histogram
		File tumor_quality_by_cycle_plot = extract_plots_tumor.quality_by_cycle_plot
		File tumor_quality_distribution_plot = extract_plots_tumor.quality_distribution_plot
	}
}
task absolute_extract{

	input{
	File absolute_summary_data_gatk_acnv
	String absolute_called_solution_capture
	String pair_id
	String analyst_id
	Int memory = 7
	Int diskGB = 50
	String docker = "us.gcr.io/tag-team-160914/neovax-tag-absolute-extract:v1"
	}

	command <<<
	abs_solution_number=`python -c 'print int("~{absolute_called_solution_capture}".split("-")[-1])'`
	Rscript /usr/local/bin/ABSOLUTE_extract_cli_start.R --solution_num $abs_solution_number --analyst_id ~{analyst_id} --rdata_modes_fn ~{absolute_summary_data_gatk_acnv} --sample_name ~{pair_id} --results_dir . --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/
	cp reviewed/SEG_MAF/~{pair_id}_ABS_MAF.txt .
	cp reviewed/SEG_MAF/~{pair_id}.segtab.txt .
	cp reviewed/samples/~{pair_id}.ABSOLUTE.~{analyst_id}.called.RData .
	cp reviewed/~{pair_id}.~{analyst_id}.ABSOLUTE.table.txt .
	cut -f4 ~{pair_id}.~{analyst_id}.ABSOLUTE.table.txt|tail -n1 > purity
	cut -f5 ~{pair_id}.~{analyst_id}.ABSOLUTE.table.txt|tail -n1 > ploidy
	>>>

	runtime {
	docker: docker
	preemptible: 1
	disks: "local-disk ~{diskGB} HDD"
	memory: "~{memory}GB"
	maxRetries: 1
	}

	output {
	File absolute_annotated_maf_capture = "~{pair_id}_ABS_MAF.txt"
	File absolute_seg_file = "~{pair_id}.segtab.txt"
	File absolute_segdat_file = "~{pair_id}.ABSOLUTE.~{analyst_id}.called.RData"
	File absolute_table = "~{pair_id}.~{analyst_id}.ABSOLUTE.table.txt"
	String purity = read_string("purity")
	String ploidy = read_string("ploidy")
	}
}

task pre_variant_review{
	input{
	String patient_id
	File mutation_validator_maf
	File neoantigen_snp_indel
	File nonsense_vars
	File? rsem_genes
	File? rsem_isoforms
	File? rnaseqc_genes
	File? gene_name_map
	Int diskGB = 50
	Float memory = 3.5
	String docker = "us.gcr.io/tag-team-160914/neovax-tag-variant-review:v4"
	}

	String output_basename = basename(mutation_validator_maf,".maf")

	command {
		python /pre_variant_review.py -p ~{patient_id} -m ~{mutation_validator_maf} -n ~{neoantigen_snp_indel} -s ~{nonsense_vars} \
		~{'-r ' + rsem_genes} ~{rsem_isoforms} ~{rnaseqc_genes} ~{'-g ' + gene_name_map}
	}

	runtime {
		docker: docker
		preemptible: 1
		disks: "local-disk ~{diskGB} HDD"
		memory: "~{memory}GB"
		maxRetries: 1
	}

	output {
	File coding_only_reviewed_maf = "~{patient_id}/results/~{output_basename}.coding_only.reviewed.maf"
	File coding_only_reviewed_maf_tsv = "~{patient_id}/results/~{output_basename}.coding_only.reviewed.maf.tsv"
	File missing_genes_list = "~{patient_id}/results/~{output_basename}.unmatched_genes.txt"
	File variant_review_log = "process_results.log"
	}
}

task extract_metric_task {

	input{
	File picard_metrics_tar
	String extract_col = "PCT_TARGET_BASES_100X"

	String docker = "us.gcr.io/tag-team-160914/neovax-tag-variant-review:v4"
	}

	command <<<
	export HYBSEL=`tar -ztf ~{picard_metrics_tar} | grep hybrid_selection_metrics`
	tar -xzvf ~{picard_metrics_tar} $HYBSEL
	python <<CODE
import os
with open(os.environ["HYBSEL"]) as metrics:
	col_index = -1
	for line in metrics:
		if line.startswith("BAIT_SET"):
			metrics_list = line.rstrip("\n").split("\t")
			col_index = metrics_list.index("~{extract_col}")
			values = next(metrics).rstrip("\n").split("\t")
			metric_value = values[col_index]
			break
outfile = open("hybsel_metric.txt", "w")
outfile.write(str(metric_value))
outfile.close()
CODE
	>>>

	runtime {
	docker: docker
	disk: "50 GB"
	}

	output {
	Float output_metric = read_float("hybsel_metric.txt")
	}
}

task extract_plots_task {
	input{
	File picard_plots_tar
	String output_basename
	String docker = "us.gcr.io/tag-team-160914/neovax-tag-variant-review:v4"
	}
	command <<<
	INSERT_PLOT=`tar -ztf ~{picard_plots_tar} | grep insert_size_histogram.pdf`
	QBC_PLOT=`tar -ztf ~{picard_plots_tar} | grep quality_by_cycle.pdf`
	QUAL_PLOT=`tar -ztf ~{picard_plots_tar} | grep quality_distribution.pdf`

	tar -xzvf ~{picard_plots_tar} $INSERT_PLOT $QBC_PLOT $QUAL_PLOT

	mv $INSERT_PLOT ~{output_basename}.insert_size_histogram.pdf
	mv $QBC_PLOT ~{output_basename}.quality_by_cycle.pdf
	mv $QUAL_PLOT ~{output_basename}.quality_distribution.pdf
	>>>

	runtime {
	docker: docker
	disk: "50 GB"
	}

	output {
	File insert_size_histogram = "~{output_basename}.insert_size_histogram.pdf"
	File quality_by_cycle_plot = "~{output_basename}.quality_by_cycle.pdf"
	File quality_distribution_plot = "~{output_basename}.quality_distribution.pdf"
	}
}