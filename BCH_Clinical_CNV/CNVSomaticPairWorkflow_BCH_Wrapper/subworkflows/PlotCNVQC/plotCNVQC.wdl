version 1.0

workflow plotCNVQC {
	call plotQCTask

	output {
		File CNV_QC_plots = plotQCTask.CNV_QC_plots
	}
}

task plotQCTask {
	input {
		File? plot_qc_script_override
		File called_copy_ratio_segments_tumor
		File called_copy_ratio_segments_normal
		File somatic_maf
		File ref_fasta_index
		String output_basename
		Float denoised_MAD_value_normal
		Float denoised_MAD_value_tumor

		String docker = "us.gcr.io/tag-public/cnv-qc-plots:1.0.0"
		Float memGb = 3.5
		Int disk_size = 50
		Int preemptible = 2
	}

	command {
		if [ -f "~{plot_qc_script_override}" ]; then
			PLOT_QC_SCRIPT="~{plot_qc_script_override}"
		else
			PLOT_QC_SCRIPT="/scripts/CNV_QC_Plot.R"
		fi
		
		Rscript $PLOT_QC_SCRIPT \
		--tumor_seg ~{called_copy_ratio_segments_tumor} \
		--normal_seg ~{called_copy_ratio_segments_normal} \
		--tumor_maf ~{somatic_maf} \
		--fasta_idx ~{ref_fasta_index} \
		--tumor_mad ~{denoised_MAD_value_tumor} \
		--normal_mad ~{denoised_MAD_value_normal} \
		--output ~{output_basename}
	}

	output {
		File CNV_QC_plots = "~{output_basename}_CNV_QC.pdf"
	}

	runtime {
			docker: docker
			bootDiskSizeGb: 12
			memory: memGb + " GB"
			disks: "local-disk " + disk_size + " HDD"
			preemptible: preemptible
	}
}