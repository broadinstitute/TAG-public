workflow MergeSvabaAnnotations{
	call mergeAnnotationsTask
}

task mergeAnnotationsTask{
	String output_basename
	File merge_annotations_script
	File svaba_vcf
	File snpsift_table
	File annot_table
	File contig_seqs

	Float? extra_mem
	Int? extra_disk
	Int? maxRetries
	Int? preemptible
	Int? cpu
    String? docker = "us.gcr.io/tag-team-160914/tag-tools:1.0.0"

	
	command{
		Rscript ${merge_annotations_script} ${svaba_vcf} ${snpsift_table} \
		${annot_table} ${contig_seqs} ${output_basename}
	}

	runtime{
        docker: select_first([docker])
        memory: 3.5 + select_first([extra_mem, 0]) + " GB"
        disks: "local-disk " + (50 + select_first([extra_disk, 0])) + " HDD"
        preemptible: select_first([preemptible, 3])
        maxRetries: select_first([maxRetries, 1])
        cpu: select_first([cpu, 2])
	}

	output{
		File final_annotation = "${output_basename}.annotated.sv.txt"
	}
}