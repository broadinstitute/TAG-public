workflow CallRevertSam {
call RevertSam
}

task RevertSam {
    File input_bam
    File? input_bam_index
    String base_name
    String? sort_order
    String? additional_args

    String? gatk_path = "/gatk/gatk"
    File? ref_fasta
    File? ref_fasta_index
    File? ref_fasta_dict

    String? docker_override
    String docker = select_first([docker_override, "us.gcr.io/broad-gatk/gatk:4.5.0.0"])
    Int? preemptible_count = 2
    Int? maxRetries = 1
    Int? threads = 4
    Float? mem = 15.0
    Int? disk_buffer = 50

    command <<<
        ${gatk_path} \
        	RevertSam \
        	--INPUT ${input_bam} \
        	--OUTPUT ${base_name}.reverted.bam \
            ${'--REFERENCE_SEQUENCE ' + ref_fasta} \
            --VALIDATION_STRINGENCY SILENT \
            ${additional_args} \
        	${"--SORT_ORDER "+ sort_order}
    >>>

    output {
        File output_bam = "${base_name}.reverted.bam"
    }

    runtime {
        docker: docker
        disks: "local-disk " + sub(((size(input_bam,"GB")+1)*2+disk_buffer),"\\..*","") + " HDD"
        memory: mem + " GB"
        cpu: threads
        maxRetries: maxRetries
        preemptible: preemptible_count
    }
}