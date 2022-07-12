version 1.0

workflow CollectSamErrorMetrics {
	input {
		File input_bam
        File vcf_file
        File vcf_index_file
        File ref_fasta
        File ref_fai
        File ref_dict
        File intervals
		String output_name
		String picard_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
	}

	call CollectSamErrorMetricsTask{
		input:
			input_bam = input_bam,
            picard_docker=picard_docker,
            vcf_file=vcf_file,
            vcf_index_file=vcf_index_file,
            ref_fasta=ref_fasta,
            ref_fai=ref_fai,
            ref_dict=ref_dict,
            intervals=intervals,
			output_name = output_name
	}
	output{
		Array[File] output_file = CollectSamErrorMetricsTask.output_metrics
	}
}

task CollectSamErrorMetricsTask {
	input {
		String output_name
		String picard_docker
		File input_bam
        File ref_fasta
        File ref_fai
        File ref_dict
        File vcf_file
        File vcf_index_file
        File intervals
		Int disk_gb = 400
		Int memory_gb = 16
	}
	command {
    java -Xms2000m -jar /usr/gitc/picard.jar CollectSamErrorMetrics I=~{input_bam} \
      O=~{output_name} R=~{ref_fasta} V=~{vcf_file} L=~{intervals} \
      ERROR_METRICS=null \
      ERROR_METRICS="ERROR:INSERT_LENGTH" \
      ERROR_METRICS="ERROR:CYCLE" \
      ERROR_METRICS="ERROR:GC_CONTENT" \
      ERROR_METRICS="ERROR" 

	}
	output {
		Array[File] output_metrics = glob("~{output_name}*")
	}
	runtime {
		docker: picard_docker
		disks: "local-disk " + disk_gb + " HDD"
		memory: memory_gb + "GB"
		cpu: "1"
	}
}