task maf_aggregator_task
	{

	String memGB
	String diskGB
	String pSetID
    String cpus
	Array[File] mafs

	command <<<

	#increase verbosity
	set -x

	#run catters
	python /usr/local/bin/tsvcat.py ${sep=" " mafs} > ${pSetID}.aggregated.maf

	>>>

	output {
		File aggregated_maf="${pSetID}.aggregated.maf"
		}

	runtime {
		docker: "broadinstitute/broadmutationcalling_filtering_beta@sha256:d2df6d9d705e90d3ee926b72a3edec5034dd5bdd64c0cf1cabd9fc8462702c79"
		memory: "${memGB} GB"
        cpu: "${cpus}"
		disks: "local-disk ${diskGB} HDD"
		}

	}

workflow maf_aggregator_workflow {

	call maf_aggregator_task

}