workflow CopyUmiFromReadName {
	call CopyUmiTask
}

task CopyUmiTask {
    String? bloodbiopsydocker = "us.gcr.io/tag-team-160914/liquidbiopsy:0.0.4.5"
    String base_name
    String? fgbio_override
    File bam_file
    File bam_index
    Boolean? remove_umi_from_read_name = true
    
    Int? preemptible = 2
    Int? maxRetries = 1
    Int? disk_pad
    Int disk_size = ceil(size(bam_file, "GB") * 5) + select_first([disk_pad,0])
    Float? extra_mem
    Float mem = 25 + select_first([extra_mem, 0])
    Int? cpu = 4
    Int compute_mem = ceil(mem) * 1000 - 500
    
    command {
        export FGBIO_LOCAL_JAR=${default="/usr/fgbio-2.0.2.jar" fgbio_override}

        ln -vs ${bam_file} ${base_name}_input.bam
        ln -vs ${bam_index} ${base_name}_input.bai

        java -Xmx${compute_mem}m -jar $FGBIO_LOCAL_JAR \
        CopyUmiFromReadName \
        -i ${base_name}_input.bam \
        -o ${base_name}.bam \
        --remove-umi ${remove_umi_from_read_name}
    }
    
    output {
        File umi_extracted_bam = "${base_name}.bam"
        File umi_extracted_bam_index = "${base_name}.bai"
    }
    
    runtime {
      docker: select_first([bloodbiopsydocker])
      disks: "local-disk " + disk_size + " HDD, /cromwell_root/tmp 500 HDD"
      memory: mem + " GB"
      maxRetries: select_first([maxRetries])
      preemptible: select_first([preemptible])
      cpu: select_first([cpu])
  	}

}