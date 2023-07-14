workflow AnnotateBreaKmer{
    call annotate_breakmer_sv

    output{
        File annotated_sv = annotate_breakmer_sv.annotated_sv
    }
}

task annotate_breakmer_sv{
    File annotation_bed
    File breakmer_sv
    String output_basename

    String? docker_override 
    String docker = select_first([docker_override, "us.gcr.io/tag-public/annotate-breakmer:v1"])
    Float? extra_mem
    Float mem = 3.5 + select_first([extra_mem, 0])
    Int? extra_disk 
    Int disk = 50 + select_first([extra_disk, 0])
    Int? preemptible = 3
    Int? maxRetries = 1
    
    command {
        python /annotate_breakmer.py -a ${annotation_bed} -b ${breakmer_sv} -o ${output_basename}
    }

    output {
        File annotated_sv = "${output_basename}.annotated.sv.tsv"
    }

    runtime {
        docker: docker
        memory: mem+" GB"
        disks: "local-disk " + disk + " HDD"
        preemptible: select_first([preemptible, 3])
        maxRetries: select_first([maxRetries, 1])
    }
}