version 1.0
workflow AnnotateBed{
    input {
        File annotate_script
        File gencode_annotation
        File bed_to_annotate
        String output_prefix
        File gene_bed
        Boolean generate_gene_base_count
        File? gene_list
        File gene_base_count_script
        Int diskGB = 50
    }
    call GenerateAnnotation {
        input:
        script = annotate_script,
        gencode_annotation = gencode_annotation, 
        bed_to_annotate = bed_to_annotate,
        output_prefix = output_prefix
    }
     scatter (i in [0]) {
        if ((generate_gene_base_count) && defined(gene_base_count_script) && defined(gene_bed)) {
            call CountGeneBases {
                input:
                script = gene_base_count_script,
                bed_to_annotate = bed_to_annotate,
                gene_bed = gene_bed,
                gene_list = gene_list,
                grouped_by_gene = GenerateAnnotation.grouped_by_gene,
                output_prefix = output_prefix
            }
        }
    }
        output {
        File ungrouped_annotation = "~{output_prefix}.ungrouped_annotation.txt"
        File annotation_per_interval = "~{output_prefix}.grouped_by_interval.annotated.txt"
        Int intergetic_base_count = GenerateAnnotation.intergetic_base_count
        Int coding_base_count = GenerateAnnotation.coding_base_count
        File grouped_by_gene = "~{output_prefix}.grouped_by_gene.txt" 
        Int genes_involved = CountGeneBases.genes_involved
        File? gene_base_count = "~{output_prefix}.grouped_by_gene.annotated.txt"

    }
}
task GenerateAnnotation {
    input {
        File script
        File gencode_annotation
        File bed_to_annotate
        String output_prefix
        Int maxRetries = 1
        Int preemptible = 3
        Int diskGB = 50
    }
    command {
        python ~{script} --annotation ~{gencode_annotation} --bed ~{bed_to_annotate} --output_prefix ~{output_prefix}
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/annotate_bed"
        preemptible: preemptible
        disks: "local-disk ~{diskGB} HDD"
        memory: "4GB"
        maxRetries: maxRetries
    }
    output {
        File ungrouped_annotation = "~{output_prefix}.ungrouped_annotation.txt"
        File annotation_per_interval = "~{output_prefix}.grouped_by_interval.annotated.txt"
        File grouped_by_gene = "~{output_prefix}.grouped_by_gene.txt"
        Int intergetic_base_count = read_int(stdout())
        Int coding_base_count = read_int(stdout())
    }
}

task CountGeneBases {
    input{
        File script
        File gene_bed
        File bed_to_annotate
        File? gene_list
        File grouped_by_gene = grouped_by_gene
        Int maxRetries = 1
        Int preemptible = 3
        Int diskGB = 50
        String output_prefix
    }
    command {
            python ~{script} --gene_bed ~{gene_bed} --bed ~{bed_to_annotate} --grouped_by_gene ~{grouped_by_gene} ~{'--gene-list ' + gene_list} --output_prefix ~{output_prefix}
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/annotate_bed"
        preemptible: preemptible
        disks: "local-disk ~{diskGB} HDD"
        memory: "4GB"
        maxRetries: maxRetries
    }
    output {
        Int genes_involved = read_int(stdout())
        File gene_base_count = "~{output_prefix}.grouped_by_gene.annotated.txt"
    }
}

