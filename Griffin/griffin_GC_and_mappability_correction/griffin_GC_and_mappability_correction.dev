version 1.0

workflow gc_and_mappability_correction{
    input {
        Boolean mappability_correction
        File bam_file
        String sample_name
        File mappability_bw
        File encode_exclude
        File centromeres
        File gaps
        File patches
        File alternative_haplotypes
        File chrom_sizes
        Int map_quality
        File mappable_regions
        File reference_genome
        Int GC_bias_size_min
        Int GC_bias_size_max
        File genome_GC_frequency
        String griffin_docker
    }
    if (mappability_correction){
        call mappability_bias {
            input:
                griffin_docker = griffin_docker,
                bam_file = bam_file,
                sample_name = sample_name,
                mappability_bw = mappability_bw,
                encode_exclude = encode_exclude,
                centromeres = centromeres,
                gaps = gaps,
                patches = patches,
                alternative_haplotypes = alternative_haplotypes,
                chrom_sizes = chrom_sizes,
                map_q = map_quality
    }
}
    call GC_counts {
        input:
            griffin_docker = griffin_docker,
            bam_file = bam_file,
            sample_name = sample_name,
            mappable_regions_path = mappable_regions,
            ref_seq = reference_genome,
            chrom_sizes = chrom_sizes,
            map_q = map_quality,
            GC_bias_size_min = GC_bias_size_min,
            GC_bias_size_max = GC_bias_size_max
}
    call GC_bias {
        input:
            griffin_docker = griffin_docker,
            GC_counts_file = GC_counts.GC_counts_file,
            mappable_regions_path = mappable_regions,
            genome_GC_frequency = genome_GC_frequency,
            GC_bias_size_min = GC_bias_size_min,
            GC_bias_size_max = GC_bias_size_max
}
    call make_samples_yaml{
        input:
            GC_bias_file = GC_bias.GC_bias_file,
            mappability_bias = mappability_bias.mappability_bias
    }
    output {
        File GC_counts_file = GC_counts.GC_counts_file
        File GC_bias_file = GC_bias.GC_bias_file
        File GC_bias_summary_plot = GC_bias.GC_bias_summary_plot
        File GC_bias_plot = GC_bias.GC_bias_plot
        File GC_bias_key_length_plot = GC_bias.GC_bias_key_length_plot
        File? mappability_bias = mappability_bias.mappability_bias
        File? mappability_plot = mappability_bias.mappability_plot
        File? mappability_plot2 = mappability_bias.mappability_plot2
        File samples_gc_yaml = make_samples_yaml.samples_gc_yaml

    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "griffin_GC_and_mappability_correction.wdl is the first task of Griffin"
    }
}
task mappability_bias{
    input {
        String griffin_docker
        File bam_file
        String sample_name
        File mappability_bw
        File encode_exclude
        File centromeres
        File gaps
        File patches
        File alternative_haplotypes
        File chrom_sizes
        Int map_q
        Int? cpu
        Int? mem
        Int? disk_space
        # If cpu, mem, and disk size were not specified, use 8 cores, 8GB, and 100 GB as default
        Int cpu_num = select_first([cpu, 8])
        Int mem_size = select_first([mem, 8])
        Int disk_size = select_first([disk_space,100])
}
    command <<<
        set -e

        # Run griffin_mappability_correction
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_mappability_correction.py \
        --bam_file ~{bam_file} \
        --bam_file_name ~{sample_name} \
        --output mappability_bias/~{sample_name}.mappability_bias.txt \
        --output_plot mappability_plots/~{sample_name}.mappability_bias.pdf \
        --mappability ~{mappability_bw} \
		--exclude_paths ~{encode_exclude} ~{centromeres} ~{gaps} ~{patches} ~{alternative_haplotypes} \
		--chrom_sizes ~{chrom_sizes} \
		--map_quality ~{map_q} \
		--CPU ~{cpu_num} \
		--tmp_dir tmp_~{sample_name}
    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
}
    output {
        File mappability_bias = "mappability_bias/~{sample_name}.mappability_bias.txt"
        File mappability_plot = "mappability_plots/~{sample_name}.mappability_bias.pdf"
        File mappability_plot2 = "mappability_plots/{samples}.mappability_bias.read_coverage_distribution.pdf"

}
}