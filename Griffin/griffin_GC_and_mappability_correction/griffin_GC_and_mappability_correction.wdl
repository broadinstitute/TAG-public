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
        String griffin_docker
        File mappable_regions
        File reference_genome
        File genome_GC_frequency
        Int GC_bias_size_min
        Int GC_bias_size_max
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
            mappable_regions = mappable_regions,
            reference_genome = reference_genome,
            chrom_sizes = chrom_sizes,
            map_q = map_quality,
            GC_bias_size_min = GC_bias_size_min,
            GC_bias_size_max = GC_bias_size_max
    }

    call GC_bias {
        input:
            griffin_docker = griffin_docker,
            sample_name = sample_name,
            GC_counts_file = GC_counts.GC_counts_file,
            mappable_name = basename(mappable_regions, ".bed"),
            genome_GC_frequency = genome_GC_frequency,
            GC_bias_size_min = GC_bias_size_min,
            GC_bias_size_max = GC_bias_size_max
    }

    output {
        File? mappability_bias_file = mappability_bias.mappability_bias_file
        File? mappability_plot = mappability_bias.mappability_plot
        File? mappability_plot2 = mappability_bias.mappability_plot2
        File GC_counts_file = GC_counts.GC_counts_file
        File GC_bias_file = GC_bias.GC_bias_file
        File GC_plots = GC_bias.GC_plots
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "griffin_GC_and_mappability_correction.wdl is the first task of Griffin"
    }
}

task mappability_bias {
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

        # Index the input bam file
        conda run --no-capture-output -n griffin_env samtools index ~{bam_file}

        # Make temporary directory
        mkdir -p results/mappability_bias/temp_~{sample_name}/

        # Create a directory for mappability plots
        mkdir -p results/mappability_plots/

        # Run griffin_mappability_correction
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_mappability_correction.py \
        --bam_file ~{bam_file} \
        --bam_file_name ~{sample_name} \
        --output results/mappability_bias/~{sample_name}.mappability_bias.txt \
        --output_plot results/mappability_plots/~{sample_name}.mappability_bias.pdf \
        --mappability ~{mappability_bw} \
        --exclude_paths ~{encode_exclude} ~{centromeres} ~{gaps} ~{patches} ~{alternative_haplotypes} \
        --chrom_sizes ~{chrom_sizes} \
        --map_quality ~{map_q} \
        --CPU ~{cpu_num} \
        --tmp_dir results/mappability_bias/temp_~{sample_name}/

        # Remove temporary directory
        rm -r results/mappability_bias/temp_~{sample_name}/
    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
}
    output {
        File mappability_bias_file = "results/mappability_bias/~{sample_name}.mappability_bias.txt"
        File mappability_plot = "results/mappability_plots/~{sample_name}.mappability_bias.pdf"
        File mappability_plot2 = "results/mappability_plots/~{sample_name}.mappability_bias.read_coverage_distribution.pdf"
    }
}

task GC_counts {
    input {
        String griffin_docker
        File bam_file
        String sample_name
        File mappable_regions
        File reference_genome
        File chrom_sizes
        Int map_q
        Int GC_bias_size_min
        Int GC_bias_size_max
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

        # Create an output directory
        mkdir results

        # Index the input bam file
        conda run --no-capture-output -n griffin_env samtools index ~{bam_file}

        # Run griffin_GC_counts
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_GC_counts.py \
        --bam_file ~{bam_file} \
        --bam_file_name ~{sample_name} \
        --mappable_regions_path ~{mappable_regions} \
        --ref_seq ~{reference_genome} \
        --chrom_sizes ~{chrom_sizes} \
        --out_dir results \
        --map_q ~{map_q} \
        --size_range ~{GC_bias_size_min} ~{GC_bias_size_max} \
        --CPU ~{cpu_num} \

    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
}
    output {
        File GC_counts_file = "results/GC_counts/~{sample_name}.GC_counts.txt"
    }
}

task GC_bias {
    input {
        String griffin_docker
        String sample_name
        String mappable_name
        File GC_counts_file
        File genome_GC_frequency
        Int GC_bias_size_min
        Int GC_bias_size_max
    }
    command <<<
        set -e

        # Uncompress the genome_GC_frequency.tar.gz into a direcotry named genome_GC_frequency
        tar -xvf ~{genome_GC_frequency} --directory ./

        # Create an output directory
        mkdir results

        # Run griffin_GC_bias
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_GC_bias.py \
        --bam_file_name ~{sample_name} \
        --mappable_name ~{mappable_name} \
        --genome_GC_frequency ./genome_GC_frequency \
        --out_dir results \
        --size_range ~{GC_bias_size_min} ~{GC_bias_size_max} \

    >>>
    runtime {
        docker: griffin_docker
        preemptible: 2
        bootDiskSizeGb: 12
        cpu: 2
        memory: "4 GB"
        disks: "local-disk " + 25 + " HDD"
}
    output {
        File GC_bias_file = "results/GC_bias/~{sample_name}.GC_bias.txt"
        File GC_plots = "results/GC_plots/~{sample_name}.GC_bias.summary.pdf"
    }
}