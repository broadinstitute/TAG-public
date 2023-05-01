version 1.0

workflow gc_correction{
    input {
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
        Array[Int] GC_bias_size_range
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
            GC_bias_size_range = GC_bias_size_range
    }

    call GC_bias {
        input:
            griffin_docker = griffin_docker,
            sample_name = sample_name,
            GC_counts_file = GC_counts.GC_counts_file,
            mappable_name = basename(mappable_regions, ".bed"),
            genome_GC_frequency = genome_GC_frequency,
            GC_bias_size_range = GC_bias_size_range
    }

    output {
        File GC_counts_file = GC_counts.GC_counts_file
        File GC_bias_file = GC_bias.GC_bias_file
        File GC_plots = GC_bias.GC_plots
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "griffin_GC_correction.wdl is the first task of Griffin"
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
        Array[Int] GC_bias_size_range
        # If cpu, mem, and disk size were not specified, use 8 cores, 8GB, and 100 GB as default
        Int cpu_num = 8
        Int mem_size = 8
        Int disk_size = 100
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
        --size_range ~{sep=" " GC_bias_size_range} \
        --CPU ~{cpu_num} \

    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 3
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
        Array[Int] GC_bias_size_range
    }
    command <<<
        set -e

        # Uncompress the genome_GC_frequency.tar.gz into a direcotry named genome_GC_frequency
        tar -xvf ~{genome_GC_frequency} --directory ./

        # Copy the GC_count_file that were generated from GC_Count task
        mkdir -p results/GC_counts/
        cp ~{GC_counts_file} results/GC_counts/

        # Run griffin_GC_bias
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_GC_bias.py \
        --bam_file_name ~{sample_name} \
        --mappable_name ~{mappable_name} \
        --genome_GC_frequency ./genome_GC_frequency \
        --out_dir results \
        --size_range ~{sep=" " GC_bias_size_range}

    >>>
    runtime {
        docker: griffin_docker
        preemptible: 2
        bootDiskSizeGb: 12
        cpu: 2
        memory: "4 GB"
        disks: "local-disk " + 25 + " HDD"
        preemptible: 2
        maxRetries: 3
}
    output {
        File GC_bias_file = "results/GC_bias/~{sample_name}.GC_bias.txt"
        File GC_plots = "results/GC_plots/~{sample_name}.GC_bias.summary.pdf"
    }
}