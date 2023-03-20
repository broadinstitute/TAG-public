version 1.0

workflow nucleosome_profiling{
    input {
        String griffin_docker
        File bam_file
        File GC_bias_file
        String sample_name
        File? mappability_bias
        File reference_genome
        File mappability_bw
        File chrom_sizes
        File sites_yaml
        String chrom_column
        String position_column
        String strand_column
        Array[String] chroms
        Array[Int] norm_window
        Array[Int] size_range
        Int map_quality
        Int number_of_sites = 0
        String sort_by = 'none'
        String ascending = 'none'
    }

    call calc_cov {
        input:
            griffin_docker = griffin_docker,
            bam_file = bam_file,
            GC_bias_file = GC_bias_file,
            sample_name = sample_name,
            mappability_bias = mappability_bias,
            reference_genome = reference_genome,
            mappability_bw = mappability_bw,
            chrom_sizes = chrom_sizes,
            sites_yaml = sites_yaml,
            chrom_column = chrom_column,
            position_column = position_column,
            strand_column = strand_column,
            chroms = chroms,
            norm_window = norm_window,
            size_range = size_range,
            map_quality = map_quality,
            number_of_sites = number_of_sites,
            sort_by = sort_by,
            ascending = ascending
    }

    output {
        File uncorrected_bw = calc_cov.uncorrected_bw
        File GC_corrected_bw = calc_cov.GC_corrected_bw
        File tmp_pybedtools = calc_cov.tmp_pybedtools
    }

    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "griffin_nucleosome_profiling.wdl is the second task of Griffin"
    }
}

task calc_cov {
    input {
        String griffin_docker
        File bam_file
        File GC_bias_file
        String sample_name
        File? mappability_bias
        File reference_genome
        File mappability_bw
        File chrom_sizes
        File sites_yaml
        String chrom_column
        String position_column
        String strand_column
        Array[String] chroms
        Array[Int] norm_window
        Array[Int] size_range
        Int map_quality
        Int number_of_sites
        String sort_by
        String ascending
        Int? cpu
        Int? mem
        Int? disk_space
        # If cpu, mem, and disk size were not specified, use 8 cores, 10GB, and 100 GB as default
        Int cpu_num = select_first([cpu, 8])
        Int mem_size = select_first([mem, 10])
        Int disk_size = select_first([disk_space,100])
        }
    command <<<
        set -e
        # Index the input bam file
        conda run --no-capture-output -n griffin_env samtools index ~{bam_file}

        # Make temporary directory
        mkdir -p results/calc_cov/temp/

        # Run griffin_coverage_script if mappability_bias path was specified
        # when mappability_bias path was specified
        # mappability_correction is True
        if [[ -f "~{mappability_bias}" ]]; then
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_coverage.py \
        --sample_name ~{sample_name} \
        --bam ~{bam_file} \
        --GC_bias ~{GC_bias_file} \
        --mappability_bias ~{mappability_bias} \
        --mappability_correction True \
        --tmp_dir results/calc_cov/temp \
        --reference_genome ~{reference_genome} \
        --mappability_bw ~{mappability_bw} \
        --chrom_sizes_path ~{chrom_sizes} \
        --sites_yaml ~{sites_yaml} \
        --griffin_scripts /BaseImage/Griffin/scripts/ \
        --chrom_column ~{chrom_column} \
        --position_column ~{position_column} \
        --strand_column ~{strand_column} \
        --chroms ~{sep=" " chroms} \
        --norm_window ~{sep=" " norm_window} \
        --size_range ~{sep=" " size_range} \
        --map_quality ~{map_quality} \
        --number_of_sites ~{number_of_sites} \
        --sort_by ~{sort_by} \
        --ascending ~{ascending} \
        --CPU ~{cpu_num}

        else
        # Run griffin_coverage_script if mappability_bias path was not specified
        # when mappability_bias path was not specified
        # mappability_correction is False
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_coverage.py \
        --sample_name ~{sample_name} \
        --bam ~{bam_file} \
        --GC_bias ~{GC_bias_file} \
        --mappability_bias none \
        --mappability_correction False \
        --tmp_dir results/calc_cov/temp \
        --reference_genome ~{reference_genome} \
        --mappability_bw ~{mappability_bw} \
        --chrom_sizes_path ~{chrom_sizes} \
        --sites_yaml ~{sites_yaml} \
        --griffin_scripts /BaseImage/Griffin/scripts/ \
        --chrom_column ~{chrom_column} \
        --position_column ~{position_column} \
        --strand_column ~{strand_column} \
        --chroms ~{sep=" " chroms} \
        --norm_window ~{sep=" " norm_window} \
        --size_range ~{sep=" " size_range} \
        --map_quality ~{map_quality} \
        --number_of_sites ~{number_of_sites} \
        --sort_by ~{sort_by} \
        --ascending ~{ascending} \
        --CPU ~{cpu_num}

        # tar zip the tmp_pybedtools
        tar -zcvf results/calc_cov/temp/tmp_pybedtools.tar.gz results/calc_cov/temp/tmp_pybedtools
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
        File uncorrected_bw = "results/calc_cov/temp/tmp_bigWig/~{sample_name}.uncorrected.bw"
        File GC_corrected_bw = "results/calc_cov/temp/tmp_bigWig/~{sample_name}.GC_corrected.bw"
        File tmp_pybedtools = "results/calc_cov/temp/tmp_pybedtools.tar.gz"
        }
}
