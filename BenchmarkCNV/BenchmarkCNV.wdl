version 1.0

workflow Benchmark_CNV_Caller {
    input {
        String gatk_docker
        File variant_callset
        String sample_name
        File eval_vcf
        File wittyer_config
        String wittyer_docker
        String wittyer_evaluation_mode
    }

    # Select vcf for specific sample
    call SelectVariant {
        input:
        gatk_docker = gatk_docker,
        vcf = variant_callset,
        sample_name = sample_name
    }

    # benchmark cnv.vcf using witty.er
    call BenchmarkCNV {
        input:
        wittyer_docker = wittyer_docker,
        truth_vcf = SelectVariant.output_vcf,
        eval_vcf = eval_vcf,
        sample_name = sample_name,
        config_file = wittyer_config,
        evaluation_mode = wittyer_evaluation_mode
    }

    # Outputs that will be retained when execution is complete
    output {
        File truth_vcf = SelectVariant.output_vcf
        File wittyer_stats = BenchmarkCNV.wittyer_stats
    }
}

    # Task 1: Select sample vcf from a large callset (e.g. 1KGP)
    task SelectVariant {

    input {
        String gatk_docker
        File vcf
        String sample_name
        Int? mem
        Int? disk_space
        # If mem and disk size were not specified, use 4GB and 100 GB as default
        Int mem_size = select_first([mem, 4])
        Int disk_size = select_first([disk_space,100])

    }
    command <<<
        set -e
        export PATH="/gatk:$PATH"
        gatk --java-options "-Xmx4g" SelectVariants \
        -V ~{vcf} \
        --sample-name ~{sample_name} \
        --exclude-non-variants \
        --remove-unused-alternates \
        -O ~{sample_name}.vcf

        # Remove Complex SV from the sample vcf because wittyer can't process CPX variants
        # Remove INV from the sample vcf because wittyer's exception
        cat ~{sample_name}.vcf | grep -v '<CPX>' | grep -v '<INV>' > ~{sample_name}_filtered.vcf
    >>>
    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
    }
    output {
        File output_vcf = "~{sample_name}_filtered.vcf"
    }

}

    # Task 2: Benchmark the large variant vcf against
    task BenchmarkCNV {

    input {
        String wittyer_docker
        File truth_vcf
        File eval_vcf
        File config_file
        String evaluation_mode
        String sample_name
        Int? mem
        Int? disk_space
        # If mem and disk size were not specified, use 4GB and 100 GB as default
        Int mem_size = select_first([mem, 4])
        Int disk_size = select_first([disk_space,100])
}
    command <<<
    set -e
    # Run Benchmarking tool wittyer
    /opt/Wittyer/Wittyer -i ~{eval_vcf} \
    -t ~{truth_vcf} \
    -em ~{evaluation_mode} \
    --configFile ~{config_file} \
    -o ~{sample_name}_wittyer_output
    >>>
    runtime {
        docker: wittyer_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
    }
    output {
        File wittyer_stats = "~{sample_name}_wittyer_output/Wittyer.Stats.json"
    }
}