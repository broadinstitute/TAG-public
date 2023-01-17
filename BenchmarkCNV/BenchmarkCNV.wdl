version 1.0

workflow Benchmark_CNV_Caller {
    input {
        String gatk_docker
        File variant_callset
        String truth_sample_name
        String query_sample_name
        String wittyer_docker
        File eval_cnv_vcf
        File wittyer_cnv_config
        String wittyer_cnv_evaluation_mode
        File eval_sv_vcf
        File wittyer_sv_config
        String wittyer_sv_evaluation_mode
    }

    # Select vcf for specific sample
    call SelectVariant {
        input:
            gatk_docker = gatk_docker,
            vcf = variant_callset,
            truth_sample_name = truth_sample_name
    }

    # benchmark cnv.vcf using witty.er
    call BenchmarkCNV {
        input:
            wittyer_docker = wittyer_docker,
            truth_vcf = SelectVariant.output_vcf,
            truth_sample_name = truth_sample_name,
            query_sample_name = query_sample_name,
            eval_cnv_vcf = eval_cnv_vcf,
            cnv_config_file = wittyer_cnv_config,
            cnv_evaluation_mode = wittyer_cnv_evaluation_mode,
            eval_sv_vcf = eval_sv_vcf,
            sv_config_file = wittyer_sv_config,
            sv_evaluation_mode = wittyer_sv_evaluation_mode,
    }

    # Outputs that will be retained when execution is complete
    output {
        File truth_vcf = SelectVariant.output_vcf
        File cnv_wittyer_stats = BenchmarkCNV.cnv_wittyer_stats
        File cnv_wittyer_annotated_vcf = BenchmarkCNV.cnv_wittyer_annotated_vcf
        File cnv_wittyer_annotated_vcf_index = BenchmarkCNV.cnv_wittyer_annotated_vcf_index
        File sv_wittyer_stats = BenchmarkCNV.sv_wittyer_stats
        File sv_wittyer_annotated_vcf = BenchmarkCNV.sv_wittyer_annotated_vcf
        File sv_wittyer_annotated_vcf_index = BenchmarkCNV.sv_wittyer_annotated_vcf_index
    }
}

    # Task 1: Select sample vcf from a large callset (e.g. 1KGP)
    task SelectVariant {

    input {
        String gatk_docker
        File vcf
        String truth_sample_name
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
        --sample-name ~{truth_sample_name} \
        --exclude-non-variants \
        --remove-unused-alternates \
        -O ~{truth_sample_name}.vcf

        # Remove Complex SV from the sample vcf because wittyer can't process CPX variants
        # Remove INV from the sample vcf because wittyer's exception
        cat ~{truth_sample_name}.vcf | grep -v '<CPX>' | grep -v '<INV>' > ~{truth_sample_name}_filtered.vcf
    >>>
    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
    }
    output {
        File output_vcf = "~{truth_sample_name}_filtered.vcf"
    }

}

    # Task 2: Benchmark the large variant vcf against
    task BenchmarkCNV {

    input {
        String wittyer_docker
        File truth_vcf
        File eval_cnv_vcf
        File cnv_config_file
        String cnv_evaluation_mode
        File eval_sv_vcf
        File sv_config_file
        String sv_evaluation_mode
        String truth_sample_name
        String query_sample_name
        Int? mem
        Int? disk_space
        # If mem and disk size were not specified, use 4GB and 100 GB as default
        Int mem_size = select_first([mem, 4])
        Int disk_size = select_first([disk_space,100])
}
    command <<<
        set -e
        # Run Benchmarking tool wittyer on dragen generated cnv.vcf
        /opt/Wittyer/Wittyer -i ~{eval_cnv_vcf} \
        -t ~{truth_vcf} \
        -em ~{cnv_evaluation_mode} \
        --configFile ~{cnv_config_file} \
        -o ~{truth_sample_name}_cnv_wittyer_output

        # Run Benchmarking tool wittyer on dragen generated sv.vcf
        /opt/Wittyer/Wittyer -i ~{eval_cnv_vcf} \
        -t ~{truth_vcf} \
        -em ~{sv_evaluation_mode} \
        --configFile ~{sv_config_file} \
        -o ~{truth_sample_name}_sv_wittyer_output

    >>>
    runtime {
        docker: wittyer_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
    }
    output {
        File cnv_wittyer_stats = "~{truth_sample_name}_cnv_wittyer_output/Wittyer.Stats.json"
        File cnv_wittyer_annotated_vcf = "~{truth_sample_name}_cnv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz"
        File cnv_wittyer_annotated_vcf_index = "~{truth_sample_name}_cnv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz.tbi"
        File sv_wittyer_stats = "~{truth_sample_name}_sv_wittyer_output/Wittyer.Stats.json"
        File sv_wittyer_annotated_vcf = "~{truth_sample_name}_sv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz"
        File sv_wittyer_annotated_vcf_index = "~{truth_sample_name}_sv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz.tbi"
    }
}