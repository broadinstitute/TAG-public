version 1.0

workflow Benchmark_CNV_Caller {
    input {
        String gatk_docker
        File variant_callset
        String sample_name
    }
    call SelectVariant {
        input:
        gatk_docker = gatk_docker,
        vcf = variant_callset,
        sample_name = sample_name
    }
    output {
        File outfile = SelectVariant.output_vcf
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