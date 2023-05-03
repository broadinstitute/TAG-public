version 1.0

workflow SelectSampleFromCallSet {
    input {
        String bcftools_docker
        File variant_callset
        String sample_name
    }

    # Select vcf for specific sample
    call SelectSampleRemoveComplexSV {
        input:
            bcftools_docker = bcftools_docker,
            vcf = variant_callset,
            sample_name = sample_name
    }


    # Outputs that will be retained when execution is complete
    output {
        File selected_vcf = SelectSampleRemoveComplexSV.output_vcf
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "SelectSampleFromCallSet.wdl is design to extract specific sample from a large callset vcf"
    }
}
    # Select sample from a vcf and remove complex SVs and INV
    task SelectSampleRemoveComplexSV {

            input {
                String bcftools_docker
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
                # Select sample using bcftools
                bcftools view -s ~{sample_name} -O v -o ~{sample_name}.vcf ~{vcf}

                # Remove Complex SV from the sample vcf because wittyer can't process CPX variants
                # Remove INV from the sample vcf because wittyer's exception
                # Remove reference allele

                bcftools view -e 'SVTYPE="INV" | SVTYPE="CPX" | GT="0/0"' ~{sample_name}.vcf -o ~{sample_name}_filtered.vcf

                >>>
            runtime {
                docker: bcftools_docker
                bootDiskSizeGb: 12
                memory: mem_size + " GB"
                disks: "local-disk " + disk_size + " HDD"
                preemptible: 2
            }
            output {
                File output_vcf = "~{sample_name}_filtered.vcf"
            }

    }