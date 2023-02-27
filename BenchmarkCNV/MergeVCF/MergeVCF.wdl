version 1.0

workflow MergeVCF {
    input {
        String bcftools_docker
        File vcf1
        File vcf2
        String sample_name
    }

    # Merge Two Input VCF Files
    call mergeVCF {
        input:
            bcftools_docker = bcftools_docker,
            vcf1 = vcf1,
            vcf2 = vcf2,
            sample_name = sample_name
    }


    # Outputs that will be retained when execution is complete
    output {
        File selected_vcf = mergeVCF.output_vcf
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "MergeVCF.wdl is design to merge two input vcfs"
    }
}
task mergeVCF {

        input {
            String bcftools_docker
            File vcf1
            File vcf2
            String sample_name
            Int? mem
            Int? disk_space
            # If mem and disk size were not specified, use 4GB and 100 GB as default
            Int mem_size = select_first([mem, 4])
            Int disk_size = select_first([disk_space,100])

        }
        command <<<
            set -e

            echo "Input VCF files:"
            echo "vcf1: ~{vcf1}"
            echo "vcf2: ~{vcf2}"

            bcftools concat ~{vcf1} ~{vcf2} -o ~{sample_name}_merged.vcf

            # Print the output file name
            echo "Merged VCF file: ~{sample_name}_merged.vcf"
            >>>
        runtime {
            docker: bcftools_docker
            bootDiskSizeGb: 12
            memory: mem_size + " GB"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: 2
        }
        output {
            File output_vcf = "~{sample_name}_merged.vcf"
        }

}