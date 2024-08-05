version 1.0

workflow CramToBamWorkflow {
    input {
        File referenceFasta
        File referenceFai
        File referenceDict
        File cram_file
        String sample_id
        String docker_image
        String gatk_docker
        String output_mode
    }
    call CramToBam {
        input: 
            referenceFasta = referenceFasta,
            referenceFai = referenceFai,
            referenceDict = referenceDict,
            cram_file = cram_file,
            sample_id = sample_id,
            docker_image = docker_image
    }

    call ValidateBam {
        input: 
            bam_file = CramToBam.output_bam,
            sample_id = sample_id,
            output_mode = output_mode,
            gatk_docker = gatk_docker
    }

    output {
        File bam_file = CramToBam.output_bam
        File bam_index_file = CramToBam.output_bam_index
        File validation_output = ValidateBam.validation_output
    }
    
    meta {
        author: "TAG Team, Patrick Schlaeger(updated script writer), Yueyao Gao(bugfixer)"
        email: "tag@broadinstitute.org, pschlaeg@broadinstitute.org, gaoyueya@broadinstitute.org"
        description: "(Bug Fix and version update in 2024) This workflow converts an inputted Cram file into a Bam file with the use of the Cram file's reference genome dictionary/fasta/fasta_index. It then validates the the newly created Bam file."
    }
}


    task CramToBam {
        input {
        File referenceFasta
        File referenceFai
        File referenceDict
        File cram_file
        String sample_id
        String docker_image = "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
        Int? mem_size = 32
        Int? disk_size = 500
        Int? preemptible_attempts = 2
        }
        command <<<
            set -e
            set -o pipefail

            ln -vs ~{referenceFasta} reference.fasta
            ln -vs ~{referenceFai} reference.fasta.fai
            ln -vs ~{referenceDict} reference.dict

            samtools view -h -T reference.fasta ~{cram_file} |
            samtools view -b -o ~{sample_id}.bam -
            samtools index -b ~{sample_id}.bam
            mv ~{sample_id}.bam.bai ~{sample_id}.bai
        >>>
        output {
            File output_bam = "~{sample_id}.bam"
            File output_bam_index = "~{sample_id}.bai"
        }
        runtime {
            docker: docker_image
            memory: mem_size + " GB"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: preemtible_attempts
        }
}


    task ValidateBam {
        input {
        File bam_file
        String sample_id
        String gatk_docker = "broadinstitute/gatk:4.6.0.0"
        String output_mode = "SUMMARY"
        Int? mem_size = 32
        Int? disk_size = 500
        Int? preemptible_attempts = 2
        }
        Int machine_mem_mb = select_first([mem_size, 7]) * 1000
        Int command_mem_mb = machine_mem_mb - 1000
        command <<<
            gatk --java-options "-Xmx~{command_mem_mb}m" ValidateSamFile \
                --INPUT ~{bam_file} \
                --OUTPUT "~{sample_id}_validation_output.txt" \
                --MODE ~{output_mode}
        >>>
        output {
            File validation_output = "~{sample_id}_validation_output.txt"
        }
        runtime {
            docker: gatk_docker
            memory: mem_size + " GB"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: preemtible_attempts
        }

}
