version 1.0

    workflow MergeBamAlignment {
        input{
            File aligned_bam
            File unmapped_bam
            File reference_fasta
            String basename
            Boolean unmap_contaminant_reads
            String? additional_args
            String? docker_override
        }
        call MergeBamAlignment {
            input:
                aligned_bam = aligned_bam,
                unmapped_bam = unmapped_bam,
                basename = basename,
                reference_fasta = reference_fasta,
                unmap_contaminant_reads = unmap_contaminant_reads
        }
        output{
            File merged_bam = MergeBamAlignment.merged_bam
        }
    }

    task MergeBamAlignment {
        input{
            File aligned_bam
            File unmapped_bam
            File reference_fasta
            String basename
            Boolean unmap_contaminant_reads
            String docker = select_first([docker_override, "us.gcr.io/broad-gatk/gatk:4.4.0.0"])
            String? docker_override
            String? additional_args
            Int? preemptible_count = 3
            Float? mem = 8
            Int? disk_size = 500
        }
        command <<<
            set -e
            mkdir -p /cromwell_root/reference

            cp ~{reference_fasta} /cromwell_root/reference/reference.fasta
            /gatk/gatk CreateSequenceDictionary \
                -R /cromwell_root/reference/reference.fasta \
                -O /cromwell_root/reference/reference.dict

            /gatk/gatk MergeBamAlignment \
                -ALIGNED_BAM ~{aligned_bam} \
                -UNMAPPED_BAM ~{unmapped_bam} \
                -O ~{basename}.merged.bam \
                -REFERENCE_SEQUENCE /cromwell_root/reference/reference.fasta \
                -UNMAP_CONTAM ~{unmap_contaminant_reads} \
                --CREATE_INDEX true \
                ${additional_args}
        >>>
        output{
            File merged_bam = "~{basename}.merged.bam"
            File merged_bam = "~{basename}.merged.bai"
        }
        runtime{
            docker: docker
            disks: "local-disk " + disk_size + " HDD"
            memory: mem + " GB"
            preemptible: preemptible_count
        }
    }